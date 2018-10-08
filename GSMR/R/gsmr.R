# gsmr: A tool for GSMR and HEIDI analysis
# The gsmr package Perform Generalized Summary-data-based Mendelian Randomization analysis (GSMR)
#  and HEterogeneity In Dependent Instruments analysis to remove pleiotropic outliers (HEIDI-outlier).
# @author Zhihong Zhu <z.zhu1@uq.edu.au>
# @author Zhili Zheng <zhili.zheng@uq.edu.au>
# @author Futao Zhang <f.zhang5@uq.edu.au>
# @author Jian Yang <jian.yang@uq.edu.au>

eps = 1e-6
library("survey");

# ************************************************** #
#         check data is missing or not               #
# ************************************************** #
check_element <- function(vals, argue) {
    vals = vals[which(is.finite(vals))]
    if(length(vals)==0) {
        stop("None of the ", argue, " is found. Please check.")
    }
}

# ************************************************** #
#          convert variable nameto string            #
# ************************************************** #
var2string <- function(argue) {
    deparse(substitute(argue))
}

# ************************************************** #
#         check parameter is defined or not          #
# ************************************************** #
check_variable <- function(argue) {
    if(!exists(argue)) {
        stop(argue, " is not defined.")
    }
}

# ************************************************** #
#       variance of bXY                              #
# ************************************************** #
var_bXY <- function(bzx, bzx_se, bzy, bzy_se) {
   varbXY = bzx_se^2*bzy^2/bzx^4 + bzy_se^2/bzx^2
   return(varbXY)
}


# ************************************************** #
#       variance-covariane matrix of bXY             #
# ************************************************** #
cov_bXY <- function(bxy_i, bxy_j, bzx, bzx_se, bzy, bzy_se, ldrho) {
    zscoreZX = bzx / bzx_se
    nsnp = dim(ldrho)[1]
    covbXY = diag(nsnp)
    if(nsnp>1) {
        zszxij = zscoreZX%*%t(zscoreZX)
        bxyij = bxy_i*bxy_j
        sezyij = bzy_se%*%t(bzy_se)
        bzxij = bzx%*%t(bzx)
        covbXY = ldrho*sezyij/bzxij + ldrho*bxyij/zszxij
    }
    return(covbXY)
}


# ************************************************** #
#       bXY by GSMR method                           #
# ************************************************** #

bxy_gsmr = function(bzx, bzx_se, bzy, bzy_se, ldrho) {
    bXY = bzy/bzx
    covbXY = cov_bXY(median(bXY), median(bXY), bzx, bzx_se, bzy, bzy_se, ldrho)
    diag(covbXY) = diag(covbXY) + eps
    # Eigen decomposition
    resbuf = eigen(covbXY, symmetric=TRUE)
    eval = as.numeric(resbuf$values)
    evec = resbuf$vectors
    if(min(eval) < eps) {
        stop("The variance-covariance matrix for bxy is not invertible!");
    }
    covbXY_inv = evec%*%diag(1/eval)%*%t(evec)
    vec_1 = rep(1, length(bzx))
    num_1_v_1 = as.numeric(solve(t(vec_1)%*%covbXY_inv%*%vec_1))
    vec_1_v = as.numeric(t(vec_1)%*%covbXY_inv)
    bXY_GLS = as.numeric(num_1_v_1*vec_1_v%*%bXY)
    varbXY_GLS = num_1_v_1
    chisqbXY_GLS = bXY_GLS^2/varbXY_GLS
    pbXY_GLS = pchisq(chisqbXY_GLS, 1, lower.tail=F)
    return(list(bxy=bXY_GLS, bxy_se=sqrt(varbXY_GLS), bxy_pval=pbXY_GLS, vec_1t_v=vec_1_v))
}

# ************************************************** #
#       standard HEIDI-outlier                       #
# ************************************************** #
std_heidi_pvalue <- function(bxy_hat, bxy_hat_se, bzx, bzx_se, bzy, bzy_se, ldrho, vec_1t_v, snp_index) {
    m = length(bzx)
    # Estimate var(bxy_i)
    var_bxy = var_bXY(bzx, bzx_se, bzy, bzy_se);
    # Estimate the cov(bxy_i, bgsmr)
    bxy = bzy/bzx
    cov_bxy = cov_bXY(bxy, bxy_hat, bzx, bzx_se, bzy, bzy_se, ldrho)
    cov_bxy = cov_bxy[,snp_index];
    cov_bxy_bgsmr = sapply(1:m, function(x) bxy_hat_se^2*vec_1t_v%*%(cov_bxy[x,]))
    # Estimate p-value
    d = bxy - bxy_hat;
    var_d = var_bxy + bxy_hat_se^2 - 2*cov_bxy_bgsmr;
    pval_het = pchisq(d^2/var_d, 1, lower.tail=F)
    return(list(d=d, var_d=var_d, pval_het=pval_het))
}

std_heidi_outlier <- function(bzx, bzx_se, bzy, bzy_se, ldrho) {
    # Estimate bxy
    bxy = bzy/bzx
    bxy_q = quantile(bxy, probs=seq(0, 1, 0.1))
    snp_index = which(bxy >= bxy_q[2] & bxy <= bxy_q[10])
    resbuf = bxy_gsmr(bzx[snp_index], bzx_se[snp_index], bzy[snp_index], bzy_se[snp_index], ldrho[snp_index,snp_index])
    bxy_hat = resbuf$bxy; bxy_hat_se = resbuf$bxy_se; vec_1t_v = resbuf$vec_1t_v;

    std_het_pval = std_heidi_pvalue(bxy_hat, bxy_hat_se, bzx, bzx_se, bzy, bzy_se, ldrho, vec_1t_v, snp_index)$pval_het
    return(std_het_pval)
}

# ************************************************** #
#      global HEIDI-outlier                          #
# ************************************************** #
global_heidi_pvalue <- function(bxy_hat, bxy_hat_se, vec_1t_v, bzx, bzx_se, bzy, bzy_se, ldrho) {
    m = length(bzx)
    # Estimate Vd matrix
    cov_bxy = cov_bXY(bxy_hat, bxy_hat, bzx, bzx_se, bzy, bzy_se, ldrho);
    cov_bxy_bgsmr = sapply(1:m, function(x) bxy_hat_se^2*vec_1t_v%*%(cov_bxy[x,]))
    var_d = cov_bxy - cov_bxy_bgsmr + bxy_hat_se^2
    var_d = t(t(var_d) - cov_bxy_bgsmr)
    # Estimate p-value
    d = bzy/bzx - bxy_hat; chival = d^2/diag(var_d)
    # Eigen decomposition
    # correlation matrix of diff
    corr_d = diag(m)
    for( i in 1 : (m-1) ) {
        for( j in (i+1) : m ) {
            corr_d[i,j] = corr_d[j,i] =
                   var_d[i,j] / sqrt(var_d[i,i]*var_d[j,j])
        }
    }
    # estimate the p value
    lambda = eigen(corr_d, symmetric=TRUE, only.values=TRUE)$values
    t = length(lambda)
    pval_het = pchisqsum(sum(chival), df=rep(1,t), a=lambda, method="sadd", lower.tail=F)

    return(pval_het);
}

global_heidi_outlier <- function(bzx, bzx_se, bzy, bzy_se, ldrho, pval_thresh = 0.05) {
    # Set pval_het
    m = length(bzx)
    pval_het = 0;
    excl_index = c()
    # Back up
    kept_index = c(1:m);
    while(1) {
        # Update summary stats
        bzx2 = bzx[kept_index]; bzx2_se = bzx_se[kept_index];
        bzy2 = bzy[kept_index]; bzy2_se = bzy_se[kept_index]; ldrho2 = ldrho[kept_index,kept_index]
        # Estimate bxy by GSMR using the given SNPs
        resbuf = bxy_gsmr(bzx2, bzx2_se, bzy2, bzy2_se, ldrho2);
        bxy_hat = resbuf$bxy; bxy_hat_se = resbuf$bxy_se; vec_1t_v = resbuf$vec_1t_v
        pval_global = global_heidi_pvalue(bxy_hat, bxy_hat_se, vec_1t_v, bzx2, bzx2_se, bzy2, bzy2_se, ldrho2);
        if(pval_global >= pval_thresh) break;
        snpindxbuf = c(1:length(bzx2))
        pval_std = std_heidi_pvalue(bxy_hat, bxy_hat_se, bzx2, bzx2_se, bzy2, bzy2_se, ldrho2, vec_1t_v, snpindxbuf)$pval_het;
        excl_index = which.min(pval_std);
        kept_index = kept_index[-excl_index];
    }
    bxy_hat_pval = pchisq((bxy_hat/bxy_hat_se)^2, 1, lower.tail=F)
    raw_index = c(1:m)
    pleio_index = raw_index[-match(kept_index, raw_index)]
    return(list(bxy=bxy_hat, bxy_se=bxy_hat_se, bxy_pval=bxy_hat_pval, pval_het=pval_global, pleio_index = pleio_index))
}


# ************************************************** #
#              Test identical elements               #
# ************************************************** #
check_vec_elements_eq <- function(dim_vec){
    return(all(dim_vec == dim_vec[1]))
}


# ************************************************** #
#                     LD pruning                     #
# ************************************************** #
# LD pruning, removing pairs of SNPs that LD r2 > threshold
# return index: the index be used for further analysis
snp_ld_prune = function(ldrho, ld_r2_thresh) {
    # initialization
    nsnp = dim(ldrho)[1]
    diag(ldrho) = 0
    ldrho[upper.tri(ldrho)]=0
    include_id = c(1:nsnp)

    # save the index which have ld r^2 > threshold
    indx = which(ldrho^2>ld_r2_thresh, arr.ind=T)
    if(length(indx)==0) return(NULL);
    indx1 = as.numeric(indx[,2])
    indx2 = as.numeric(indx[,1])
    slct_indx = c(indx1, indx2)
    indxbuf = unique(sort(slct_indx))

    # count how many SNPs in high LD
    nproc = length(indxbuf)
    n_slct_snp = as.numeric()
    for( i in 1 : nproc) {
        n_slct_snp[i] = length(which(slct_indx==indxbuf[i]))
    }  

    # decide the index to remove
    nproc = length(indx1)
    if(nproc==0) return(NULL);
    for( i in 1 : nproc ) {
        n1 = n_slct_snp[which(indxbuf==indx1[i])]
        n2 = n_slct_snp[which(indxbuf==indx2[i])]
        if(n1 < n2) {
            t = indx1[i]; indx1[i] = indx2[i]; indx2[i] = t
        }
    }
    indx1 = unique(sort(indx1))
    return(indx1)
}


# ************************************************** #
#              Filter the summary data               #
# ************************************************** #
# filter the input to se and zscore
# parameters see HEIDI or SMR_Multi functions
# return index: the index be used for further analysis
# NOTICE: the parameters will be revised in place after run!!!!! TAKE CARE
filter_summdat <- function(snp_id, bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, n_ref, nsnps_thresh, pvalue_thresh, ld_r2_thresh, fdr_thresh){
    # make sure they are numeric numbers
    m = length(bzx)
    message("There are ", m, " SNPs in the dataset.")
    remain_index = seq(1,m)
    snpIDp = as.character(snp_id);
    bZXp = as.numeric(as.character(bzx))
    seZXp = as.numeric(as.character(bzx_se));
    pZXp = as.numeric(as.character(bzx_pval));
    bZYp = as.numeric(as.character(bzy))
    seZYp = as.numeric(as.character(bzy_se));
    ldrhop = matrix(as.numeric(as.character(unlist(ldrho))), m, m)

    # remove SNPs with missing value
    na_snps=c()
    indx = which(!is.finite(bZXp) | !is.finite(seZXp) | !is.finite(pZXp) | !is.finite(bZYp) | !is.finite(seZYp) | 
                 is.na(bZXp) | is.na(seZXp) | is.na(pZXp) | is.na(bZYp) | is.na(seZYp))
    if(length(indx)>0) {
        na_snps = snpIDp[indx];
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        snpIDp = snpIDp[-indx] 
        warning(length(indx), " SNPs were removed due to missing estimates in the summary data.")
    }
    m = length(bZXp)
    if(m < nsnps_thresh) stop("At least ", nsnps_thresh, " SNPs are required. Note: this hard limit can be changed by the \"nsnps_thresh\".");

    # remove SNPs with very small SE
    indx = which(seZXp < eps | seZYp < eps )
    if(length(indx)>0) {
        na_snps = c(na_snps, snpIDp[indx]);
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        snpIDp = snpIDp[-indx]
        warning(length(indx), " SNPs were removed due to extremely small standard error. Please check that data.")
    }
    m = length(bZXp)
    if(m < nsnps_thresh) stop("At least ", nsnps_thresh, " SNPs are required. Note: this hard limit can be changed by the \"nsnps_thresh\"."); 
    # remove SNPs with missing LD
    indx = which(is.na(ldrhop[upper.tri(ldrhop)]))
    if(length(indx) > 0)
        stop("LD correlations between ", length(indx), " pairs of SNPs are missing. Please check the MAF of the SNPs and the missingness rate in the reference sample.")

    # z score of bzx
    weak_snps = c()
    indx = which(pZXp > pvalue_thresh)
    if(length(indx)>0) {
        weak_snps = snpIDp[indx];
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx];
        snpIDp = snpIDp[-indx];
        warning(length(indx), " non-significant SNPs were removed.")
    }
    m = length(bZXp)
    if(m < nsnps_thresh) stop("At least ", nsnps_thresh, " SNPs are required. Note: this hard limit can be changed by the \"nsnps_thresh\"."); 

    # check LD r
    linkage_snps = c()
    indx = snp_ld_prune(ldrhop, ld_r2_thresh)
    if(length(indx) > 0) {
        linkage_snps = snpIDp[indx]
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        snpIDp = snpIDp[-indx]
        warning("There were SNPs in high LD. After LD pruning with a LD r2 threshold of ", ld_r2_thresh, ", ", length(indx), " SNPs were removed. Note: The threshold of LD can be changed by the \"ld_r2_thresh\".")
    }
    m = length(bZXp)
    if(m < nsnps_thresh) stop("At least ", nsnps_thresh, " SNPs are required. Note: this hard limit can be changed by the \"nsnps_thresh\"."); 
    
    # update LD correlation matrix
    var_rho = 1/n_ref
    pval_rho = pchisq(ldrhop[upper.tri(ldrhop)]^2/var_rho, 1, lower.tail=F)
    qval_rho = p.adjust(pval_rho, method = "fdr")
    qval_mat = matrix(0, m, m)
    qval_mat[upper.tri(qval_mat)] = qval_rho; qval_mat = t(qval_mat); qval_mat[upper.tri(qval_mat)] = qval_rho;
    rho_index = which(qval_mat >= fdr_thresh, arr.ind=T)
    ldrhop[rho_index] = 0 

    message(length(remain_index), " SNPs were retained after filtering.")

    # replace the parameters in place
    eval.parent(substitute(bzx <- bZXp))
    eval.parent(substitute(bzy <- bZYp))
    eval.parent(substitute(bzx_se <- seZXp))
    eval.parent(substitute(bzy_se <- seZYp))
    eval.parent(substitute(bzx_pval <- pZXp))
    eval.parent(substitute(ldrho <- ldrhop))
    return(list(remain_index=remain_index, na_snps=na_snps, weak_snps=weak_snps, linkage_snps=linkage_snps))
}

# ************************************************** #
#         standardization of b and s.e.              #
# ************************************************** #
#' @title Standardization of effect size and its standard error
#' @description Standardization of SNP effect and its standard error using z-statistic, allele frequency and sample size
#' @usage std_effect(snp_freq, b, se, n)
#' @param snp_freq vector, allele frequencies
#' @param b vector, SNP effects on risk factor
#' @param se vector, standard errors of b
#' @param n vector, per-SNP sample sizes for GWAS of the risk factor
#' @examples
#' data("gsmr")
#' std_effects = std_effect(gsmr_data$a1_freq, gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_n)
#'
#' @return Standardised effect (b) and standard error (se)
#' @export
std_effect <- function(snp_freq, b, se, n) {
    # double check the counts
    # if data set is empty
    check_element(snp_freq,"Minor Allele Frequency")
    check_element(b,"Effect Size")
    check_element(se,"Standard Error")
    check_element(n,"Sample Size")

    message("std effect: ", length(b), " instruments loaded.")

    # length is different or not
    len_vec <- c(length(snp_freq),length(b),length(se),length(n))
    if (!check_vec_elements_eq(len_vec)){
        stop("Lengths of the input vectors are different. Please check.");
    }
  
    # check missing values
    indx = which(!is.finite(snp_freq) | !is.finite(b) | !is.finite(se) | !is.finite(n) | is.na(snp_freq) | is.na(b) | is.na(se) | is.na(n))
    if (length(indx)>0) {
        stop("There are ", length(indx), " SNPs with missing estimates in the summary data. Please check.");
    } 

    # make sure they are numeric numbers
    snpfreq = as.numeric(as.character(snp_freq))
    b = as.numeric(as.character(b))
    se = as.numeric(as.character(se))
    n = as.numeric(as.character(n))

    zscore = b / se
    b_p = zscore / sqrt(2*snp_freq*(1-snp_freq)*(n+zscore^2))
    se_p = 1 / sqrt(2*snp_freq*(1-snp_freq)*(n+zscore^2))

    return(list(b=b_p,se=se_p))
}


# ************************************************** #
#                   GSMR analysis                    #
# ************************************************** #
#' @title Generalized Summary-data-based Mendelian Randomization analysis
#' @description GSMR (Generalised Summary-data-based Mendelian Randomisation) is a flexible and powerful approach that utilises multiple genetic instruments to test for causal association between a risk factor and disease using summary-level data from independent genome-wide association studies.
#' @usage gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10)
#' @param bzx vector, SNP effects on risk factor
#' @param bzx_se vector, standard errors of bzx
#' @param bzx_pval vector, p values for bzx
#' @param bzy vector, SNP effects on disease
#' @param bzy_se vector, standard errors of bzy
#' @param ldrho LD correlation matrix of the SNPs
#' @param snpid genetic instruments
#' @param n_ref sample size of the reference sample
#' @param heidi_outlier_flag flag for HEIDI-outlier analysis
#' @param gwas_thresh threshold p-value to select instruments from GWAS for risk factor
#' @param heidi_outlier_thresh HEIDI-outlier threshold 
#' @param nsnps_thresh the minimum number of instruments required for the GSMR analysis (we do not recommend users to set this number smaller than 10)
#' @param ld_r2_thresh LD r2 threshold to remove SNPs in high LD
#' @param ld_fdr_thresh FDR threshold to remove the chance correlations between SNP instruments 
#' @examples
#' data("gsmr")
#' gsmr_result = gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, ldrho, gsmr_data$SNP, n_ref, T, 5e-8, 0.01, 10, 0.1, 0.05) 
#'
#' @return Estimate of causative effect of risk factor on disease (bxy), the corresponding standard error (bxy_se), p-value (bxy_pval), SNP index (used_index), SNPs with missing values, with non-significant p-values and those in LD.
#' @export
gsmr <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, n_ref,
                 heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_thresh=0.05, 
                 nsnps_thresh=10, ld_r2_thresh=0.05, ld_fdr_thresh=0.05) {
    global_heidi_thresh = heidi_thresh;
    # subset of LD r matrix
    len1 = length(Reduce(intersect, list(snpid, colnames(ldrho))))
    len2 = length(snpid)
    len_vec <- c(len1, len2)
    if (!check_vec_elements_eq(len_vec)){
        stop(paste(len2 - len1, " SNPs are missing in the LD correlation matrix. Please check.", sep=""))
    }
    ldrho = ldrho[snpid, snpid]
    # double check the counts
    len_vec <- c(length(bzx),length(bzx_se),length(bzy),length(bzy_se),dim(ldrho)[1],dim(ldrho)[2])
    if (!check_vec_elements_eq(len_vec)){
        stop("Lengths of the input vectors are different. Please check.");
    }
    message("GSMR analysis: ", length(bzx), " instruments loaded.")
    # filter dataset
    resbuf <- filter_summdat(snpid, bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, n_ref, nsnps_thresh, gwas_thresh, ld_r2_thresh, ld_fdr_thresh)
    remain_index<-resbuf$remain_index;
    na_snps<-resbuf$na_snps; weak_snps<-resbuf$weak_snps; linkage_snps<-resbuf$linkage_snps;
    if(length(remain_index) < nsnps_thresh) {
      stop("Not enough SNPs for the GSMR analysis. At least ", nsnps_thresh, " SNPs are required. Note: this hard limit can be changed by the \"nsnps_thresh\".");
    }
    pleio_snps=NULL;
    if(heidi_outlier_flag ) {
        # Standard HEIDI-outlier
        nsnp = length(bzx); kept_index = c(1:nsnp)
        while(1) {
            nsnp_iter = length(kept_index)
            std_het_pval = std_heidi_outlier(bzx[kept_index], bzx_se[kept_index], bzy[kept_index], bzy_se[kept_index], ldrho[kept_index,kept_index])$pval_het
            indi_heidi_thresh = min(c(0.01, 0.05/nsnp_iter))
            if(min(std_het_pval) >= indi_heidi_thresh) break;
            excl_index = which.min(std_het_pval)
            pleio_snps = c(pleio_snps, remain_index[kept_index[excl_index]])
            kept_index = kept_index[-excl_index]
        }
        # update summary stats
        bzx = bzx[kept_index]; bzx_se = bzx_se[kept_index]; bzx_pval = bzx_pval[kept_index];
        bzy = bzy[kept_index]; bzy_se = bzy_se[kept_index]; ldrho = ldrho[kept_index, kept_index];

        # Global HEIDI-outlier
        resbuf = global_heidi_outlier(bzx, bzx_se, bzy, bzy_se, ldrho, global_heidi_thresh)
        bxy_hat = resbuf$bxy; bxy_hat_se = resbuf$bxy_se; bxy_hat_pval = resbuf$bxy_pval
        global_het_pval = resbuf$pval_het;
        if(length(resbuf$pleio_index) > 0) pleio_snps = c(pleio_snps, remain_index[resbuf$pleio_index])
        if(length(pleio_snps) > 0) {
            remain_index = remain_index[-pleio_snps]
            pleio_snps = snpid[pleio_snps]
        }
    } else {
        # Estimate bxy by GSMR
        resbuf = bxy_gsmr(bzx, bzx_se, bzy, bzy_se, ldrho);
        bxy_hat = resbuf$bxy; bxy_hat_se = resbuf$bxy_se; bxy_hat_pval = resbuf$bxy_pval;
        vec_1t_v = resbuf$vec_1t_v;
        global_het_pval = global_heidi_pvalue(bxy_hat, bxy_hat_se, vec_1t_v, bzx, bzx_se, bzy, bzy_se, ldrho);
    }
    return(list(bxy=bxy_hat, bxy_se=bxy_hat_se, bxy_pval=bxy_hat_pval, used_index=remain_index,
               global_het_pval=global_het_pval, na_snps=na_snps, weak_snps=weak_snps, linkage_snps=linkage_snps,
               pleio_snps=pleio_snps))
}

# ************************************************** #
#            Bi-directional GSMR analysis            #
# ************************************************** #
#' @title Bi-directional GSMR analysis
#' @description Bi-directional GSMR analysis is composed of a forward-GSMR analysis and a reverse-GSMR analysis that uses SNPs associated with the disease (e.g. at  < 5e-8) as the instruments to test for putative causal effect of the disease on the risk factor.
#' @usage bi_gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snpid, heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10)
#' @param bzx vector, SNP effects on risk factor
#' @param bzx_se vector, standard errors of bzx
#' @param bzx_pval vector, p values for bzx
#' @param bzy vector, SNP effects on disease
#' @param bzy_se vector, standard errors of bzy
#' @param bzy_pval vector, p values for bzy
#' @param ldrho LD correlation matrix of the SNPs
#' @param snpid genetic instruments
#' @param n_ref sample size of the reference sample
#' @param heidi_outlier_flag flag for HEIDI-outlier analysis
#' @param gwas_thresh threshold p-value to select  instruments from GWAS for risk factor
#' @param heidi_outlier_thresh HEIDI-outlier threshold 
#' @param nsnps_thresh the minimum number of instruments required for the GSMR analysis (we do not recommend users to set this number smaller than 10)
#' @param ld_r2_thresh LD r2 threshold to remove SNPs in high LD
#' @param ld_fdr_thresh FDR threshold to remove the chance correlations between SNP instruments
#' @examples
#' data("gsmr")
#' gsmr_result = bi_gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval, ldrho, gsmr_data$SNP, n_ref, T, 5e-8, 0.01, 10, 0.1, 0.05) 
#'
#' @return Estimate of causative effect of risk factor on disease (forward_bxy), the corresponding standard error (forward_bxy_se), p-value (forward_bxy_pval) and SNP index (forward_index), and estimate of causative effect of disease on risk factor (reverse_bxy), the corresponding standard error (reverse_bxy_se), p-value (reverse_bxy_pval), SNP index (reverse_index), SNPs with missing values, with non-significant p-values and those in LD.
#' @export
bi_gsmr <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snpid, n_ref,
               heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, 
               nsnps_thresh=10, ld_r2_thresh=0.05, ld_fdr_thresh=0.05) {
    ## Forward GSMR
    message("Forward GSMR analysis...")   
    gsmr_result=gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, n_ref, heidi_outlier_flag, gwas_thresh, heidi_outlier_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)
    bxy1 = gsmr_result$bxy; bxy1_se = gsmr_result$bxy_se; bxy1_pval = gsmr_result$bxy_pval;
    bxy1_index = gsmr_result$used_index;
    na_snps1 = gsmr_result$na_snps; weak_snps1 = gsmr_result$weak_snps; 
    linkage_snps1 = gsmr_result$linkage_snps; pleio_snps1 = gsmr_result$pleio_snps;

    ## Reverse GSMR
    message("Reverse GSMR analysis...")           
    gsmr_result=gsmr(bzy, bzy_se, bzy_pval, bzx, bzx_se, ldrho, snpid, n_ref, heidi_outlier_flag, gwas_thresh, heidi_outlier_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)
    bxy2 = gsmr_result$bxy; bxy2_se = gsmr_result$bxy_se; bxy2_pval = gsmr_result$bxy_pval;
    bxy2_index = gsmr_result$used_index;
    na_snps2 = gsmr_result$na_snps; 
    weak_snps2 = gsmr_result$weak_snps; 
    linkage_snps2 = gsmr_result$linkage_snps;
    pleio_snps2 = gsmr_result$pleio_snps;
    return(list(forward_bxy=bxy1, forward_bxy_se=bxy1_se, 
                forward_bxy_pval=bxy1_pval, forward_index=bxy1_index,
                reverse_bxy=bxy2, reverse_bxy_se=bxy2_se,             
                reverse_bxy_pval=bxy2_pval, reverse_index=bxy2_index,
                forward_na_snps=na_snps1, forward_weak_snps=weak_snps1, 
                forward_linkage_snps=linkage_snps1, forward_pleio_snps=pleio_snps1,
                reverse_na_snps=na_snps2, reverse_weak_snps=weak_snps2, 
                reverse_linkage_snps=linkage_snps2, reverse_pleio_snps=pleio_snps2))
}

bi_gsmr_v2_beta <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snpid, n_ref,
               heidi_outlier_flag=T, gwas_thresh=5e-8, global_heidi_thresh=0.05, indi_heidi_thresh=0.01, 
               nsnps_thresh=10, ld_r2_thresh=0.05, ld_fdr_thresh=0.05) {
    ## Forward GSMR
    message("Forward GSMR analysis...")
    gsmr_result=gsmr_v2_beta(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, n_ref, heidi_outlier_flag, gwas_thresh, global_heidi_thresh, indi_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)
    bxy1 = gsmr_result$bxy; bxy1_se = gsmr_result$bxy_se; bxy1_pval = gsmr_result$bxy_pval;
    bxy1_index = gsmr_result$used_index;
    na_snps1 = gsmr_result$na_snps; 
    weak_snps1 = gsmr_result$weak_snps; 
    linkage_snps1 = gsmr_result$linkage_snps; 
    pleio_snps1 = gsmr_result$pleio_snps;

    ## Reverse GSMR
    message("Reverse GSMR analysis...")
    gsmr_result=gsmr_v2_beta(bzy, bzy_se, bzy_pval, bzx, bzx_se, ldrho, snpid, n_ref, heidi_outlier_flag, gwas_thresh, global_heidi_thresh, indi_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh)
    bxy2 = gsmr_result$bxy; bxy2_se = gsmr_result$bxy_se; bxy2_pval = gsmr_result$bxy_pval;
    bxy2_index = gsmr_result$used_index;
    na_snps2 = gsmr_result$na_snps;
    weak_snps2 = gsmr_result$weak_snps;
    linkage_snps2 = gsmr_result$linkage_snps;
    pleio_snps2 = gsmr_result$pleio_snps;
    return(list(forward_bxy=bxy1, forward_bxy_se=bxy1_se,
                forward_bxy_pval=bxy1_pval, forward_index=bxy1_index,
                reverse_bxy=bxy2, reverse_bxy_se=bxy2_se,
                reverse_bxy_pval=bxy2_pval, reverse_index=bxy2_index,
                forward_na_snps=na_snps1, forward_weak_snps=weak_snps1, 
                forward_linkage_snps=linkage_snps1, forward_pleio_snps=pleio_snps1,
                reverse_na_snps=na_snps2, reverse_weak_snps=weak_snps2,
                reverse_linkage_snps=linkage_snps2, reverse_pleio_snps=pleio_snps2))
}

