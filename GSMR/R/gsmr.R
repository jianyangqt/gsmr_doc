# gsmr: A tool for GSMR and HEIDI analysis
# The gsmr package Perform Generalized Summary-data-based Mendelian Randomization analysis (GSMR)
#  and HEterogeneity In Dependent Instruments analysis to remove pleiotropic outliers (HEIDI-outlier).
# @author Zhihong Zhu <z.zhu1@uq.edu.au>
# @author Zhili Zheng <zhilizheng@outlook.com>
# @author Futao Zhang <f.zhang5@uq.edu.au>
# @author Jian Yang <jian.yang@uq.edu.au>

eps = 1e-6;

# ************************************************** #
#         check data is missing or not               #
# ************************************************** #
check_element <- function(vals, argue) {
    vals = vals[which(is.finite(vals))]
    if(length(vals)==0) {
        stop("None of the ", argue, " is found. Please double check.")
    }
}

# ************************************************** #
#       variance-covariane matrix of bXY             #
# ************************************************** #
cov_bXY <- function(bzx, bzx_se, bzy, bzy_se, ldrho) {
    bXY = bzy / bzx
    zscoreZX = bzx / bzx_se
    nsnp = dim(ldrho)[1]
    covbXY = diag(nsnp)
    if(nsnp>1) {
        zszxij = zscoreZX%*%t(zscoreZX)
        bxyij = bXY%*%t(bXY)
        sezyij = bzy_se%*%t(bzy_se)
        bzxij = bzx%*%t(bzx)
        covbXY = ldrho*sezyij/bzxij + ldrho*bxyij/zszxij - bxyij/zszxij^2
    }
    return(covbXY)
}

# ************************************************** #
#                 HEIDI test                         #
# ************************************************** #
#' @importFrom survey pchisqsum
heidi <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho,
                 heidi_thresh = pchisq(10, 1, lower.tail=F),
                 nSNPs_thresh=10, maxid=integer(0), gwas_flag=TRUE) {

    # double check the counts
    len_vec <- c(length(bzx),length(bzx_se),length(bzy),length(bzy_se),dim(ldrho)[1],dim(ldrho)[2])
    if (!check_vec_elements_eq(len_vec)){
        stop("Lengths of input variables are different. Please double check.");
    }
    if(len_vec[1] < nSNPs_thresh) {
        stop("More than ", nSNPs_thresh, " genetic instruments are required in the HEIDI test.")
    }
    
    remain_index = seq(1, length(bzx))
    if(gwas_flag) {
        remain_index <- filter_summdat(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, heidi_thresh)
    }
    m = length(remain_index)

    # recaculate the zscore_ZX for vector has been changed
    zscore_ZX = bzx / bzx_se
    bXY = bzy / bzx
    seSMR = sqrt( (bzy_se^2*bzx^2 + bzx_se^2*bzy^2) / bzx^4 )
    # remap the maxid according to filtered data
    maxid = which(remain_index==maxid)
    if(length(maxid) != 1) {
        stop("The reference SNP for HEIDI test is missing.")
    }
    # diff = bXY_top - bXY_-i
    dev = bXY[maxid] - bXY[-maxid]
    # v matrix
    covbXY = cov_bXY(bzx, bzx_se, bzy, bzy_se, ldrho)
    tmp1 = diag(covbXY)[maxid]
    tmp2 = covbXY[-maxid, -maxid]
    tmp3 = covbXY[maxid, -maxid]
    vdev = tmp1 + tmp2 - tmp3
    vdev = t(t(vdev) - tmp3)
    diag(vdev) = diag(vdev) + eps
    # variance of diff
    vardev = diag(covbXY)[-maxid] + tmp1 - 2*tmp3
    vardev = vardev + eps
    chisq_dev = dev^2 / vardev
    if(m>2) {
        # correlation matrix of diff
        corr_dev = diag(m-1)
        # more than 2 instruments
        for( i in 1 : (m-2) ) {
            for( j in (i+1) : (m-1) ) {
                corr_dev[i,j] = corr_dev[j,i] =
                       vdev[i,j] / sqrt(vdev[i,i]*vdev[j,j])
            }
        }

        # estimate the p value
        lambda = eigen(corr_dev, symmetric=TRUE, only.values=TRUE)$values
        t = length(lambda)
        pHet = pchisqsum(sum(chisq_dev)+eps, df=rep(1,t), a=lambda, method="sadd", lower.tail=F)

    } else {
        pHet = pchisq(chisq_dev, 1, lower.tail=F)
    }

    return(list(pheidi=pHet, nsnps=m, used_index=remain_index))
}

# ************************************************** #
#          Iterations for HEIDI-ouliter              #
# ************************************************** #

heidi_outlier_iter <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh, heidi_thresh, remain_index) {
    # remove pleiotropic instruments
    m = length(remain_index)
    # grab reference instrument
    bxy = bzy/bzx
    bxy_q = quantile(bxy, probs = seq(0, 1, 0.2))
    min_bxy = as.numeric(bxy_q[3]); max_bxy = as.numeric(bxy_q[4]);
    slctindx = which(bxy <= max_bxy & bxy >= min_bxy)
    refid = slctindx[which.min(bzx_pval[slctindx])]
    pheidi = as.numeric()
    if(length(slctindx)==0) {
      stop("None SNPs were retained by HEIDI-outlier test.");
    }
    for( i in 1 : m ) {
      if( i==refid ) next
      heidi_result = heidi(bzx[c(refid,i)], bzx_se[c(refid,i)], bzx_pval[c(refid,i)],
                          bzy[c(refid,i)], bzy_se[c(refid,i)],
                          ldrho[c(refid,i), c(refid,i)],
                          gwas_thresh, 2, 1, FALSE)
      pheidi[i] = as.numeric(heidi_result$pheidi)
    }
    remain_index = remain_index[sort(c(refid,which(pheidi>=heidi_thresh)))]
    return(remain_index)
}

# ************************************************** #
#              Test identical elements               #
# ************************************************** #
check_vec_elements_eq <- function(dim_vec){
    return(all(dim_vec == dim_vec[1]))
}


# ************************************************** #
#              Filter the summary data               #
# ************************************************** #
# filter the input to se and zscore
# parameters see HEIDI or SMR_Multi functions
# return index: the index be used for further analysis
# NOTICE: the parameters will be revised in place after run!!!!! TAKE CARE
filter_summdat <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, pvalue_thresh){
    # make sure they are numeric numbers
    m = length(bzx)
    message("There are ", m, " SNPs in the dataset.")
    remain_index = seq(1,m)
    bZXp = as.numeric(as.character(bzx))
    seZXp = as.numeric(as.character(bzx_se));
    pZXp = as.numeric(as.character(bzx_pval));
    bZYp = as.numeric(as.character(bzy))
    seZYp = as.numeric(as.character(bzy_se));
    ldrhop = matrix(as.numeric(as.character(unlist(ldrho))), m, m)
    # remove SNPs with missing SE
    indx = which(!is.finite(seZXp) | !is.finite(seZYp))
    if(length(indx)>0) {
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        warning(length(indx), " SNPs were removed due to missing standard error.")
    }
    # remove SNPs with very small SE
    indx = which(seZXp < eps | seZYp < eps )
    if(length(indx)>0) {
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        warning(length(indx), " SNPs were removed due to extremely small standard error.")
    }
    # z score of bzx
    indx = which(pZXp > pvalue_thresh)
    if(length(indx)>0) {
        bZXp = bZXp[-indx]; seZXp = seZXp[-indx]; pZXp = pZXp[-indx];
        bZYp = bZYp[-indx]; seZYp = seZYp[-indx];
        ldrhop = ldrhop[-indx, -indx];
        remain_index = remain_index[-indx]
        warning(length(indx), " SNPs were removed due to insignficant associations between instruments and risk factor.")
    }
    message(length(remain_index), " SNPs were retained with filtering of weak instruments.")
    # replace the parameters in place
    eval.parent(substitute(bzx <- bZXp))
    eval.parent(substitute(bzy <- bZYp))
    eval.parent(substitute(bzx_se <- seZXp))
    eval.parent(substitute(bzy_se <- seZYp))
    eval.parent(substitute(bzx_pval <- pZXp))
    eval.parent(substitute(ldrho <- ldrhop))
    return(remain_index)
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
#' std_effects = std_effect(gsmr_data$freq, gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_n)
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
        stop("Lengths of input variables are different. Please double check.");
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
#             HEIDI-outlier analysis                 #
# ************************************************** #
#' @title HEIDI-outlier analysis
#' @description An analysis to detect and eliminate from the analysis instruments that show significant pleiotropic effects on both risk factor and disease
#' @usage heidi_outlier(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10)
#' @param bzx vector, SNP effects on risk factor
#' @param bzx_se vector, standard errors of bzx
#' @param bzx_pval vector, p values for bzx
#' @param bzy vector, SNP effects on disease
#' @param bzy_se vector, standard errors of bzy
#' @param ldrho LD correlation matrix of the SNPs
#' @param snpid genetic instruments
#' @param gwas_thresh threshold p-value to select instruments from GWAS for risk factor
#' @param heidi_outlier_thresh threshold p-value to remove pleiotropic outliers (the default value is 0.01)
#' @param nsnps_thresh the minimum number of instruments required for the GSMR analysis (we do not recommend users to set this number smaller than 10)
#' @examples
#' data("gsmr")
#' filtered_index = heidi_outlier(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, ldrho, gsmr_data$SNP, 5e-8, 0.01, 10)
#'
#' @return Retained index of genetic instruments
#' @export
heidi_outlier <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid,
                      gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10) {
    # Subset of LD r matrix
    len1 = length(Reduce(intersect, list(snpid, colnames(ldrho))))
    len2 = length(snpid)
    len_vec <- c(len1, len2)    
    if (!check_vec_elements_eq(len_vec)){ 
        stop(paste(len2 - len1, " SNPs are missing in the LD correlation matrix. Please double check.", sep=""))
    }
    ldrho = ldrho[snpid, snpid]
    # double check the counts
    len_vec <- c(length(bzx),length(bzx_se),length(bzx_pval),length(bzy),length(bzy_se),dim(ldrho)[1],dim(ldrho)[2])
    if (!check_vec_elements_eq(len_vec)){
      stop("Lengths of input variables are different. Please double check.");
    }
    message("HEIDI-outlier: ", length(bzx), " instruments loaded.")
    # filter dataset
    remain_index <- filter_summdat(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh)
    if(length(remain_index) < nsnps_thresh) {
      stop("HEIDI-outlier stopped because the number of instruments < ", nsnps_thresh, ".");
    }
    # Perform HEIDI-outlier
    remain_index <- heidi_outlier_iter(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh, heidi_outlier_thresh, remain_index)
    message(length(remain_index), " SNPs were retained by HEIDI-outlier test.")
    return(remain_index)
}

# ************************************************** #
#                   GSMR analysis                    #
# ************************************************** #
#' @title Generalized Summary-data-based Mendelian Randomization analysis
#' @description GSMR (Generalised Summary-data-based Mendelian Randomisation) is a flexible and powerful approach that utilises multiple genetic instruments to test for putative causal association between a risk factor and disease using summary-level data from independent genome-wide association studies.
#' @usage gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10)
#' @param bzx vector, SNP effects on risk factor
#' @param bzx_se vector, standard errors of bzx
#' @param bzx_pval vector, p values for bzx
#' @param bzy vector, SNP effects on disease
#' @param bzy_se vector, standard errors of bzy
#' @param ldrho LD correlation matrix of the SNPs
#' @param snpid genetic instruments
#' @param heidi_outlier_flag flag for HEIDI-outlier analysis
#' @param gwas_thresh threshold p-value to select instruments from GWAS for risk factor
#' @param heidi_outlier_thresh HEIDI-outlier threshold 
#' @param nsnps_thresh the minimum number of instruments required for the GSMR analysis (we do not recommend users to set this number smaller than 10)
#' @examples
#' data("gsmr")
#' gsmr_result = gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, ldrho, gsmr_data$SNP, T, 5e-8, 0.01, 10) 
#'
#' @return Estimate of causative effect of risk factor on disease (bxy), the corresponding standard error (bxy_se), p-value (bxy_pval) and SNP index (used_index).
#' @export
gsmr <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se,
                ldrho, snpid, heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10) {
    # subset of LD r matrix
    len1 = length(Reduce(intersect, list(snpid, colnames(ldrho))))
    len2 = length(snpid)
    len_vec <- c(len1, len2)    
    if (!check_vec_elements_eq(len_vec)){
        stop(paste(len2 - len1, " SNPs are missing in the LD correlation matrix. Please double check.", sep=""))
    }
    ldrho = ldrho[snpid, snpid]
    # double check the counts
    len_vec <- c(length(bzx),length(bzx_se),length(bzy),length(bzy_se),dim(ldrho)[1],dim(ldrho)[2])
    if (!check_vec_elements_eq(len_vec)){
        stop("Lengths of input variables are different. Please double check.");
    }
    message("GSMR analysis: ", length(bzx), " instruments loaded.")
    # filter dataset
    remain_index <- filter_summdat(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh)
    if(length(remain_index) < nsnps_thresh) {
      stop("GSMR analysis stopped because the number of instruments < ", nsnps_thresh, ".");
    }
    if(heidi_outlier_flag) {
        # Perform HEIDI-outlier
        remain_index2 <- heidi_outlier_iter(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, gwas_thresh, heidi_outlier_thresh, remain_index)
        if(length(remain_index2) < nsnps_thresh) {
            stop("GSMR analysis stopped because the number of instruments < ", nsnps_thresh, ".");
        } else {
            message(length(remain_index2), " SNPs were retained by HEIDI-outlier test.") 
        }
        # Update estimates
        remain_index2 = match(remain_index2, remain_index)
        bzx = bzx[remain_index2]; bzx_se = bzx_se[remain_index2]; 
        bzy = bzy[remain_index2]; bzy_se = bzy_se[remain_index2];
        ldrho = ldrho[remain_index2,remain_index2];
    }
    # do the SMR test with multiple instruments
    message("Computing estimate of bxy at each instrument.")
    bXY = bzy/bzx

    message("Estimating variance-covariance matrix of bxy.")
    covbXY = cov_bXY(bzx, bzx_se, bzy, bzy_se, ldrho)
    diag(covbXY) = diag(covbXY) + eps
    # Eigen decomposition
    resbuf = eigen(covbXY, symmetric=TRUE)
    eval = as.numeric(resbuf$values)
    evec = resbuf$vectors
    if(min(abs(eval)) < eps) {
      stop("The variance-covariance matrix of bxy is not invertable!");
    }
    covbXY_inv = evec%*%diag(1/eval)%*%t(evec)
    message("Estimating bxy using all the instruments.")
    vec_1 = rep(1, length(bzx))
    num_1_v_1 = as.numeric(solve(t(vec_1)%*%covbXY_inv%*%vec_1))
    vec_1_v = as.numeric(t(vec_1)%*%covbXY_inv)
    bXY_GLS = num_1_v_1*vec_1_v%*%bXY
    varbXY_GLS = num_1_v_1
    chisqbXY_GLS = bXY_GLS^2/varbXY_GLS
    pbXY_GLS = pchisq(chisqbXY_GLS, 1, lower.tail=F)
    message("GSMR analysis is completed.")
    return(list(bxy=bXY_GLS, bxy_se=sqrt(varbXY_GLS), bxy_pval=pbXY_GLS,
                used_index=remain_index))
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
#' @param heidi_outlier_flag flag for HEIDI-outlier analysis
#' @param gwas_thresh threshold p-value to select  instruments from GWAS for risk factor
#' @param heidi_outlier_thresh HEIDI-outlier threshold 
#' @param nsnps_thresh the minimum number of instruments required for the GSMR analysis (we do not recommend users to set this number smaller than 10)
#' @examples
#' data("gsmr")
#' gsmr_result = bi_gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval, ldrho, gsmr_data$SNP, T, 5e-8, 0.01, 10) 
#'
#' @return Estimate of causative effect of risk factor on disease (forward_bxy), the corresponding standard error (forward_bxy_se), p-value (forward_bxy_pval) and SNP index (forward_index), and estimate of causative effect of disease on risk factor (reverse_bxy), the corresponding standard error (reverse_bxy_se), p-value (reverse_bxy_pval) and SNP index (reverse_index).
#' @export
bi_gsmr <- function(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval,
                ldrho, snpid, heidi_outlier_flag=T, gwas_thresh=5e-8, heidi_outlier_thresh=0.01, nsnps_thresh=10) {
    ## Forward GSMR
    message("Forward GSMR analysis...")   
    gsmr_result=gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, ldrho, snpid, heidi_outlier_flag, gwas_thresh, heidi_outlier_thresh, nsnps_thresh)
    bxy1 = gsmr_result$bxy; bxy1_se = gsmr_result$bxy_se; bxy1_pval = gsmr_result$bxy_pval;
    bxy1_index = gsmr_result$used_index;

    ## Reverse GSMR
    message("Reverse GSMR analysis...")           
    gsmr_result=gsmr(bzy, bzy_se, bzy_pval, bzx, bzx_se, ldrho, snpid, heidi_outlier_flag, gwas_thresh, heidi_outlier_thresh, nsnps_thresh)
    bxy2 = gsmr_result$bxy; bxy2_se = gsmr_result$bxy_se; bxy2_pval = gsmr_result$bxy_pval;
    bxy2_index = gsmr_result$used_index;
    return(list(forward_bxy=bxy1, forward_bxy_se=bxy1_se, 
                forward_bxy_pval=bxy1_pval, forward_index=bxy1_index,
                reverse_bxy=bxy2, reverse_bxy_se=bxy2_se,             
                reverse_bxy_pval=bxy2_pval, reverse_index=bxy2_index))
}


