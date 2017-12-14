
echo "#########################################################"
echo "Init R and essential packages"
echo "--------------------------------"
echo "Only need to be performed once"

echo ""
echo "Check R and pandoc"
r_path=`command -v R`
echo $r_path
if [[ -z "$r_path" ]]; then
    echo "Error: R not fould";
    exit
fi

rscript_path=`command -v Rscript`
if [[ -z "$rscript_path" ]]; then
    echo "Error: Rscript not found";
    exit
fi

pandoc_path=`command -v pandoc`
if [[ -z "$pandoc_path" ]]; then
    echo "Error: pandoc not found";
    exit
fi

pdf_path=`command -v wkhtmltopdf`
if [[ -z "$pdf_path" ]]; then
    echo "Error: wkhtmltopdf not found";
    exit
fi

echo "Install R pacakges"
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE); install.packages(c("roxygen2", "knitr", "rmarkdown", "survey"), Sys.getenv("R_LIBS_USER"), repos="http://cloud.r-project.org");'

echo ""
echo "Check finished"
echo "#########################################################"
echo ""

