mkdir -p build
mkdir -p build/GSMR
mkdir -p build/GSMR_web
mkdir -p build/GSMR_web/static/img

cp -r ./GSMR build/
cp -r ./template/* build/GSMR_web/
mv build/GSMR/*.zip build/GSMR_web/static/

cd build
Rscript -e "roxygen2::roxygenize('./GSMR', roclets=c('rd', 'collate', 'namespace'))"
node ../depends.js
Rscript .custom_install.R

R CMD build GSMR
GSMR_pack=`ls -t *.tar.gz | head -1`
R CMD install $GSMR_pack
# install.packages("./gsmr_doc-master/gsmr_1.0.4.tar.gz", repos = NULL, type = "source")
mv $GSMR_pack ../previous/package/

ver=`sed -n 's/gsmr_\(.*\).tar.gz/\1/p' <<< ${GSMR_pack}`

Rscript -e 'rmarkdown::render("GSMR/vignettes/GSMR-intro.Rmd", "html_fragment")'
mv GSMR/vignettes/GSMR-intro.html GSMR_web/


node ../parseDoc.js
sed "s/\${ver}/${ver}/g" ./GSMR_web/index.html > ./GSMR_web/index.html.tmp
mv ./GSMR_web/index.html.tmp ./GSMR_web/index.html
npm run mini

wkhtmltopdf --zoom 10 ./GSMR_web/index.html ../previous/package/gsmr_doc_${ver}.pdf

# zip
rm ./GSMR_web/GSMR-intro.html

cp ../previous/package/* GSMR_web/static/

zip -r GSMR_web_${ver}.zip GSMR_web
cp GSMR_web_${ver}.zip ../previous/web/

echo "######################################"
echo "Build success"
echo "Deploy the web using build/GSMR_web_${ver}.zip" 

