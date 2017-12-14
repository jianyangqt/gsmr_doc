# gsmr_doc

## Install dependencies
* [R](https://cloud.r-project.org/)
* [pandoc](https://pandoc.org/installing.html)
* [nodejs](https://nodejs.org/en/download/)
* [wkhtmltopdf](https://wkhtmltopdf.org/downloads.html)

These packages have to be installed once.

## Install
```
npm install
npm run initR
```
Run in terminal and run these command once. **The above steps have to be run at least once.** 

## How to compile document

### Update code

GSMR/R/\*.R: R codes. The package document is inlined with each function. It will update to R package and web page by build system.

GSMR/vignettes/GSMR-intro.Rmd: Main document of the web site locates all here. It is same with package introduction.  If you want to attach additional files, just put the files as zip format in GSMR/ and put link in the document ./static/xx.zip. The build system will recognize them. 

GSMR/DESCRPTION:  package description, pay attention to version number, dates. Version number is the only ID to update package. **The version number shall be updated at first when update the package**.

GSMR/data/gsmr.rdata:  demo data for document. We shall not put too much data, as it will make a huge size R package.

GSMR/test_data.zip:  toy data, it will move to web page directly.

### Run build
```
npm run build
```

The document is located in build/GSMR_web_version.zip. This folder is not pushed into the repository. All built document can be found in previous/web/gsmr_web_version.zip. We shall not remove the /previous folder, it keeps track of all previous document, and is used when build a new document.  

### Update github
```
git add .
git commit -a -m "version message"
git push
```
