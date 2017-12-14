# gsmr_doc

## Install dependencies
* R
* pandoc
* nodejs
* wkhtmltopdf

## Install
```
npm install
npm run initR
```
If there are some error, install those dependencies

## How to compile document

### Update code

GSMR/R/*.R: R codes

GSMR/DESCRPTION:  package description, version number, dates

GSMR/data/gsmr.rdata:  demo data for document

GSMR/test_data.zip:  toy data, it will move to web page directly


### Run build
```
npm run build
```

The document is located in build/GSMR_web_version.zip


### Update github
```
git add .
git commit -a -m "version message"
git push
```
