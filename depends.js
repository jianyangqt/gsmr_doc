let fs = require("fs");
console.log("#############################")
console.log("Update depend packages");

let ns = fs.readFileSync("./GSMR/NAMESPACE", "utf8");
let imports = ns.match(/import.*?\(([a-zA-Z0-9.]*).*\)/g);
var import_packages = [];
for (var import_index in imports) {
    var import_item = imports[import_index];
    import_packages.push(import_item.match(/import.*?\(([a-zA-Z0-9.]*).*\)/)[1]);
}
import_packages = [...new Set(import_packages)];
var import_string = import_packages.join(", ");

let des_file = './GSMR/DESCRIPTION';
let des = fs.readFileSync("./GSMR/DESCRIPTION", "utf8");
let contents = des.replace(/^Imports:.*$/gm ,"Imports: " + import_string);

fs.writeFileSync(des_file, contents);

let version_number = contents.match(/Version:\s*(.*)\s*\n/)[1];

var install_string = "";
var doc_install_string = "";
if(import_packages.length > 0){
    var packs = import_packages.map(function(item){
        return "'" + item + "'";
    }).join(", ");

    install_string = `install.packages(c(${packs}), Sys.getenv('R_LIBS_USER'), repos='http://cloud.r-project.org');`;
    doc_install_string = `install.packages(c(${packs}));`;
}

fs.writeFileSync(".custom_install.R", install_string); 

let intro_file = "./GSMR/vignettes/GSMR-intro.Rmd";
let intro_doc = fs.readFileSync(intro_file, "utf8");
intro_doc = intro_doc.replace(/\$\{VER\}/g, version_number).replace("${INSTALL_PACK}", doc_install_string);
console.log(version_number);
console.log(doc_install_string);

fs.writeFileSync(intro_file, intro_doc);

console.log("");
console.log("Success");
console.log("#############################")
