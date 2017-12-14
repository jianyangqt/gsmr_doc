let fs = require("fs");
let base64ToImage = require('base64-to-image');

console.log("===========================");
console.log("Generating HTML document");
console.log("---------------------------");

const image_path = './GSMR_web/static/img/';

let htmlSource = fs.readFileSync("GSMR_web/GSMR-intro.html", "utf8");
let pre_pos = 0;
const mark = '<img src="';
let num_image = 0;
let out_html = "";
while(true){
    let pos = htmlSource.indexOf('<img src="', pre_pos);
    if(pos == -1){
        out_html += htmlSource.substring(pre_pos);
        break;
    }
    console.log("Find a image at " + pos);
    out_html += htmlSource.substring(pre_pos, pos + mark.length);
    out_html += "./static/img/plot" + num_image + ".png";
    //console.log(out_html);
    pre_pos = pos + mark.length;
    let end_pos = htmlSource.indexOf('"', pre_pos);
    let temp_image = htmlSource.substring(pre_pos, end_pos);
    base64ToImage(temp_image, image_path, {'fileName': 'plot' + num_image + '.png', 'type':'png'});
    num_image++;
    //console.log(temp_image);
    pre_pos = end_pos;

}

const getRDItem = (mark, doc) => {
    let pos = doc.indexOf(mark);
    let value = doc.substring(pos + mark.length, doc.indexOf("}", pos)); 
    return value.replace(/^\s*[\r\n]/gm, "").replace(/$\n/g, "");
}

const getRD2Item = (mark, doc) => {
    let items = [];
    let pre_pos = 0;
    while(true){
        let pos = doc.indexOf(mark, pre_pos);
        if(pos == -1){
            break
        }
        let pos1 = doc.indexOf("}", pos);
        let item1 = doc.substring(pos + mark.length, pos1);
        let pos2 = doc.indexOf("}", pos1 + 1);
        let item2 = doc.substring(pos1 + 2, pos2);
        items.push({item_name: item1, item_description: item2});
        pre_pos = pos2;
    }
    return items;
}

let doc = fs.readFileSync("../template/index.html", "utf8");

doc = doc.replace("${content}", out_html);

const rd_path = './GSMR/man';
let rd_files = fs.readdirSync(rd_path);
console.log("Find " + rd_files.length + " document");

let pack_doc = "";

const template = ({name, description, usage, items, value, code}) => `
                    <hr>
                    <div id="${name}" class="section level3">
                        <h3>${name}</h3>
                        <blockquote>
                            <p>${description}</p>
                        </blockquote>
                        <h4>Usage</h4>
                        <pre>${usage}</pre>
                        <h4>Arguments</h4>
                        <table summary="R argblock">
                            ${items}
                       </table>

                        <h4>Value</h4>
                        <p>${value}</p>

                        <h4>Examples</h4>
                        <pre><code>${code}</code></pre>
                    </div>
`;

const table_item = ({item_name, item_description}) => `
                            <tr valign="top">
                                <td><code>${item_name}</code></td>
                                <td>
                                    <p>${item_description}</p>
                                </td>
                            </tr>
`;
 
 

for (let rd_index in rd_files){
    let rd_file = rd_files[rd_index];
    
    if(rd_file.endsWith(".Rd")){
        if(rd_file.indexOf("-package") == -1){
            console.log(rd_file);
            let rd_doc = fs.readFileSync(rd_path + "/" + rd_file, "utf8");
            //console.log(rd_doc);
            console.log("Get parameters");
            let name = getRDItem("\\name{", rd_doc);
            let usage = getRDItem("\\usage{", rd_doc);
            let description = getRDItem("\\description{", rd_doc);
            let example = getRDItem("\\examples{", rd_doc);
            let value = getRDItem("\\value{", rd_doc);
            let items = getRD2Item("\\item{", rd_doc);

            let table_doc = items.map(table_item).join('');

            pack_doc += template({name: name, description: description, usage: usage, 
                                  items: table_doc, value: value, code: example});

        }
    }
}

//console.log(pack_doc);

doc = doc.replace('${pack_doc}', pack_doc);

fs.writeFileSync("./GSMR_web/index.html", doc);

//console.log(doc);
console.log("Create HTML document successfully");

