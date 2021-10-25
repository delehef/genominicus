const WINDOW = 15;
const GENE_WIDTH = 20;
const GENE_HEIGHT = 15;
const GENE_SPACING = 5;

MAIN_GENE = "?";

function drawGenes(it, landscape) {
    MAIN_GENE = landscape.me.name;
    var xoffset = WINDOW - landscape.lefts.length;
    for (const [i, g] of landscape.lefts.entries()) {
        it.rect(GENE_WIDTH, GENE_HEIGHT)
            .fill(g.color)
            .attr({class: "iam-" + g.name + " gene"})
            .move((i + xoffset)*(GENE_WIDTH + GENE_SPACING), 2)
        ;
    }
    it.rect(GENE_WIDTH, GENE_HEIGHT)
        .fill(landscape.me.color)
        .attr({class: "main-gene"})
        .move((WINDOW)*(GENE_WIDTH + GENE_SPACING), 2)
    ;
    for (const [i, g] of landscape.rights.entries()) {
        it.rect(GENE_WIDTH, GENE_HEIGHT)
            .fill(g.color)
            .attr({class: "iam-" + g.name + " gene"})
            .move((i + WINDOW + 1)*(GENE_WIDTH + GENE_SPACING), 2);
    }
}

function drawCluster(it, cluster) {
    var x = 0
    for (const [i, poly] of cluster.entries()) {
        var y = 0;
        let width = poly.genes.reduce((ax, g) => ax + g[1], 0);
        for (const g of poly.genes) {
            let thickness = g[1]*GENE_HEIGHT;
            if (g[0].name == MAIN_GENE) {
                it.rect(width*GENE_WIDTH, thickness)
                    .fill(g[0].color)
                    .attr({class: "main-gene"})
                    .move(x, 2 + GENE_HEIGHT-y-thickness);
            } else {
                it.rect(width*GENE_WIDTH, thickness)
                    .fill(g[0].color)
                    .attr({"stroke": "#323232", "stroke-width": 0.5, class: "gene iam-" + g[0].name})
                    .move(x, 2 + GENE_HEIGHT-y-thickness);
            }
            y += thickness;
        }
        x += width * (GENE_WIDTH + GENE_SPACING)
    }
    return x;
}

function insertAt(it, root, depth) {
    var div = document.createElement("div");
    div.classList.add("node")

    var content = document.createElement("div");
    content.classList.add("intrinsic");
    div.appendChild(content);

    if (it.isDuplication) {
        div.classList.add("D")
    }


    if (it.gene != "") {
        var geneName = document.createElement("a");
        geneName.classList.add("nametag");
        geneName.style.color = it.color;
        geneName.appendChild(document.createTextNode(it.gene))
        geneName.href = "https://www.ensembl.org/Multi/Search/Results?q="+it.gene+";site=ensembl"
        content.appendChild(geneName);

        var locus = document.createElement("a");
        locus.classList.add("nametag");
        locus.style.color = it.color;
        locus.appendChild(document.createTextNode(it.species + "/" + it.chr))
        locus.href = "https://www.genomicus.bio.ens.psl.eu/genomicus-104.02/cgi-bin/search.pl?query="+it.gene+"&view=default&nocache=$random"
        content.appendChild(locus);
    }
    if (it.clustered) {
        div.classList.add("node-container")
        var cluster = document.createElement("div");
        cluster.classList.add("cluster");
        var l = SVG().addTo(cluster);
        let width = drawCluster(l, it.clustered);
        l.size(width, GENE_HEIGHT + 4);
        content.appendChild(cluster);
    } else {
        var l = SVG().addTo(content).size(2*WINDOW*(GENE_WIDTH + GENE_SPACING) + GENE_WIDTH + 2*GENE_SPACING, GENE_HEIGHT + 4);
        drawGenes(l, it.repr);
    }


    if (it.children.length > 0) {
        div.addEventListener("click", function(e) { this.classList.toggle("condensed"); e.stopPropagation(); });
        div.style.cursor = "pointer";
        it.children.forEach((i) => insertAt(i, div, depth-1));
    }


    root.appendChild(div);
}

window.onload = () => {
    let root = document.getElementById("root")
    data.children.forEach((i) => insertAt(i, root, 10));

    Array.prototype.forEach.call(document.getElementsByClassName("gene"), function(el) {
        el.classList.forEach(className => {
            if (className.startsWith("iam-")) {
                el.addEventListener('mouseenter', e => {
                    Array.prototype.forEach.call(
                        document.getElementsByClassName(className),
                        (ee) => ee.classList.toggle("highlight-gene")
                    );
                });
                el.addEventListener('mouseleave', e => {
                    Array.prototype.forEach.call(
                        document.getElementsByClassName(className),
                        (ee) => ee.classList.toggle("highlight-gene")
                    );
                });
            }
        })
    });
}
