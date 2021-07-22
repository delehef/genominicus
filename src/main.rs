use colorsys::{Rgb, Hsl};
use crate::nhx::*;
use clap::*;
use postgres::{Client, NoTls};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use svarog::*;

const WINDOW: i64 = 15;
const GENE_WIDTH: f32 = 15.;
const GENE_SPACING: f32 = 5.;
const BRANCH_WIDTH: f32 = 20.;
const FONT_SIZE: f32 = 10.;

const ANCESTRAL_QUERY: &str = concat!(
    "select name, ancestral, annotations.species, direction, chr, start from annotations ",
    "left join mapping on name=modern where name=",
    "(select parent from annotations where name=(select parent from annotations where name=$1 limit 1))"
);
const LEFTS_QUERY: &str = concat!(
    "select coalesce(ancestral, 'NA') as ancestral, direction from ",
    "(select * from annotations left join mapping on name=modern where ",
    "annotations.species=$1 and chr=$2 and type='gene' and start<$3 order by start desc) as x limit $4"
);
const RIGHTS_QUERY: &str = concat!(
    "select coalesce(ancestral, 'NA') as ancestral, direction from ",
    "(select * from annotations left join mapping on name=modern where ",
    "annotations.species=$1 and chr=$2 and type='gene' and start>$3 order by start asc) as x limit $4"
);

mod nhx;

fn draw_background(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &Tree,
    n_: usize,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;

    for &n in tree[n_].children().iter() {
        let child = &tree[n];
        let new_y = if child.is_leaf() {
            y + 20.
        } else {
            draw_background(
                svg,
                depth,
                tree,
                n,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
            )
        };

        if tree[n_].is_duplication() {
            let d = xoffset / depth;
            svg.polygon()
                .from_pos_dims(
                    xoffset + BRANCH_WIDTH / 2.,
                    y - 6.,
                    width - xoffset - d * BRANCH_WIDTH,
                    new_y - y - 6.,
                )
                .style(|s| {
                    s.fill_color(StyleColor::Percent(0.5, 0.5, 1.))
                        .fill_opacity(0.1 + 0.9 * d)
                });
        }
        y = new_y;
    }
    y
}

fn name2color<S: AsRef<str>>(name: S) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let rgb = Rgb::from((bytes[0] as f32, bytes[1] as f32, bytes[2] as f32));

    let mut hsl: Hsl = rgb.into();
    hsl.set_lightness(hsl.lightness().clamp(30., 40.));

    let rgb: Rgb = hsl.into();
    StyleColor::Percent(rgb.red() as f32/255., rgb.green() as f32/255., rgb.blue() as f32/255.)
}
fn gene2color<S: AsRef<str>>(name: S) -> StyleColor {
    let bytes: [u8; 16] = md5::compute(name.as_ref().as_bytes()).into();
    let r = (bytes[0] as f32 / 255.).clamp(0.1, 0.9);
    let g = (bytes[1] as f32 / 255.).clamp(0.1, 0.9);
    let b = (bytes[2] as f32 / 255.).clamp(0.1, 0.9);
    StyleColor::Percent(r, g, b)
}

fn draw_tree(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
    links: &mut Vec<(Vec<String>, String, Vec<String>)>,
) -> f32 {
    fn draw_gene(
        svg: &mut SvgDrawing,
        x: f32,
        y: f32,
        right: bool,
        color: StyleColor,
    ) -> &mut Polygon {
        if right {
            svg.polygon()
                .add_point(x, y - 5.)
                .add_point(x + GENE_WIDTH - 3., y - 5.)
                .add_point(x + GENE_WIDTH, y)
                .add_point(x + GENE_WIDTH - 3., y + 5.)
                .add_point(x, y + 5.)
                .style(|s| {
                    s.fill_color(color)
                        .stroke_width(0.5)
                        .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
                })
        } else {
            svg.polygon()
                .add_point(x, y)
                .add_point(x + 3., y - 5.)
                .add_point(x + GENE_WIDTH, y - 5.)
                .add_point(x + GENE_WIDTH, y + 5.)
                .add_point(x + 3., y + 5.)
                .style(|s| {
                    s.fill_color(color)
                        .stroke_width(0.5)
                        .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
                })
        }
    }

    let mut pg = Client::connect(
        "host=localhost user=franklin dbname=duplications",
        NoTls,
    )
    .unwrap();
    let mut y = yoffset;
    let mut old_y = 0.;
    for (i, child) in node.children().iter().map(|i| &tree[*i]).enumerate() {
        if i > 0 {
            svg.line()
                .from_coords(xoffset, old_y, xoffset, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        }
        old_y = y;

        if child.is_leaf() {
            // Leaf branch
            svg.line()
                .from_coords(xoffset, y, depth, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            // Landscape support line
            svg.line()
                .start(xlabels - 5., y)
                .end(
                    xlabels + (GENE_WIDTH + GENE_SPACING) * (2. * WINDOW as f32 + 1.)
                        - GENE_SPACING
                        + 5.,
                    y,
                )
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            if let Some(name) = &child.name {
                let protein_name = name.split('_').next().unwrap();
                if let Ok(r) = pg.query_one(ANCESTRAL_QUERY, &[&protein_name]) {
                    let gene_name: &str = r.get("name");
                    let ancestral_name: &str = r.get("ancestral");
                    let species: &str = r.get("species");
                    let chr: &str = r.get("chr");
                    let pos: i32 = r.get("start");
                    let direction: &str = r.get("direction");

                    let proto_lefts = pg
                        .query(LEFTS_QUERY, &[&species, &chr, &pos, &WINDOW])
                        .unwrap()
                        .into_iter()
                        .map(|row| {
                            let ancestral: String = row.get("ancestral");
                            let direction: String = row.get("direction");
                            (ancestral, direction)
                        })
                        .collect::<Vec<_>>();
                    let proto_rights = pg
                        .query(RIGHTS_QUERY, &[&species, &chr, &pos, &WINDOW])
                        .unwrap()
                        .into_iter()
                        .map(|row| {
                            let ancestral: String = row.get("ancestral");
                            let direction: String = row.get("direction");
                            (ancestral, direction)
                        })
                        .collect::<Vec<_>>();
                    let (lefts, rights) = if direction == "+" {
                        (proto_lefts, proto_rights)
                    } else {
                        (proto_rights, proto_lefts)
                    };

                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {}/{}", protein_name, species, chr))
                        .style(|s| s.fill_color(name2color(species)));

                    // Left tail
                    let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in lefts.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == "+", gene2color(&g.0));
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                        y,
                        true,
                        gene2color(&ancestral_name),
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in rights.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        draw_gene(svg, xstart, y, g.1 == "+", gene2color(&g.0));
                    }
                    links.push((
                        lefts.iter().map(|x| x.0.clone()).collect(),
                        ancestral_name.into(),
                        rights.iter().map(|x| x.0.clone()).collect(),
                    ));
                } else {
                    // The node was not found in the database
                    eprintln!("{} not found", name);
                    links.push((Vec::new(), name.into(), Vec::new()));
                }
            } else { // The node does not have a name
            }
            y += 20.;
        } else {
            svg.line()
                .from_coords(xoffset, y, xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_tree(
                svg,
                depth,
                tree,
                child,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
                links,
            );
        }
    }

    if node.is_duplication() {
        svg.polygon()
            .from_pos_dims(xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(0.8, 0., 0.)));
    }
    y
}

fn draw_links(
    svg: &mut SvgDrawing,
    links: &[(Vec<String>, String, Vec<String>)],
    yoffset: f32,
    xlabels: f32,
) {
    let mut y = yoffset;
    for w in links.windows(2) {
        let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].0.iter().enumerate() {
            let x1 = xbase - i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].0
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase - j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .start(x1, y + 5.)
                    .end(x2, y + 20. - 5.)
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }

        let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].2.iter().enumerate() {
            let x1 = xbase + i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].2
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase + j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .start(x1, y + 5.)
                    .end(x2, y + 20. - 5.)
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }

        y += 20.;
    }
}

fn process_file(filename: &str) {
    println!("Processing {}", filename);
    let t = Tree::from_filename(filename).unwrap();

    let depth = BRANCH_WIDTH * (t.topological_depth().1 + 1.);
    let longest_name = t
        .leaf_names()
        .iter()
        .map(|(_, name)| name.map_or(0, |s| s.len()))
        .max()
        .unwrap() as f32
        * FONT_SIZE;
    let xlabels = 0.85 * (10. + depth + longest_name + 50.);
    let width = xlabels + (2. * WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING) + 60.;
    let mut svg = SvgDrawing::new();
    svg.text()
        .pos(FONT_SIZE, FONT_SIZE)
        .text(Path::new(filename).file_stem().unwrap().to_str().unwrap());

    draw_background(&mut svg, depth, &t, 0, 10.0, 50.0, xlabels, width);
    let mut links = Vec::new();
    draw_tree(
        &mut svg, depth, &t, &t[0], 10.0, 50.0, xlabels, width, &mut links,
    );
    draw_links(&mut svg, &links, 50.0, xlabels);
    svg.auto_fit();

    let mut out = File::create(&format!("{}.svg", filename)).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}

fn main() {
    let args = App::new("Genominicus")
        .version(clap::crate_version!())
        .author(clap::crate_authors!())
        .arg(
            Arg::with_name("FILE")
                .help("Sets the input file to use")
                .required(true)
                .multiple(true),
        )
        .get_matches();

    for filename in values_t!(args, "FILE", String).unwrap().iter() {
        process_file(filename);
    }
}
