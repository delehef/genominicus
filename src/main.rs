use crate::nhx::*;
use std::fs::File;
use std::io::prelude::*;
use svarog::*;

const WINDOW: usize = 15;
const GENE_WIDTH: f32 = 15.;
const GENE_SPACING: f32 = 5.;
const BRANCH_WIDTH: f32 = 20.;
const FONT_SIZE: f32 = 10.;

mod nhx;

fn drawDs(
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
            drawDs(
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
            svg.push(Box::new(
                Polygon::from_pos_dims(
                    xoffset + BRANCH_WIDTH / 2., y - 6.,
                    width - xoffset - d * BRANCH_WIDTH, new_y - y - 6.,
                )
                .style(|s| {
                    s.fill_color(StyleColor::Percent(0.5, 0.5, 1.))
                        .fill_opacity(0.1 + 0.9 * d)
                }),
            ));
        }
        if child.is_leaf() {
            svg.push(Box::new(Text::from_pos(depth, y+5.).text(&child.name.as_ref().unwrap_or(&"".into()))));
        }
        y = new_y;
    }
    y
}

fn main() {
    let t = Tree::from_filename(
        "/home/franklin/work/duplications/done/prims_065/final_trees/76-profilenj.nhx",
    )
    .unwrap();
    t.print();

    let depth = 280.;//BRANCH_WIDTH * (t.topological_depth().1 + 1.);
    let longest_name = t
        .leaf_names()
        .iter()
        .map(|(_, name)| name.map_or(0, |s| s.len()))
        .max()
        .unwrap() as f32
        * FONT_SIZE;
    let xlabels = 0.85 * (10. + depth + longest_name + 50.);
    let width = xlabels + (2. * WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING) + 60.;
    let height = 20. * (t.leaves().count() as f32 + 2.) + 20.;
    dbg!(longest_name);
    dbg!(xlabels);
    dbg!(width);
    dbg!(height);
    dbg!(depth);
    let mut svg = SvgDrawing::with_size(width as usize, height as usize);
    svg.push(Box::new(
        Text::new().pos(FONT_SIZE, FONT_SIZE).text("Title"),
    ));
    drawDs(&mut svg, depth, &t, 0, 10.0, 50.0, xlabels, width);
    let mut out = File::create("out.svg").unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
