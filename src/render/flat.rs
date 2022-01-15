use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;

fn draw_background(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;

    let mut children = node
        .children
        .as_ref()
        .map(|children| children.iter().map(|i| &tree[*i]).collect::<Vec<_>>())
        .unwrap_or_default();
    children.sort_by_cached_key(|x| {
        if let Some(DbGene { species, .. }) = x
            .name
            .as_ref()
            .and_then(|name| name.split('#').next())
            .and_then(|gene_name| genes.get(&gene_name.to_string()))
        {
            species.as_str()
        } else {
            "zzz"
        }
    });

    for &child in children.iter() {
        let new_y = if child.is_leaf() {
            y + 20.
        } else {
            draw_background(
                svg,
                genes,
                depth,
                tree,
                child,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
            )
        };

        if node.is_duplication() {
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
fn draw_gene<'a>(
    svg: &'a mut SvgDrawing,
    x: f32,
    y: f32,
    right: bool,
    color: &StyleColor,
    name: &str,
) -> &'a mut Polygon {
    if right {
        svg.polygon()
            .add_point(x, y - 5.)
            .add_point(x + GENE_WIDTH - 3., y - 5.)
            .add_point(x + GENE_WIDTH, y)
            .add_point(x + GENE_WIDTH - 3., y + 5.)
            .add_point(x, y + 5.)
            .set_hover(name)
            .style(|s| {
                s.fill_color(color.clone())
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
            .set_hover(name)
            .style(|s| {
                s.fill_color(color.clone())
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            })
    }
}

fn draw_tree(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    colormap: &ColorMap,
    depth: f32,
    tree: &Tree,
    node: &Node,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
    links: &mut Vec<(Vec<String>, String, Vec<String>)>,
) -> f32 {
    let mut y = yoffset;
    let mut old_y = 0.;
    let mut children = node
        .children
        .as_ref()
        .map(|children| children.iter().map(|i| &tree[*i]).collect::<Vec<_>>())
        .unwrap_or_default();
    children.sort_by_key(|c| c.name.as_deref().unwrap_or("Z"));
    for (i, child) in children.iter().enumerate() {
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
                .from_points([
                    (xlabels - 5., y),
                    (
                        xlabels + (GENE_WIDTH + GENE_SPACING) * (2. * WINDOW as f32 + 1.)
                            - GENE_SPACING
                            + 5.,
                        y,
                    ),
                ])
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));

            if let Some(name) = &child.name {
                let gene_name = name.split('#').next().unwrap();
                if let Some(DbGene {
                    ancestral,
                    species,
                    chr,
                    left_tail,
                    right_tail,
                    ..
                }) = genes.get(gene_name)
                {
                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {}/{}", gene_name, species, chr))
                        .style(|s| s.fill_color(name2color(species)));

                    // Left tail
                    let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in left_tail.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            xstart,
                            y,
                            true,
                            colormap
                                .get(&g.clone())
                                .unwrap_or(&StyleColor::String("#aaa".to_string())),
                            g,
                        );
                        if g == ancestral {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                        y,
                        true,
                        &gene2color(&ancestral),
                        ancestral,
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, g) in right_tail.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            xstart,
                            y,
                            true, // g.1.to_string() == direction,
                            colormap
                                .get(&g.clone())
                                .unwrap_or(&StyleColor::String("#aaa".to_string())),
                            g,
                        );
                        if g == ancestral {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }
                    links.push((left_tail.to_vec(), ancestral.into(), right_tail.to_vec()));
                } else {
                    // The node was not found in the database
                    eprintln!("{} -- {} not found", name, gene_name);
                    links.push((Vec::new(), name.into(), Vec::new()));
                }
                y += 20.;
            }
        } else {
            svg.line()
                .from_coords(xoffset, y, xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_tree(
                svg,
                genes,
                colormap,
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
        let dcs = node.data.get("DCS").and_then(|s| s.parse::<f32>().ok());

        svg.text()
            .pos(xoffset - FONT_SIZE, yoffset - FONT_SIZE)
            .text(
                dcs.map(|s| format!("{}%", (s * 100.) as u8))
                    .unwrap_or("?".to_string()),
            );
        let dcs = dcs.unwrap_or(0.0);
        svg.polygon()
            .from_pos_dims(xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(1.0 - dcs, dcs, 0.)));
    } else if !node.is_leaf() {
        svg.polygon()
            .from_pos_dims(xoffset - 3., yoffset - 3., 6., 6.)
            .style(|s| s.fill_color(StyleColor::Percent(0., 0., 0.)));
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
                    .from_points([(x1, y + 5.), (x2, y + 20. - 5.)])
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
                    .from_points([(x1, y + 5.), (x2, y + 20. - 5.)])
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

pub fn render(t: &Tree, genes: &GeneCache, colormap: &ColorMap, out_filename: &str) {
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
    draw_background(
        &mut svg, genes, depth, t, &t[0], 10.0, 50.0, xlabels, width,
    );
    let mut links = Vec::new();
    draw_tree(
        &mut svg, genes, colormap, depth, t, &t[0], 10.0, 50.0, xlabels, width, &mut links,
    );
    draw_links(&mut svg, &links, 50.0, xlabels);
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
