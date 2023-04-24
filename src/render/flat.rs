use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;

const MARGIN_TOP: f32 = 100.0;

fn draw_background(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    depth: f32,
    tree: &NewickTree,
    node: usize,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
) -> f32 {
    let mut y = yoffset;

    let mut children = tree.children(node).to_vec();
    children.sort_by_key(|c| tree.name(*c).cloned().unwrap_or_else(|| "Z".to_string()));

    if children.is_empty() {
        return y + 20.;
    }

    for &child in children.iter() {
        let new_y = if tree[child].is_leaf() {
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

        if tree.is_duplication(node) {
            let d = xoffset / depth;
            svg.polygon()
                .from_pos_dims(
                    xoffset + BRANCH_WIDTH / 2.,
                    y - 6.,
                    width - xoffset - d * BRANCH_WIDTH,
                    new_y - y - 6.,
                )
                .style(|s| {
                    s.fill_color(Some(StyleColor::Percent(0.5, 0.5, 1.)))
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
                s.fill_color(Some(color.clone()))
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
                s.fill_color(Some(color.clone()))
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
    tree: &NewickTree,
    n: usize,
    xoffset: f32,
    yoffset: f32,
    xlabels: f32,
    width: f32,
    links: &mut Vec<(f32, Vec<String>, String, Vec<String>)>,
    render: &RenderSettings,
) -> f32 {
    let mut y = yoffset;
    let mut old_y = 0.;
    let mut children = tree[n].children().to_vec();
    children.sort_by_key(|c| tree.name(*c).cloned().unwrap_or_else(|| "Z".to_string()));
    if children.is_empty() {
        return y + 20.;
    }

    for (i, child) in children.iter().enumerate() {
        if i > 0 {
            svg.line()
                .from_coords(xoffset, old_y, xoffset, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        }
        old_y = y;

        if tree[*child].is_leaf() {
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

            if let Some(gene_name) = tree.name(*child).as_ref() {
                if let Some(DbGene {
                    ancestral,
                    species,
                    chr,
                    left_tail,
                    right_tail,
                    ..
                }) = genes.get(gene_name.as_str())
                {
                    // Gene/protein name
                    svg.text()
                        .pos(depth, y + 5.)
                        .text(format!("{} {}/{}", gene_name, species, chr))
                        .style(|s| s.fill_color(Some(name2color(species))));

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
                    links.push((y, left_tail.to_vec(), ancestral.into(), right_tail.to_vec()));
                } else {
                    // The node was not found in the database
                    eprintln!("{} not found", gene_name);
                    links.push((y, Vec::new(), gene_name.to_string(), Vec::new()));
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
                *child,
                xoffset + BRANCH_WIDTH,
                y,
                xlabels,
                width,
                links,
                render,
            );
        }
    }

    let grafting_method = tree.attrs(n).get("METHOD").cloned().unwrap_or_default();
    fn caret<'a>(
        svg: &'a mut SvgDrawing,
        xoffset: f32,
        yoffset: f32,
        w: f32,
        dcs: Option<f32>,
        method: &'a str,
    ) {
        match method {
            "ELC" => {
                let _ = svg.circle().x(xoffset).y(yoffset).radius(w).style(|s| {
                    s.fill_color(Some(if let Some(dcs) = dcs {
                        StyleColor::Percent(1.0 - dcs, dcs, 0.)
                    } else {
                        StyleColor::Percent(0., 0., 0.)
                    }))
                });
            }
            "SEQ" => {
                let _ = svg
                    .polygon()
                    .from_pos_dims(xoffset, yoffset, w, w)
                    .transform(|c| c.rotate(45.))
                    .style(|s| {
                        s.fill_color(Some(if let Some(dcs) = dcs {
                            StyleColor::Percent(1.0 - dcs, dcs, 0.)
                        } else {
                            StyleColor::Percent(0., 0., 0.)
                        }))
                    });
            }
            "SYN" => {
                let _ = svg
                    .polygon()
                    .from_pos_dims(xoffset, yoffset, w, w)
                    .style(|s| {
                        s.fill_color(Some(if let Some(dcs) = dcs {
                            StyleColor::Percent(1.0 - dcs, dcs, 0.)
                        } else {
                            StyleColor::Percent(0., 0., 0.)
                        }))
                    });
            }
            _ => {}
        };
    }

    if tree.is_duplication(n) {
        let dcs = tree.attrs(n).get("DCS").and_then(|s| s.parse::<f32>().ok());
        let elc = tree.attrs(n).get("ELC").and_then(|s| s.parse::<i32>().ok());
        let ellc = tree
            .attrs(n)
            .get("ELLC")
            .and_then(|s| s.parse::<i32>().ok());

        let pretty_dcs = dcs
            .map(|s| format!("{:2.0}%", (s * 100.)))
            .unwrap_or_else(|| "?".to_string());
        let pretty_elc = elc
            .map(|elc| {
                if elc == 0 {
                    "".to_owned()
                } else {
                    format!("L:{}", elc)
                }
            })
            .unwrap_or_else(|| "?".to_string());
        let pretty_ellc = ellc
            .map(|ellc| {
                if ellc == 0 {
                    "".to_owned()
                } else {
                    format!("L:{}", ellc)
                }
            })
            .unwrap_or_else(|| "?".to_string());

        let mut label_offset = 0.;
        for (doit, text) in [
            (render.cs, pretty_dcs),
            (render.elc, pretty_elc),
            (render.ellc, pretty_ellc),
        ]
        .iter()
        {
            if *doit {
                svg.text()
                    .pos(
                        xoffset - FONT_SIZE,
                        yoffset + FONT_SIZE + 1.1 * label_offset * FONT_SIZE,
                    )
                    .text(text);
                label_offset += 1.;
            }
        }

        caret(svg, xoffset - 3., yoffset - 3., 6., dcs, &grafting_method);
    } else if !tree[n].is_leaf() {
        caret(svg, xoffset - 3., yoffset - 3., 6., None, &grafting_method);
    }

    if render.inner_nodes {
        tree.attrs(n).get("S").map(|name| {
            svg.text()
                .pos(xoffset, yoffset - FONT_SIZE)
                .transform(|t| t.rotate_from(-30., xoffset, yoffset - FONT_SIZE))
                .text(name)
        });
    }

    y
}

fn draw_links(
    svg: &mut SvgDrawing,
    links: &[(f32, Vec<String>, String, Vec<String>)],
    xlabels: f32,
) {
    for w in links.windows(2) {
        let xbase = xlabels + (WINDOW as f32 - 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].1.iter().enumerate() {
            let x1 = xbase - i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].1
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase - j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .from_points([(x1, w[0].0 + 5.), (x2, w[1].0 - 5.)])
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }

        let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
        for (i, ancestral) in w[0].3.iter().enumerate() {
            let x1 = xbase + i as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
            for j in
                w[1].3
                    .iter()
                    .enumerate()
                    .filter_map(|(j, name)| if name == ancestral { Some(j) } else { None })
            {
                let x2 = xbase + j as f32 * (GENE_WIDTH + GENE_SPACING) + GENE_WIDTH / 2.;
                svg.line()
                    .from_points([(x1, w[0].0 + 5.), (x2, w[1].0 - 5.)])
                    .style(|s| {
                        s.stroke_color(StyleColor::String("#000".into()))
                            .stroke_width(1.0)
                            .dashed(&[2, 2])
                    });
            }
        }
    }
}

pub fn render(
    t: &NewickTree,
    genes: &GeneCache,
    colormap: &ColorMap,
    out_filename: &str,
    render: &RenderSettings,
) {
    let depth = BRANCH_WIDTH * (t.topological_depth().1 as f32 + 1.);
    let longest_name = (t.leaf_names().map(|name| name.len()).max().unwrap() as f32
        + t.leaves()
            .map(|l| t.attrs(l).get("S").map(|s| s.len()).unwrap_or(0))
            .max()
            .unwrap() as f32
        + 20.)
        * FONT_SIZE;
    let xlabels = 0.85 * (10. + depth + longest_name + 20.);
    let width = xlabels + (2. * WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING) + 60.;
    let mut svg = SvgDrawing::new();
    draw_background(
        &mut svg,
        genes,
        depth,
        t,
        t.root(),
        10.0,
        MARGIN_TOP,
        xlabels,
        width,
    );
    let mut links = Vec::new();
    draw_tree(
        &mut svg,
        genes,
        colormap,
        depth,
        t,
        t.root(),
        10.0,
        MARGIN_TOP,
        xlabels,
        width,
        &mut links,
        render,
    );
    if render.links {
        draw_links(&mut svg, &links, xlabels);
    }
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
