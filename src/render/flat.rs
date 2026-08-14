use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;
use syntesuite::genebook::FamilyId;
use syntesuite::genebook::Gene;
use syntesuite::Strand;

const MARGIN_TOP: f32 = 100.0;

struct DrawState {
    node: NodeHandle,
    xoffset: f32,
    yoffset: f32,
}

fn draw_background(
    svg: &mut SvgDrawing,
    depth: f32,
    tree: &NewickTree,
    current: DrawState,
    width: f32,
) -> f32 {
    let mut y = current.yoffset;

    let mut children = tree.children(current.node).unwrap().to_vec();
    children.sort_by_key(|c| tree.name(*c).cloned().unwrap_or_else(|| "Z".to_string()));

    if children.is_empty() {
        return y + 20.;
    }

    for &child in children.iter() {
        let new_y = if tree.is_leaf(child) {
            y + 20.
        } else {
            draw_background(
                svg,
                depth,
                tree,
                DrawState {
                    node: child,
                    xoffset: current.xoffset + BRANCH_WIDTH,
                    yoffset: y,
                },
                width,
            )
        };

        if tree.is_duplication(current.node) {
            let d = current.xoffset / depth;
            svg.polygon()
                .from_pos_dims(
                    current.xoffset + BRANCH_WIDTH / 2.,
                    y - 6.,
                    width - current.xoffset - d * BRANCH_WIDTH,
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

struct GeneGlyph<'a> {
    x: f32,
    y: f32,
    strand: Strand,
    color: &'a StyleColor,
    name: &'a str,
}
fn draw_gene<'a>(svg: &'a mut SvgDrawing, g: GeneGlyph<'_>) -> &'a mut Polygon {
    match g.strand {
        Strand::Direct => svg
            .polygon()
            .add_point(g.x, g.y)
            .add_point(g.x + 3., g.y - 5.)
            .add_point(g.x + GENE_WIDTH, g.y - 5.)
            .add_point(g.x + GENE_WIDTH, g.y + 5.)
            .add_point(g.x + 3., g.y + 5.)
            .set_hover(g.name)
            .style(|s| {
                s.fill_color(Some(g.color.clone()))
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            }),
        Strand::Reverse => svg
            .polygon()
            .add_point(g.x, g.y - 5.)
            .add_point(g.x + GENE_WIDTH - 3., g.y - 5.)
            .add_point(g.x + GENE_WIDTH, g.y)
            .add_point(g.x + GENE_WIDTH - 3., g.y + 5.)
            .add_point(g.x, g.y + 5.)
            .set_hover(g.name)
            .style(|s| {
                s.fill_color(Some(g.color.clone()))
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            }),
        Strand::Unknown => svg
            .polygon()
            .add_point(g.x + 1.5, g.y - 5.)
            .add_point(g.x + GENE_WIDTH - 1.5, g.y - 5.)
            .add_point(g.x + GENE_WIDTH - 1.5, g.y + 5.)
            .add_point(g.x + 1.5, g.y + 5.)
            .set_hover(g.name)
            .style(|s| {
                s.fill_color(Some(g.color.clone()))
                    .stroke_width(0.5)
                    .stroke_color(StyleColor::Percent(0.2, 0.2, 0.2))
            }),
    }
}

fn draw_tree(
    svg: &mut SvgDrawing,
    genes: &GeneCache,
    colormap: &ColorMap,
    petmap: &PetnameMap,
    depth: f32,
    tree: &NewickTree,
    current: DrawState,
    xlabels: f32,
    links: &mut Vec<(f32, Vec<FamilyId>, FamilyId, Vec<FamilyId>)>,
    render: &RenderSettings,
) -> f32 {
    let mut y = current.yoffset;
    let mut old_y = 0.;
    let mut children = tree.children(current.node).unwrap().to_vec();
    children.sort_by_key(|c| tree.name(*c).cloned().unwrap_or_else(|| "Z".to_string()));
    if children.is_empty() {
        return y + 20.;
    }

    for (i, child) in children.iter().enumerate() {
        if i > 0 {
            svg.line()
                .from_coords(current.xoffset, old_y, current.xoffset, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
        }
        old_y = y;

        if tree.is_leaf(*child) {
            // Leaf branch
            svg.line()
                .from_coords(current.xoffset, y, depth, y)
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
                if let Some(Gene {
                    family,
                    species,
                    chr,
                    strand,
                    left_landscape,
                    right_landscape,
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
                    for (k, tg) in left_landscape.iter().enumerate() {
                        let xstart = xbase - (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            GeneGlyph {
                                x: xstart,
                                y,
                                strand: tg.strand,
                                color: colormap
                                    .get(&tg.family)
                                    .unwrap_or(&StyleColor::String("#aaa".to_string())),
                                name: &petmap[&tg.family],
                            },
                        );
                        if tg.family == *family {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }

                    // The Gene
                    draw_gene(
                        svg,
                        GeneGlyph {
                            x: xlabels + WINDOW as f32 * (GENE_WIDTH + GENE_SPACING),
                            y,
                            strand: *strand,
                            color: &gene2color(&family.to_ne_bytes()),
                            name: &petmap[family],
                        },
                    )
                    .style(|s| {
                        s.stroke_width(2.)
                            .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                    });

                    // Right tail
                    let xbase = xlabels + (WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING);
                    for (k, tg) in right_landscape.iter().enumerate() {
                        let xstart = xbase + (k as f32) * (GENE_WIDTH + GENE_SPACING);
                        let drawn = draw_gene(
                            svg,
                            GeneGlyph {
                                x: xstart,
                                y,
                                strand: tg.strand,
                                color: colormap
                                    .get(&tg.family)
                                    .unwrap_or(&StyleColor::String("#aaa".to_string())),
                                name: &petmap[&tg.family],
                            },
                        );
                        if tg.family == *family {
                            drawn.style(|s| {
                                s.stroke_width(2.)
                                    .stroke_color(StyleColor::Percent(0.1, 0.1, 0.1))
                            });
                        }
                    }
                    links.push((
                        y,
                        left_landscape
                            .iter()
                            .map(|tg| tg.family)
                            .collect::<Vec<_>>(),
                        *family,
                        right_landscape
                            .iter()
                            .map(|tg| tg.family)
                            .collect::<Vec<_>>(),
                    ));
                } else {
                    // The node was not found in the database
                    eprintln!("{} not found", gene_name);
                    links.push((y, Vec::new(), 0.into(), Vec::new()));
                }
            }
            y += 20.;
        } else {
            svg.line()
                .from_coords(current.xoffset, y, current.xoffset + BRANCH_WIDTH, y)
                .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(0.5));
            y = draw_tree(
                svg,
                genes,
                colormap,
                petmap,
                depth,
                tree,
                DrawState {
                    node: *child,
                    xoffset: current.xoffset + BRANCH_WIDTH,
                    yoffset: y,
                },
                xlabels,
                links,
                render,
            );
        }
    }

    let grafting_method = tree
        .attrs(current.node)
        .get("METHOD")
        .cloned()
        .unwrap_or_default();
    struct CaretGlyph<'a> {
        xoffset: f32,
        yoffset: f32,
        width: f32,
        dcs: Option<f32>,
        method: &'a str,
    }
    fn caret(svg: &mut SvgDrawing, g: CaretGlyph<'_>) {
        match g.method {
            "ELC" => {
                let _ = svg
                    .circle()
                    .x(g.xoffset)
                    .y(g.yoffset)
                    .radius(g.width / 2.)
                    .style(|s| {
                        s.fill_color(Some(if let Some(dcs) = g.dcs {
                            StyleColor::Percent(1.0 - dcs, dcs, 0.)
                        } else {
                            StyleColor::Percent(0., 0., 0.)
                        }))
                    });
            }
            "SEQ" => {
                let _ = svg
                    .polygon()
                    .from_pos_dims(
                        g.xoffset - g.width / 2.,
                        g.yoffset - g.width / 2.,
                        g.width,
                        g.width,
                    )
                    .transform(|c| c.rotate_from(45., g.xoffset, g.yoffset))
                    .style(|s| {
                        s.fill_color(Some(if let Some(dcs) = g.dcs {
                            StyleColor::Percent(1.0 - dcs, dcs, 0.)
                        } else {
                            StyleColor::Percent(0., 0., 0.)
                        }))
                    });
            }
            "SYN" => {
                let _ = svg
                    .polygon()
                    .from_pos_dims(
                        g.xoffset - g.width / 2.,
                        g.yoffset - g.width / 2.,
                        g.width,
                        g.width,
                    )
                    .style(|s| {
                        s.fill_color(Some(if let Some(dcs) = g.dcs {
                            StyleColor::Percent(1.0 - dcs, dcs, 0.)
                        } else {
                            StyleColor::Percent(0., 0., 0.)
                        }))
                    });
            }
            _ => {
                if let Some(dcs) = g.dcs {
                    let _ = svg
                        .polygon()
                        .from_pos_dims(
                            g.xoffset - g.width / 2.,
                            g.yoffset - g.width / 2.,
                            g.width,
                            g.width,
                        )
                        .style(|s| {
                            s.stroke_color(StyleColor::Percent(1.0 - dcs, dcs, 0.))
                                .fill_color(None)
                                .stroke_width(2.)
                        });
                }
            }
        };
    }

    for (label_offset, annotation) in render.node_annotations.iter().enumerate() {
        if let Some(annotation) = tree.attrs(current.node).get(annotation) {
            svg.text()
                .pos(
                    current.xoffset - FONT_SIZE,
                    current.yoffset + FONT_SIZE + 1.1 * label_offset as f32,
                )
                .text(annotation);
        }
    }

    caret(
        svg,
        CaretGlyph {
            xoffset: current.xoffset,
            yoffset: current.yoffset,
            width: 6.,
            dcs: tree[current.node]
                .data()
                .attrs
                .get("DCS")
                .and_then(|dcs| str::parse::<f32>(dcs).ok()),
            method: &grafting_method,
        },
    );

    if render.inner_tags {
        tree.attrs(current.node).get("S").map(|name| {
            svg.text()
                .pos(current.xoffset, current.yoffset - FONT_SIZE)
                .transform(|t| t.rotate_from(-30., current.xoffset, current.yoffset - FONT_SIZE))
                .text(name)
        });
    }

    y
}

fn draw_links(
    svg: &mut SvgDrawing,
    links: &[(f32, Vec<FamilyId>, FamilyId, Vec<FamilyId>)],
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
    petmap: &PetnameMap,
    out_filename: &str,
    render: &RenderSettings,
) {
    let depth = BRANCH_WIDTH * (t.topological_depth().unwrap().1 as f32 + 1.);
    let longest_name = (t.leaf_names().map(|name| name.len()).max().unwrap_or(0) as f32
        + t.leaves()
            .map(|l| t.attrs(l).get("S").map(|s| s.len()).unwrap_or(0))
            .max()
            .unwrap_or(0) as f32
        + 20.)
        * FONT_SIZE;
    let xlabels = 0.85 * (10. + depth + longest_name + 20.);
    let width = xlabels + (2. * WINDOW as f32 + 1.) * (GENE_WIDTH + GENE_SPACING) + 60.;
    let mut svg = SvgDrawing::new();
    draw_background(
        &mut svg,
        depth,
        t,
        DrawState {
            node: t.root(),
            xoffset: 10.0,
            yoffset: MARGIN_TOP,
        },
        width,
    );
    let mut links = Vec::new();
    draw_tree(
        &mut svg,
        genes,
        colormap,
        petmap,
        depth,
        t,
        DrawState {
            node: t.root(),
            xoffset: 10.0,
            yoffset: MARGIN_TOP,
        },
        xlabels,
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
