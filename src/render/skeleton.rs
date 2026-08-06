use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;

const STEP_FORWARD: f32 = 20.;

struct TreeState {
    node: NodeHandle,
    xoffset: f32,
    yoffset: f32,
}
fn draw_tree(
    svg: &mut SvgDrawing,
    t: &NewickTree,
    current: TreeState,
    render: &RenderSettings,
) -> f32 {
    let mut y = current.yoffset;
    let leaves_count = t.leaves_of(current.node).len() as f32;
    let size = 10. * leaves_count.log10();
    let step_forward = STEP_FORWARD
        + if t.is_duplication(current.node) {
            size
        } else {
            0.
        };

    if t.descendants(current.node)
        .iter()
        .any(|&d| t.is_duplication(d))
        || t.is_duplication(current.node)
    {
        for &c in t.children(current.node).unwrap() {
            let thickness = (t.leaves_of(c).len() as f32).sqrt();
            let leaves_count = t.leaves_of(c).len() as f32;

            if t.descendants(c).iter().any(|&d| t.is_duplication(d)) || t.is_duplication(c) {
                svg.line()
                    .from_points([
                        (current.xoffset, current.yoffset),
                        (current.xoffset, y),
                        (current.xoffset + step_forward, y),
                    ])
                    .style(|s| {
                        s.stroke_color(StyleColor::RGB(0, 0, 0))
                            .stroke_width(thickness)
                            .fill_color(None)
                    });
                y = draw_tree(
                    svg,
                    t,
                    TreeState {
                        node: c,
                        xoffset: current.xoffset + step_forward,
                        yoffset: y,
                    },
                    render,
                );
            } else if t.is_duplication(current.node) {
                svg.line()
                    .from_points([
                        (current.xoffset, current.yoffset),
                        (current.xoffset, y),
                        (current.xoffset + step_forward, y),
                    ])
                    .style(|s| {
                        s.stroke_color(StyleColor::RGB(0, 0, 0))
                            .stroke_width(thickness)
                            .fill_color(None)
                    });
                svg.polygon()
                    .from_coords([
                        (current.xoffset + step_forward, y),
                        (current.xoffset + step_forward + leaves_count, y),
                        (
                            current.xoffset + step_forward + leaves_count,
                            y + leaves_count,
                        ),
                    ])
                    .style(|s| s);
                y += 5. + leaves_count;
            }
        }
    } else {
        svg.polygon()
            .from_coords([
                (current.xoffset + 2. * size, y),
                (current.xoffset + 2. * size + leaves_count, y),
                (current.xoffset + 2. * size + leaves_count, y + leaves_count),
            ])
            .style(|s| s);
        y += leaves_count;
    }

    if t.is_duplication(current.node) {
        let dcs = t
            .attrs(current.node)
            .get("DCS")
            .and_then(|s| s.parse::<f32>().ok());

        for (label_offset, annotation) in render.node_annotations.iter().enumerate() {
            let label_offset = label_offset as f32;
            if let Some(annotation) = t.attrs(current.node).get(annotation) {
                svg.text()
                    .pos(
                        current.xoffset - FONT_SIZE,
                        current.yoffset + FONT_SIZE + 1.1 * label_offset * FONT_SIZE,
                    )
                    .text(annotation);
            }
        }

        let dcs = dcs.unwrap_or_default();
        svg.polygon()
            .from_pos_dims(
                current.xoffset - size / 2.,
                current.yoffset - size / 2.,
                size,
                size,
            )
            .style(|s| s.fill_color(Some(StyleColor::Percent(1.0 - dcs, dcs, 0.))));
        if render.inner_tags {
            t.attrs(current.node).get("S").map(|name| {
                svg.text()
                    .pos(current.xoffset, current.yoffset - FONT_SIZE)
                    .transform(|t| {
                        t.rotate_from(-30., current.xoffset, current.yoffset - FONT_SIZE)
                    })
                    .text(name)
            });
        }
        y + size
    } else {
        y + 5.
    }
}

pub fn render(t: &NewickTree, out_filename: &str, render: &RenderSettings) {
    let mut svg = SvgDrawing::new();
    draw_tree(
        &mut svg,
        t,
        TreeState {
            node: t.root(),
            xoffset: 80.,
            yoffset: 80.,
        },
        render,
    );
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
