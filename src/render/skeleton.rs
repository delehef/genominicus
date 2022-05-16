use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;

const STEP_FORWARD: f32 = 20.;

fn draw_tree(
    svg: &mut SvgDrawing,
    t: &NewickTree,
    n: usize,
    xoffset: f32,
    yoffset: f32,
    render: &RenderSettings,
) -> f32 {
    let mut y = yoffset;
    let leaves_count = t.leaves_of(n).len() as f32;
    let size = 10. * leaves_count.log10();
    let step_forward = STEP_FORWARD + if t.is_duplication(n) { size } else { 0. };

    if t.descendants(n).iter().any(|&d| t.is_duplication(d)) || t.is_duplication(n) {
        for &c in t[n].children() {
            let thickness = (t.leaves_of(c).len() as f32).sqrt();
            let leaves_count = t.leaves_of(c).len() as f32;

            if t.descendants(c).iter().any(|&d| t.is_duplication(d)) || t.is_duplication(c) {
                svg.line()
                    .from_points([
                        (xoffset, yoffset),
                        (xoffset, y),
                        (xoffset + step_forward, y),
                    ])
                    .style(|s| {
                        s.stroke_color(StyleColor::RGB(0, 0, 0))
                            .stroke_width(thickness)
                            .fill_color(None)
                    });
                y = draw_tree(svg, t, c, xoffset + step_forward, y, render);
            } else if t.is_duplication(n) {
                svg.line()
                    .from_points([
                        (xoffset, yoffset),
                        (xoffset, y),
                        (xoffset + step_forward, y),
                    ])
                    .style(|s| {
                        s.stroke_color(StyleColor::RGB(0, 0, 0))
                            .stroke_width(thickness)
                            .fill_color(None)
                    });
                svg.polygon()
                    .from_coords([
                        (xoffset + step_forward, y),
                        (xoffset + step_forward + leaves_count, y),
                        (xoffset + step_forward + leaves_count, y + leaves_count),
                    ])
                    .style(|s| s);
                y += 5. + leaves_count;
            }
        }
    } else {
        svg.polygon()
            .from_coords([
                (xoffset + 2. * size, y),
                (xoffset + 2. * size + leaves_count, y),
                (xoffset + 2. * size + leaves_count, y + leaves_count),
            ])
            .style(|s| s);
        y += leaves_count;
    }

    if t.is_duplication(n) {
        let dcs = t[n]
            .data
            .attrs
            .get("DCS")
            .and_then(|s| s.parse::<f32>().ok());
        let elc = t[n]
            .data
            .attrs
            .get("ELC")
            .and_then(|s| s.parse::<i32>().ok());
        let ellc = t[n]
            .data
            .attrs
            .get("ELLC")
            .and_then(|s| s.parse::<i32>().ok());

        let pretty_dcs = dcs
            .map(|s| format!("{:2.0}%", (s * 100.)))
            .unwrap_or("?".to_string());
        let pretty_elc = elc
            .map(|elc| {
                if elc == 0 {
                    "".to_owned()
                } else {
                    format!("L:{}", elc)
                }
            })
            .unwrap_or("?".to_string());
        let pretty_ellc = ellc
            .map(|ellc| {
                if ellc == 0 {
                    "".to_owned()
                } else {
                    format!("L:{}", ellc)
                }
            })
            .unwrap_or("?".to_string());

        let mut label_offset = 0.;
        if render.cs {
            svg.text()
                .pos(
                    xoffset - FONT_SIZE,
                    yoffset + FONT_SIZE + 1.1 * label_offset * FONT_SIZE,
                )
                .text(pretty_dcs);
            label_offset += 1.;
        }
        if render.elc {
            svg.text()
                .pos(
                    xoffset - FONT_SIZE,
                    yoffset + FONT_SIZE + 1.1 * label_offset * FONT_SIZE,
                )
                .text(pretty_elc);
            label_offset += 1.;
        }
        if render.ellc {
            svg.text()
                .pos(
                    xoffset - FONT_SIZE,
                    yoffset + FONT_SIZE + 1.1 * label_offset * FONT_SIZE,
                )
                .text(pretty_ellc);
            label_offset += 1.;
        }
        let dcs = dcs.unwrap_or(0.0);
        svg.polygon()
            .from_pos_dims(xoffset - size / 2., yoffset - size / 2., size, size)
            .style(|s| s.fill_color(Some(StyleColor::Percent(1.0 - dcs, dcs, 0.))));
        if render.inner_nodes {
            t[n].data.name.as_ref().map(|name| {
                svg.text()
                    .pos(xoffset, yoffset - FONT_SIZE)
                    .transform(|t| t.rotate_from(-30., xoffset, yoffset - FONT_SIZE))
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
    draw_tree(&mut svg, t, t.root(), 80., 80., render);
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
