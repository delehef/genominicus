use std::fs::File;
use std::io::prelude::*;

use crate::utils::*;
use newick::*;
use svarog::*;

fn draw_tree(
    svg: &mut SvgDrawing,
    t: &Tree,
    n: usize,
    xoffset: f32,
    yoffset: f32,
    render: &RenderSettings,
) -> f32 {
    let mut y = yoffset;
    let leaves_count = t.leaves_of(n).len() as f32;
    let size = 10. * leaves_count.log10();

    if t.descendants(n).iter().any(|d| t[*d].is_duplication()) {
        if let Some(ref cs) = t[n].children {
            for c in cs {
                let thickness = 2.*t.leaves_of(*c).len() as f32;

                if t.descendants(*c).iter().any(|d| t[*d].is_duplication()) {
                    svg.line()
                        .from_coords(xoffset, yoffset, xoffset, y)
                        .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(1.0));
                    svg.line()
                        .from_coords(xoffset, y, xoffset + 2. * size, y)
                        .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(1.0));
                    y = draw_tree(svg, t, *c, xoffset + 2. * size, y, render);
                    // y += 20.;
                } else if t[n].is_duplication() {
                    svg.line()
                        .from_coords(xoffset, yoffset, xoffset, y)
                        .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(1.0));
                    svg.line()
                        .from_coords(xoffset, y, xoffset + 2.*size, y)
                        .style(|s| s.stroke_color(StyleColor::RGB(0, 0, 0)).stroke_width(1.0));
                    svg.polygon()
                        .from_coords([
                            (xoffset + 2. * size, y),
                            (xoffset + 2. * size + thickness, y - thickness / 2.),
                            (xoffset + 2. * size + thickness, y + thickness / 2.),
                        ])
                        .style(|s| s);
                    y += 5. + thickness;
                }
            }
        }
    } else {
        svg.polygon()
            .from_coords([
                (xoffset + 2. * size, y),
                (xoffset + 2. * size + leaves_count, y - leaves_count / 2. - 1.),
                (xoffset + 2. * size + leaves_count, y + leaves_count / 2. + 1.),
            ])
            .style(|s| s);
    }

    if t[n].is_duplication() {
        let dcs = t[n].data.get("DCS").and_then(|s| s.parse::<f32>().ok());
        let elc = t[n].data.get("ELC").and_then(|s| s.parse::<i32>().ok());
        let ellc = t[n].data.get("ELLC").and_then(|s| s.parse::<i32>().ok());

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
            .style(|s| s.fill_color(StyleColor::Percent(1.0 - dcs, dcs, 0.)));
    }
    y + 5.
}

pub fn render(t: &Tree, out_filename: &str, render: &RenderSettings) {
    let mut svg = SvgDrawing::new();
    draw_tree(&mut svg, t, 0, 80., 80., render);
    svg.auto_fit();
    let mut out = File::create(out_filename).unwrap();
    out.write_all(svg.render_svg().as_bytes()).unwrap();
}
