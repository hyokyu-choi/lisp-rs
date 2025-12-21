use plotters::prelude::*;

pub fn plot_one(
    title: &str,
    ts: Vec<f64>,
    datas: Vec<Vec<f64>>,
    labels: Vec<&str>,
    save_file: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(&save_file, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .caption(title, ("sans-serif", 50).into_font())
        .build_cartesian_2d(0.0f32..10.0f32, -2.0f32..2.0f32)?;

    chart.configure_mesh().draw()?;

    for idx in 0..(datas.len()) {
        let ts_32 = ts.iter().map(|t| *t as f32);
        let data_32 = datas[idx].iter().map(|x| *x as f32);
        chart
            .draw_series(LineSeries::new(
                std::iter::zip(ts_32, data_32),
                &Palette99::pick(idx),
            ))?
            .label(labels[idx])
            .legend(move |(x, y)| {
                PathElement::new(vec![(x, y), (x + 20, y)], &Palette99::pick(idx))
            });
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    Ok(())
}
