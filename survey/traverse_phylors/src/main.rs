use phylo::prelude::*;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("reading tree...");
    let u1k = File::open("/home/emma/Documents/LDD3/Stage/project/survey/ultra1k.nwk")?;
    let reader = BufReader::new(u1k);

    let mut trees = Vec::new(); 

    for line in reader.lines(){
        let line = line?;
        if line.trim().is_empty(){
            continue;
        }

        let tree = PhyloTree::from_newick(line.as_bytes())?;
        trees.push(tree);
    }

    println!("loaded {} trees!", trees.len());

    // quick counts check on the first tree to verify traversal node counts
    if let Some(first) = trees.first() {
        println!(
            "counts (first tree): pre={}, post={}, bfs={}",
            first.preord_ids(first.get_root_id()).count(),
            first.postord_ids(first.get_root_id()).count(),
            first.bfs_ids(first.get_root_id()).count(),
        );
    }

    if trees.is_empty() {
        return Err("no trees were loaded from ultra1k.nwk".into());
    }

    println!("starting benchmarks...");
    let start = Instant::now();

    let preorder_times = benchmark_sys(&trees, traverse_preorder, trees.len());
    let postorder_times = benchmark_sys(&trees, traverse_postorder, trees.len());
    let levelorder_times = benchmark_sys(&trees, traverse_levelorder, trees.len());
    let ultrametric_times = benchmark_sys(&trees, is_ultrametric, trees.len());

    println!("benchmark completed in {:.2} seconds", start.elapsed().as_secs_f64());

    println!("\n=== phylo-rs Benchmark Statistics ===");
    print_summary("preorder", &preorder_times);
    print_summary("postorder", &postorder_times);
    print_summary("levelorder", &levelorder_times);
    print_summary("is_ultrametric", &ultrametric_times);

    // create violin plot similar to other traverse_* scripts
    let all_methods = vec![
        ("preorder", preorder_times),
        ("postorder", postorder_times),
        ("levelorder", levelorder_times),
        ("is_ultrametric", ultrametric_times),
    ];

    if let Err(e) = draw_violin_plot(&all_methods, "traverse_phylors.png") {
        eprintln!("failed to draw violin plot: {}", e);
    } else {
        println!("\nPlot saved to 'traverse_phylors.png'");
    }

    Ok(())
}

fn benchmark_sys<T, F>(trees: &[PhyloTree], f: F, times: usize) -> Vec<f64>
where
    F: Fn(&PhyloTree) -> T,
{
    let measured = times.min(trees.len());
    let mut ttaken = Vec::with_capacity(measured);

    if let Some(first_tree) = trees.first() {
        black_box(f(first_tree));
    }

    for tree in trees.iter().take(measured) {
        let tstart = Instant::now();
        black_box(f(tree));
        ttaken.push(tstart.elapsed().as_secs_f64() * 1_000_000.0);
    }

    ttaken
}

fn traverse_preorder(tree: &PhyloTree) -> usize {
    tree.preord_ids(tree.get_root_id()).count()
}

fn traverse_postorder(tree: &PhyloTree) -> usize {
    tree.postord_ids(tree.get_root_id()).count()
}

fn traverse_levelorder(tree: &PhyloTree) -> usize {
    tree.bfs_ids(tree.get_root_id()).count()
}

fn is_ultrametric(tree: &PhyloTree) -> bool {
    let mut distances = Vec::new();
    collect_leaf_distances(tree, tree.get_root_id(), 0.0, &mut distances);

    if distances.len() <= 1 {
        return true;
    }

    let min_distance = distances
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max_distance = distances
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    (max_distance - min_distance).abs() < 1e-10
}

fn collect_leaf_distances(
    tree: &PhyloTree,
    node_id: TreeNodeID<PhyloTree>,
    distance: f64,
    distances: &mut Vec<f64>,
) {
    let node = tree.get_node(node_id).unwrap();

    if node.is_leaf() {
        distances.push(distance);
        return;
    }

    for child_id in node.get_children() {
        let child = tree.get_node(child_id).unwrap();
        let child_distance = distance + child.get_weight().unwrap_or(0.0) as f64;
        collect_leaf_distances(tree, child_id, child_distance, distances);
    }
}

fn print_summary(label: &str, times: &[f64]) {
    println!("\n{}:", label);
    println!("  Mean time (µs): {:.4}", mean(times));
    println!("  Median time (µs): {:.4}", median(times));
    println!("  Std dev (µs): {:.4}", std_dev(times));
    println!("  Min time (µs): {:.4}", times.iter().copied().fold(f64::INFINITY, f64::min));
    println!("  Max time (µs): {:.4}", times.iter().copied().fold(f64::NEG_INFINITY, f64::max));
}

fn mean(times: &[f64]) -> f64 {
    times.iter().sum::<f64>() / times.len() as f64
}

fn median(times: &[f64]) -> f64 {
    let mut sorted = times.to_vec();
    sorted.sort_by(|left, right| left.partial_cmp(right).unwrap());

    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

fn std_dev(times: &[f64]) -> f64 {
    let average = mean(times);
    let variance = times
        .iter()
        .map(|value| {
            let diff = value - average;
            diff * diff
        })
        .sum::<f64>()
        / times.len() as f64;

    variance.sqrt()
}

fn draw_violin_plot(data: &Vec<(&str, Vec<f64>)>, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    // collect global y range and max density
    let mut global_min = f64::INFINITY;
    let mut global_max = f64::NEG_INFINITY;
    for (_, times) in data.iter() {
        if times.is_empty() { continue; }
        let min = times.iter().copied().fold(f64::INFINITY, f64::min);
        let max = times.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if min < global_min { global_min = min; }
        if max > global_max { global_max = max; }
    }
    if global_min.is_infinite() || global_max.is_infinite() {
        return Err("no data to plot".into());
    }

    let y_min = global_min.max(1.0);
    let y_lower = (y_min / 1.8).max(1.0);

    // KDE grid
    let grid_n = 120usize;
    let mut densities: Vec<Vec<f64>> = Vec::new();
    let mut max_density = 0f64;
    for (_, times) in data.iter() {
        let den = kde(times, grid_n);
        max_density = max_density.max(den.iter().copied().fold(0./0., f64::max).max(0.0));
        densities.push(den);
    }
    if max_density == 0.0 { max_density = 1.0; }

    // prepare backend
    let root = BitMapBackend::new(filename, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .caption("Phylo-rs: Tree traversal time", ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(80)
        .build_cartesian_2d(0f64..(data.len() as f64 + 1.0), (y_lower..global_max).log_scale())?;

    chart.configure_mesh()
        .disable_x_mesh()
        .x_labels(0)
        .y_desc("Time (microseconds)")
        .draw()?;

    // draw violins
    let xs: Vec<f64> = (0..grid_n).map(|i| i as f64 / (grid_n - 1) as f64).collect();
    for (idx, (_name, times)) in data.iter().enumerate() {
        if times.is_empty() { continue; }
        let den = &densities[idx];
        let x_pos = (idx + 1) as f64;
        let max_half_width = 0.35f64; // in x units

        // create points: left side then right side reversed
        let mut points: Vec<(f64, f64)> = Vec::new();
        for (j, &g) in xs.iter().enumerate() {
            let gx = global_min + (global_max - global_min) * g;
            let w = (den[j] / max_density) * max_half_width;
            points.push((x_pos - w, gx));
        }
        for (j, &g) in xs.iter().enumerate().rev() {
            let gx = global_min + (global_max - global_min) * g;
            let w = (den[j] / max_density) * max_half_width;
            points.push((x_pos + w, gx));
        }

        chart.draw_series(std::iter::once(Polygon::new(
            points,
            &RGBAColor(0, 0, 255, 0.3),
        )))?;

        // median line
        let med = median(&times);
        chart.draw_series(std::iter::once(PathElement::new(vec![(x_pos-0.25, med),(x_pos+0.25, med)], &BLACK)))?;

        chart.draw_series(std::iter::once(Text::new(
            data[idx].0.to_string(),
            (x_pos, y_lower),
            ("sans-serif", 18).into_font(),
        )))?;
    }

    Ok(())
}

fn kde(samples: &Vec<f64>, grid_n: usize) -> Vec<f64> {
    let n = samples.len();
    if n == 0 { return vec![0.0; grid_n]; }
    let mean = samples.iter().sum::<f64>() / n as f64;
    let std = (samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
    let h = if std > 0.0 { 1.06 * std * (n as f64).powf(-0.2) } else { 1.0 };
    let min = samples.iter().copied().fold(f64::INFINITY, f64::min);
    let max = samples.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut densities = Vec::with_capacity(grid_n);
    let two_pi_h = (2.0 * std::f64::consts::PI).sqrt() * h * n as f64;
    for i in 0..grid_n {
        let x = min + (max - min) * (i as f64 / (grid_n - 1) as f64);
        let mut sum = 0.0;
        for &s in samples.iter() {
            let u = (x - s) / h;
            sum += (-0.5 * u * u).exp();
        }
        densities.push(sum / two_pi_h.max(1e-12));
    }
    densities
}
