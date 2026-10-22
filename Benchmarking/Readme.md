# Benchmarking Predictive Simulations of Human Gait

This repository provides tools to **benchmark** predictive simulations of human gait against experimental data.  
The goal is to assess how well your model reproduces experimentally observed adaptations under various gait conditions (e.g., walking on slopes, walking with ankle loads, etc.).

For more background and methodology, please refer to our paper:  
[Afschrift et al. (2025), *Benchmarking the predictive capability of human gait simulations*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012713).

---

## Overview

This repository contains MATLAB scripts to:

1. Run predictive simulations for several experimental gait conditions.
2. Compare your simulation results to experimental data. (Automatically download and link benchmark datasets for analysis and visualization.)

Experimental datasets are hosted in a separate repository:  
👉 [**predsim_benchmark_data**](https://github.com/MaartenAfschrift/predsim_benchmark_data)

---

## Getting Started

### 1. Running Benchmark Simulations

Use the MATLAB function **`benchmark_predsim.m`** to execute predictive simulations for multiple gait conditions.

You can extend or modify the benchmark by adding new gait conditions directly in this function.

Example script: **`benchmark_falisse2022.m`**

This example shows how to call the benchmarking function.

#### Required Inputs

- **`S`** – Default settings structure for *PredSim* (see PredSim manual).  
- **`osim_path`** – Path to your OpenSim model (standard PredSim input).  
- **`S_benchmark`** – Structure defining which simulations to run and where to save them.

Example configuration:

```matlab
S_benchmark.studies = {'vanderzee2022', 'browning2008'};
S_benchmark.out_folder = fullfile(pathRepo, 'Results', 'Benchmark_Falisse2022');

S_benchmark.gait_speeds = true;
S_benchmark.gait_speeds_selection = 0.6:0.2:2;
```

Supported studies are `vanderzee2022`, `browning2008`, `koelewijn2019`,
`gomenuka2014`, and `schertzer2014`. An empty or omitted `studies` list runs
only the requested gait-speed sweep. With `gait_speeds = true` and no speed
selection or range, the sweep is 0.6:0.1:2 m/s. The example uses MATLAB's
default parallel cluster; a custom profile is optional.

---

### 2. Comparing Simulations with Experimental Data

After running your simulations, you can compare them to experimental results using data from the [**predsim_benchmark_data**](https://github.com/MaartenAfschrift/predsim_benchmark_data) repository.

You can either:

- Write your own comparison scripts, **or**
- Use the provided helper functions for automatic comparison and visualization.

#### Key Functions

**`add_benchmarkdata_to_simresults.m`**

- Automatically downloads the required experimental datasets from GitHub and stores them in `PredSim/Benchmarking/Data`.
- Loops through all simulations in your results folder (`S_benchmark.out_folder`).
- Finds corresponding experimental data (using unique identifiers).
- Adds a `benchmark` structure to each simulation result file.

To visualize results:

```matlab
add_benchmarkdata_to_simresults(..., 'BoolPlot', true)
```

Cached experimental data is reused. To explicitly refresh it, pass
`'OverwriteData', true`. Exact experimental IDs are preferred. Schertzer
data contains small differences between digitized and prescribed speeds;
matching allows at most 0.005 m/s only when load location and total mass
match and there is a single candidate. Its zero-load torso model is compared
with the torso baseline; ankle and knee baselines are distinct measurements.
Unmatched or ambiguous IDs are reported rather than assigned arbitrarily.

Plot normalization uses the unloaded mass from `S.subject.mass` or the
original OpenSim model. If that model has moved, pass `'SubjectMass', massKg`.
The default normalization leg length is 0.85 m for the Falisse example;
pass `'LegLength', lengthM` for a different subject. For example:

```matlab
add_benchmarkdata_to_simresults(resultsFolder, 'BoolPlot', true, ...
    'SubjectMass', 70, 'LegLength', 0.90);
```

Partial runs can be plotted: missing matches and explicitly unconverged
results are excluded from per-study comparisons. Existing nonempty IDs are
preserved when rerunning a benchmark; conflicting IDs require checking the
original simulation condition before migrating old results. Reattaching
data clears an old experimental attachment if the current ID has no unique
match.

This generates default comparison figures for quick inspection.
Requested coordinates absent from a simulation model produce a warning and
blank curves; the remaining comparisons are still plotted.

The metabolic plots integrate the simulated muscle metabolic power over a
complete gait cycle; they do not add a basal metabolic offset. Interpret
absolute experimental comparisons with care: the
[benchmark paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012713)
primarily compares changes from reference conditions because resting-rate
measurements differ between studies. The quick inspection figures do not
reproduce that complete reference-condition analysis.

---

## Data Structure

Each benchmarked study has:

- A unique identifier linking simulation and experimental data.
- Consistent metadata for gait conditions and measured variables.

More details about the data format and organization are available in the [**predsim_benchmark_data**](https://github.com/MaartenAfschrift/predsim_benchmark_data) repository.

---

## Citation

If you use this code or data in your research, please cite:

> Afschrift, M., et al. (2025). *Benchmarking the predictive capability of human gait simulations.*

---

## Contributing

Contributions and extensions are welcome!  
You can:

- Add new gait conditions to `benchmark_predsim.m`
- Contribute new experimental datasets or analysis tools to [**predsim_benchmark_data**](https://github.com/MaartenAfschrift/predsim_benchmark_data)

Please submit a pull request or open an issue for discussion.

---

## License

This project is released under the [MIT License](LICENSE) (or specify your actual license).
