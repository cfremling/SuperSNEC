# SuperSNEC

SuperSNEC is a fast 1D Eulerian radiation-hydrodynamics code for producing stripped-envelope supernova light curves. It extends the SuperNova Explosion Code ([SNEC 1.01](https://stellarcollapse.org/SNEC)) with runtime adaptive gridding, adaptive cadence control for radioactive heating, optimized gamma-ray deposition, and flexible nickel mixing. A 100-zone SuperSNEC model achieves an RMS magnitude residual of 0.022 mag relative to a 1000-zone SNEC reference, at a runtime of ~1.6 s per model on an Apple M1 Pro (~410x speedup).

## Attribution

SuperSNEC is derived from SNEC 1.01 (Morozova, Piro, & Ott; [stellarcollapse.org/SNEC](https://stellarcollapse.org/SNEC)). This repository contains substantial modifications and additions relative to the original code. The original SNEC is licensed under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).

Original SNEC reference:
> Morozova, V., Piro, A. L., & Ott, C. D. 2015, "Light curves of core-collapse supernovae with substantial mass loss using the new open-source SuperNova Explosion Code (SNEC)", *The Astrophysical Journal*, 814, 63. [doi:10.1088/0004-637X/814/1/63](https://doi.org/10.1088/0004-637X/814/1/63)

## License

This work is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). See [LICENSE](LICENSE) for the full legal text.

## Key Changes from SNEC 1.01

As required by the CC BY-NC-SA 4.0 license, the following summarizes the modifications made in SuperSNEC:

- **Runtime adaptive gridding** --- surface-concentrated initial grid with photosphere-tracking rezone during the simulation
- **Adaptive Ni heating cadence** --- fractional-change trigger with configurable period bounds to skip unnecessary recomputations
- **Optimized gamma-ray deposition** --- precomputed opacities and cache-efficient ray tracing
- **Tuned solver convergence thresholds** --- optimized EPSTOL and ITMAX for stripped-envelope models
- **Flexible Ni mixing** --- box, exponential, and two-component convolution kernels with independent smoothing controls
- **Smooth photosphere luminosity correction** --- eliminates artifacts from discrete zone crossings during photosphere recession at low zone counts
- **Adaptive output cadence** --- three-band (fast/mid/late) output schedule with independent profile and scalar cadences
- **Minimal output mode** --- for production grid runs

## Requirements

- `gfortran` (tested with GCC 13+)
- macOS Accelerate framework or LAPACK/BLAS
- Python 3.10+ with `numpy`, `matplotlib`, `scipy` (for the analysis pipeline)

## Quick Start

```bash
make clean && make
./snec                    # uses ./parameters
# Output written to Data/
```

## Repository Structure

```
src/            Fortran 90 source code
profiles/       Progenitor stellar models (MESA)
tables/         OPAL opacity, bolometric correction, grid pattern
Analysis/       Python analysis pipeline (sweeps, tables, figures)
Paper/          LaTeX source and auto-generated tables
plotting/       Gnuplot scripts for quick visualization of run output
benchmark/      Multi-core throughput benchmarking tools
docs/           Technical notes (adaptive grid, output cadence)
tmp/            Reference light curves for validation
```

## Plotting

The `plotting/` directory contains gnuplot scripts for quick inspection of SNEC output. Run from the repo root:

```bash
./plot.sh       # opens bolometric LC, band magnitudes, and Ni profile plots
```

## Benchmarking

The `benchmark/` directory contains tools for measuring parallel throughput:

```bash
benchmark/benchmark_snec.sh [MAX_JOBS]
```

Generates isolated run copies in `benchmark/runs/` (gitignored). Typical result: ~9 runs/s at 10+ parallel jobs with fast settings (imax=60).

## How to Cite

If you use SuperSNEC, please cite both:

1. Fremling, Christoffer & Hinds, K-Ryan (2026), "SuperSNEC: Fast and Accurate Light Curve Production for Large Hydrodynamic Model Grids Using Adaptive Gridding", submitted to *The Astrophysical Journal*
2. Morozova, V., Piro, A. L., & Ott, C. D. (2015), *The Astrophysical Journal*, 814, 63 --- the original SNEC

## Configuration

All runtime behavior is controlled by the `./parameters` file. See the paper appendix for the full parameter table. Key modes:

- `grid_mode = "adaptive_runtime"` --- SuperSNEC default (adaptive gridding)
- `grid_mode = "legacy_pattern"` --- original SNEC behavior (fixed grid from `tables/GridPattern.dat`)
