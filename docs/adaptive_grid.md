# Runtime adaptive grid

SNEC now supports two grid modes:

1. `adaptive_runtime` (default): start from a surface-concentrated baseline
   grid and remesh during evolution.
2. `legacy_pattern`: reproduce the historical fixed grid from
   `grid_pattern_file` (`tables/GridPattern.dat` by default).

Modes such as `fixed_pattern`, `uniform_in_mass`, and
`from_file_by_mass` are removed from the active code path.

## Adaptive monitor and remeshing

When `grid_mode = "adaptive_runtime"`, the code rebuilds the mass grid every
`grid_update_interval_days` (simulation days). If the interval is `<= 0`, no
runtime remesh is performed and the initial adaptive baseline grid remains
fixed.

The remesh monitor combines:

- optical-depth weighting in outer layers (for breakout/cooling phases),
- Ni-rich core emphasis,
- composition-interface signals (H/He/Ni gradients).

Cell weights are equidistributed to obtain new mass coordinates. A cap can be
applied with `early_bias_max_weight_ratio` to prevent a few zones from taking
extreme weight and collapsing the timestep.

## Conservative remap

After each accepted remesh:

- boundary fields (`r`, `vel`) are interpolated in mass,
- cell-centered scalars and composition are conservatively remapped,
- geometry, density, viscosity terms, EOS, opacity, luminosity, and optical
  depth are rebuilt.

The remap is skipped when the proposed new grid is effectively unchanged, which
avoids unnecessary interpolation diffusion.

## User-facing controls

- `grid_mode`: `adaptive_runtime` or `legacy_pattern`.
- `grid_update_interval_days`: remesh cadence in runtime adaptive mode.
- `grid_pattern_file`: required only for `legacy_pattern`.
- `grid_inner_frac_early`, `grid_inner_frac_late`: early/late inner-region
  fractions used by runtime target blending.
- `grid_outer_ratio_early`, `grid_outer_ratio_late`: early/late geometric
  tapering of outer cells.
- `grid_thin_fraction_threshold`, `grid_thin_fraction_scale`,
  `grid_relax_alpha`: transition control and smoothing.
- `early_bias_max_weight_ratio`: optional cap on monitor weight contrast.
