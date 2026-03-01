# Adaptive output cadence parameters

SNEC now derives its profile dumps, scalar diagnostics, and radioactive heating
updates from the adaptive scheduler implemented in
[`src/adaptive_timing.F90`](../src/adaptive_timing.F90). The scheduler reads the
`parameters` file and maps the current simulation time onto "fast", "mid", or
"late" cadence bands before the main loop advances the solution in
[`src/snec.F90`](../src/snec.F90). Each knob in the timing block of
`parameters` has a distinct role:

| Parameter | Meaning |
|-----------|---------|
| `dtout_fast`, `dtout_mid`, `dtout` | Target spacing between full profile dumps (density/temperature profiles saved by `output_all(0:1)`). They correspond to the early/transition/late cadences used by `compute_profile_period`. |
| `dtout_scalar_fast`, `dtout_scalar_mid`, `dtout_scalar` | Target spacing between scalar diagnostics (bolometric and filtered light curves via `output_all(2)` as well as conservation ledgers). They feed `compute_scalar_period`. |
| `output_mid_transition_days`, `output_late_transition_days` | Simulation-time cutoffs (expressed in days) that determine when the scheduler switches from fast → mid → late cadences. Before `output_mid_transition_days` the "fast" targets are used. Between the mid and late thresholds the code uses the "mid" cadence, and after the late threshold it falls back to the slowest targets. |
| `dtout_check` | Fixed interval for `tdump_check`, which forces periodic diagnostics that are not tied to the adaptive scheduler (e.g., extra sanity outputs). Because it bypasses `compute_adaptive_periods`, it remains constant throughout the run. |

The scheduler first picks a target period from the appropriate column in the
table above and then **quantizes** that period so it is an integer multiple of
the current hydrodynamic time step. This behavior is handled by
`quantize_to_step` so the code never schedules an output between two hydro
steps. As the simulation clock advances past the configured transition times,
the cadence automatically migrates from the fast to mid to late targets without
oversampling or skipping physical evolution.

### Practical interpretation of the default values

Using the defaults

```
dtout_fast = 1.728d4        ≈ 0.2 day
output_mid_transition_days = 1.0d0
output_late_transition_days = 1.0d1
```

the solver emits roughly five profile dumps per day throughout the first day of
evolution. Between day 1 and day 10 the cadence relaxes to `dtout_mid = 1.7d5`
(about one dump every two days). Finally, after day 10 the "late" cadence
(`dtout = 1.7d6`, ≈20 days) kicks in. The scalar cadence behaves the same way
but with faster defaults so ugriz magnitudes follow the light curve more closely
at early times. `dtout_check = 1.7d9` (~20000 days) is intentionally long, so
its watchdog output only fires a few times in extremely long runs.
