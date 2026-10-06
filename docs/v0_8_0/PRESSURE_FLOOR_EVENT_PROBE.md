# Chunc Pressure-Floor Event Probe

This committed extract makes the first fine-grid pressure-floor event available
to the team without requiring an ignored local artifact directory. The
machine-readable rows are in
`docs/v0_8_0/data/chunc_pressure_floor_event_probe_806.json`.

## Scope

The probe ran the frozen Chunc MTV startup configuration at target 800 / actual
806 cells, CFL 0.3, and `t_max=0.05 s`. It retained the production 1000 Pa
pressure floor and every physics default. Commit `55b19cb` adds diagnostic-only
event recording immediately before that floor is applied.

The ignored source artifact is
`artifacts/verification_chunc_startup_point/2026-10-05T16-32-34_55b19cb/`.
The committed extract records SHA-256 values for its `point.json` and
`numerical_limits_by_cell.csv`, code revision, compatible openMotor revision,
motor checksum, and baseline configuration checksum.

## Observations

The largest event occurs in gap cell 803 at `x=0.8163576 m`, two cells upstream
of the nozzle. It is directly after the geometry changes from a 22.85 mm bore
cell to a 47.63 mm aft-gap cell. Its local burn mass source is effectively zero.

The first PISO pressure correction changes the cell from the prior 1000 Pa
floor to -1161.60 Pa. The second correction is also negative and yields the
pre-clamp pressure -3882.32 Pa. The first and second local continuity right-hand
sides are negative because modeled eastward mass flow exceeds westward flow.

The next gap cell has a smaller maximum deficit (269.34 Pa) and is adjacent to
the choked nozzle state. The largest event is instead observed while the nozzle
is in the subsonic-outflow state. This rules out a simple statement that only a
choked nozzle boundary causes the first extreme event.

## Interpretation limits

This is a local ledger for the event with each cell's largest recorded deficit.
The pressure-correction solves are global tridiagonal systems, so the local row
does not uniquely attribute a correction to one source cell or prove a boundary
condition is wrong. The data support the narrower conclusion that the event is
associated with the aft bore-to-gap expansion, zero local source, and strong
eastward flux divergence. No physics default, floor threshold, or solver
equation was changed.

The next discriminating experiments were matching 806-cell CFL 0.15 and 0.075
probes, completed below. Cadence and broader-grid testing remain required before
a counterfactual boundary-condition experiment.

## CFL sensitivity

The matching CFL runs are committed in
`docs/v0_8_0/data/chunc_pressure_floor_cfl_probe_806.json`. It reached the
same configured 0.05 s duration and retained every setting except CFL.

| CFL | Floor-active cells | Floor activations | Minimum pressure | Peak pressure | Peak time | Max fill Mach |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.30 | 2 | 110 | 1.000 kPa | 12.415426 MPa | 23.7352 ms | 0.4463842 |
| 0.15 | 0 | 0 | 10.096875 kPa | 12.419761 MPa | 23.6210 ms | 0.4463722 |
| 0.075 | 0 | 0 | 12.080966 kPa | 12.421929 MPa | 23.5637 ms | 0.4463780 |

The no-floor result at CFL 0.15 means the aft undershoot is timestep-sensitive.
The corresponding changes in peak pressure, peak time, and first-3-ms maximum
Mach are 0.0349%, 0.4836%, and 0.00267%, respectively. The floor also remains
inactive at CFL 0.075; relative to CFL 0.075, the CFL 0.15 changes are
0.0174%, 0.2429%, and 0.00130%. This is limited evidence that these selected
startup aggregates have stabilized after the floor disappears. It is not a
full-solution convergence claim, a physics validation, or a reproduction of a
research-note result. The next baseline work is cadence and broader-grid tests
from the floor-free CFL 0.075 configuration.

## Floor-free spatial probe

The 406-, 806-, and 1606-cell CFL 0.075 runs are recorded in
`docs/v0_8_0/data/chunc_floor_free_spatial_probe.json`. All three reached the
configured 0.05 s duration with no pressure-floor activation. Burn-rate and
geometry update intervals followed the existing solver default, which scales
as `max(10, actual_cells // 5)`.

| Actual cells | Peak pressure | Peak time | Max fill Mach | Minimum pressure |
| ---: | ---: | ---: | ---: | ---: |
| 406 | 12.377546 MPa | 23.5142 ms | 0.4588228 | 15.5386 kPa |
| 806 | 12.421929 MPa | 23.5637 ms | 0.4463780 | 12.0810 kPa |
| 1606 | 12.440641 MPa | 23.8922 ms | 0.4402303 | 9.9656 kPa |

From 806 to 1606 cells, peak pressure, peak time, and first-3-ms maximum Mach
change by 0.150%, 1.375%, and 1.396%. These pass the provisional 2% screen.
Minimum pressure changes by 21.23%, however, despite remaining above the floor,
and therefore is not spatially converged. The peak-time differences also grow
between the last two refinements, so no formal order or extrapolated value is
claimed. This study remains startup-only numerical evidence, not physics
validation or reproduction of a research-note result.

## Floor-free cadence probe

The complete four-point result is committed as
`docs/v0_8_0/data/chunc_cadence_sensitivity_806.json` (SHA-256
`8ac6a260fd9d77a26bfda46c4fbb19cc87f65e4893c24520d9438355754789e6`).
Every point used 806 actual cells, CFL 0.075, and the same 0.05 s startup
configuration. The only changes were explicit burn-rate/source and geometry
update intervals. All four runs passed the health gate and had zero pressure-
floor activations.

| Burn steps | Geometry steps | Peak pressure | Peak time | Minimum pressure | Max ignition refresh delay |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 161 | 161 | 12.421929 MPa | 23.5637 ms | 12.0810 kPa | 8.827 us |
| 80 | 161 | 12.423295 MPa | 23.5379 ms | 12.4893 kPa | 4.358 us |
| 161 | 80 | 12.422284 MPa | 23.5643 ms | 12.0810 kPa | 8.827 us |
| 80 | 80 | 12.422979 MPa | 23.5386 ms | 12.4893 kPa | 4.358 us |

Halving only the burn/source interval changes peak pressure by 0.0110%, peak
time by 0.109%, clipping-correction energy by 0.532%, and first-3-ms maximum
Mach by 0%. Minimum pressure changes by 3.380%, exceeding the provisional 2%
screen. Halving only the geometry interval changes every reported physical
metric by less than 0.003%. The combined result closely follows the burn/source
refinement, so no material cadence interaction appears in these startup
metrics. The approximately halved ignition refresh delay is the intended
numerical consequence of the finer burn/source interval, not a physical
validation result.

Together with the spatial probe, this supports using the CFL 0.075 default-
cadence case as a bounded startup reference for peak pressure, peak time, and
early Mach. It does not establish convergence of transient minimum pressure,
full-solution convergence, full-burn cadence convergence, or physics validity.
