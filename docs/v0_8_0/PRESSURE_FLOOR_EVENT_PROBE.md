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

The next discriminating experiment was a matching 806-cell, CFL 0.15 event
probe, completed below. CFL 0.075 remains required before testing source cadence
or a counterfactual boundary condition.

## CFL sensitivity

The matching CFL 0.15 run is committed in
`docs/v0_8_0/data/chunc_pressure_floor_cfl_probe_806.json`. It reached the
same configured 0.05 s duration and retained every setting except CFL.

| CFL | Floor-active cells | Floor activations | Minimum pressure | Peak pressure | Peak time | Max fill Mach |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.30 | 2 | 110 | 1.000 kPa | 12.415426 MPa | 23.7352 ms | 0.4463842 |
| 0.15 | 0 | 0 | 10.096875 kPa | 12.419761 MPa | 23.6210 ms | 0.4463722 |

The no-floor result at CFL 0.15 means the aft undershoot is timestep-sensitive.
The corresponding changes in peak pressure, peak time, and first-3-ms maximum
Mach are 0.0349%, 0.4836%, and 0.00267%, respectively. This is not a
floor-free convergence claim: CFL 0.075 remains the next required matched run.
