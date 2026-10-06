# 04 — Opt-in & Experimental Add-ins

Everything here is **off or neutral by default**. None of it is part of the
validated canonical run (Hasegawa A at the v0.7.5 knobs). These are research
features and diagnostic switches — knowing they exist (and that they're *not*
load-bearing) is the point of this chapter. The README's feature table is the
one-line index; this is the mechanism-level detail.

**Golden rule:** if you're learning the core or debugging a normal run, leave
all of these at their defaults. Flip **one at a time** and re-run the suite
when experimenting.

Most of these exist because of one hard, still-open problem: the code
**over-predicts the ignition pressure spike on high-L/D motors**. That's a
documented, accepted limitation of applying Ma's *quasi-steady* erosive model
during a fast ignition transient — see
[`../v0_7_4/IGNITION_SPIKE_CLOSEOUT.md`](../v0_7_4/IGNITION_SPIKE_CLOSEOUT.md)
and `project_ignition_model_audit`. Groups 1–2 below are the levers that were
tried against it.

---

## Group 1 — Ignition-transient shaping

### Phase F — flame-spread front  `flame_front_enabled` (False), `flame_front_velocity`
[`_advance_flame_front`](../../srm_1d/simulation.py#L794). Normally the bore
fills with hot gas and heats the *whole* grain nearly simultaneously, so many
cells ignite at once → a large synchronized burn onset → pressure spike. Phase F
instead makes ignition **spread as a physical front**: a cell is only
`ignitable` if the front (advancing at the fixed `flame_front_velocity`) has
reached it, or it's in the igniter/cartridge induction zone. The caller
**withholds all surface heating** (convective + pyrogen + radiation) from cells
the front hasn't reached, so they can't light early.

Key design point: `flame_front_velocity` is a **velocity** (~1–10 m/s for
AP/HTPB lateral flame spread; literature default ≈3 m/s), held constant across
motors and **grid-independent**. An earlier formulation derived the speed from
a `q''/(ρ·Cps·ΔT)` flux (that's a *regression* velocity ~mm/s — the wrong
quantity) and was abandoned. Cut the spike ~35 % alone but doesn't fully close
it.

### Phase Z — Zeldovich–Novozhilov dynamic burn rate  `zn_enabled` (False), `kappa_zn`
[`_advance_zn_burn_rate`](../../srm_1d/simulation.py#L846). The core burn rate
is **quasi-steady** — it responds instantly to pressure/flow. Real propellant
has thermal inertia in the solid: the burn rate *lags* changes. Phase Z relaxes
a **dynamic** rate `r_dyn` toward the quasi-steady Ma target with a first-order
lag
```
τ = kappa_zn · alpha_solid / max(r_dyn, floor)²
dr_dyn/dt = (r_qs − r_dyn) / τ
```
(`alpha_solid` = solid thermal diffusivity; `τ` ≈ few ms). It runs every step
(dt-accurate) and overwrites `r_total`/`r_erosive` so both the mass-source and
the regression see the lagged rate.

> **Important discipline:** Phase Z *lags* Ma's output; it never replaces the
> Ma model. `r_qs` (the quasi-steady target) is always the Ma-2020 rate from
> `02`. Don't refactor this into a substitute burn law
> (`feedback_keep_ma_erosive_model`).

Also cut the spike ~35 % alone, but **F and Z don't stack** (F+Z ≈ F) — both
address the same synchronization, so a residual ~1.25× overshoot remains. This
is why the spike is documented as a limitation rather than "fixed."

### Burn establishment lag  `tau_establishment` (0.0)
After a cell ignites, ramp its combustion output linearly 0→1 over
`tau_establishment` seconds instead of switching on at full rate. A softer
onset. Default 0 = instant.

### `port_mach_cap` (0.0 = off)
[`solver.py` Step 3c](../../srm_1d/solver.py#L574). Clips interior face
velocities to `≤ Mach·a`. During the fill, the pressure-based PISO (no interior
choking limit, not shock-capturing) can produce grid-divergent supersonic
contact velocities; this bounds them. Verified **decoupled from `P_peak`** (it
does not change the pressure result — it's solver hygiene, not a spike fix).
See [`../v0_7_4/IGNITION_SPIKE_REOPENED.md`](../v0_7_4/IGNITION_SPIKE_REOPENED.md).

### (Historical) Phase B-v2 h_c augment  — shipped DISABLED
[`_compute_flame_front_augment`](../../srm_1d/simulation.py#L742). An earlier
attempt to boost `h_c` on the cell just downstream of a fresh ignition. It
**double-counted** with PISO's local-Re tracking and *amplified* the spike, so
it ships off. Left in the tree as a documented negative result; don't wire it
back in without reading `docs/v0_7_2/`.

---

## Group 2 — Extra heat paths to the grain surface (ignition assist)

Beyond plain gas convection (`h_c`, `03`), two optional paths can heat a cold
surface toward `T_ignition`:

### DeMar pyrogen surface heating  (heat-delivery mode)
The igniter plume doesn't just add gas — it can **directly heat** the grain
near the igniter (a DeMar-style surface flux). Controlled by the pyrogen's
`heat_flux_cal_cm2_s` and the heat-delivery mode; the plenum's radiated share
is debited from the gas so energy isn't double-counted (`03` Part C). Isolate
with `diagnostic_disable_pyrogen_surface_heating`.

### Adjacent-cell radiation  `propellant.radiation_emissivity` (0.0)
When positive, a **burning cell radiates** to its unignited neighbors, helping
them reach `T_ignition` (models flame radiation spreading ignition). The
emitting gas cell is debited for the radiated energy (unless
`diagnostic_disable_radiation_gas_sink`). At `emissivity = 0` (default) this is
entirely inert.

---

## Group 3 — Igniter topology  `injection_topology` ('forward_plenum')

Where/how the pyrogen enters the bore (`03` Part A):

- **`forward_plenum`** (default) — 0-D plenum vents through an orifice into cell
  0; the validated path.
- **`head_basket` / `aft_basket`** — **uncontained** loose pellets in the bore
  that burn at the **local** bore pressure (no plenum, no orifice, no momentum
  jet), injecting into their host cells. For basket-igniter studies. These
  needed extra heat-flux completeness work to ignite at all (they start at
  atmospheric bore P); see `docs/v0_7_3/`. `aft_basket` in particular is
  geometry-sensitive (cartridge too near the nozzle can stall). Not for
  production runs.

---

## Group 4 — Diagnostic isolation switches (all `False`)

Six `diagnostic_disable_*` flags each turn **off** one physics term so you can
measure its contribution by A/B. They are experiment tools, **not** run
options — a "disabled" run is deliberately non-physical.

| Flag | Turns off |
|---|---|
| `diagnostic_disable_erosive` | the Ma erosive increment (leaves Saint-Robert `r₀`) — the cleanest way to see erosion's effect |
| `diagnostic_disable_endfaces` | end-face regression + end-face mass source |
| `diagnostic_disable_momentum` | pyrogen axial momentum injection (keeps mass + enthalpy) |
| `diagnostic_disable_pyrogen_surface_heating` | DeMar direct surface heating (keeps mass/enthalpy/momentum) |
| `diagnostic_disable_adjacent_radiation` | adjacent-cell radiation (keeps the emissivity in the summary) |
| `diagnostic_disable_radiation_gas_sink` | the gas debit for radiation (keeps the receiver heating) — pure isolation |

`diagnostic_history_capacity` and `initial_gas_temperature` are similar
diagnostic-only overrides (early-terminate a probe run; override the initial
bore-gas temperature) — they don't change the equations.

`history_memmap_path` is an opt-in output-storage override for long or tightly
stepped runs. It keeps every scalar-history row at full float64 resolution but
backs the 42 preallocated channels with one caller-owned disk file instead of
RAM. The default remains in-memory storage. The target file must not already
exist, persists after the run, and must remain available while returned mapped
arrays are in use. Axial snapshots remain in memory and are still controlled by
`snapshot_interval`.

The numerical-verification single-point CLI creates and manifests this file in
its timestamped artifact directory when passed `--history-memmap`; use an
explicit `--history-capacity` when the required row count is already bounded.

---

## Also technically optional: throat erosion
`nozzle.erosion_coeff` / `slag_coeff` (Step 3b). If nonzero, the throat
diameter evolves with the pressure history (erosion widens it, slag narrows
it), recomputed each step. Zero = fixed throat. This is normal SRM physics
(you know it), just flagged here because it's off unless the nozzle provides
the coefficients. Note the vocabulary trap: this **throat-material** erosion is
unrelated to **erosive burning** (`02`) — different physics, no shared
parameter (`feedback_erosion_vocabulary_discipline`).

Next: `05_IO_AND_OPENMOTOR.md` — how motors get in (`.ric`), how results get
out (CSV/plots), and the openMotor plugin boundary.
