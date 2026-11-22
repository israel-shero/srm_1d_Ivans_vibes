# Provenance and *-flagged propellant transport data

The machine-readable register is [data/propellant_source_register.yaml](data/propellant_source_register.yaml). It deliberately separates **documented manufacturer values** from **calculated estimates*** and **unverified historical composition hypotheses***.

## Evidence and limitations

1. **RCS 8225W White Lightning and RCS 8125RC Redline data sheets:** Lists typical product-level chamber temperature, gamma, molecular weight and propellant density, among other properties. The sheet warns that its values are not specifications. These are documented manufacturer claims, **not** independently measured values.
2. **AeroTech/RCS SDS (2023):** Identifies broad chemical families and explicitly withholds some identities/concentrations as trade secrets. It does not authenticate a particular recipe or quantify combustion-gas transport.
3. **TQC/TRF repost (2006):** Includes claimed Red/White/Blue formulations with a third-party attribution. It is **not proof** that they correspond to modern RCS/AeroTech Redline, White Lightning, or Blue Thunder products. Their inclusion in this register tracks the historical hypothesis, not the formula itself.
4. **Southwest Research Institute / DOT report (2020):** Identifies Blue Thunder chemistry as proprietary, consistent with treating any external recipe identification as tentative.
5. **Rutter et al., *Aerosol and Air Quality Research* (2007):** Measured zinc- and chloride-rich fine-particle emissions from tested White Lightning, Blue Thunder and Black Jack motors. This adds independent *exhaust* chemistry evidence, but the authors' discussion of *unburned propellant* composition was explicitly an assumption; it is neither a current product-specification assay nor a unique reverse identification of a commercial formulation.

## What the example estimate actually is

For **RCS 8225W White Lightning and RCS 8125RC Redline**, the register includes explicitly asterisked, gas-phase **Sutherland/Eucken fallback** values using each RCS-reported chamber temperature, molecular weight and gamma. It reproduces the form already present in `srm_1d.propellant.create_gas_properties_estimated`:

- Specific gas constant: `R_s = R_universal / M`
- Ideal-gas approximation: `Cp* = gamma * R_s / (gamma - 1)`
- Viscosity*: temperature-scaled Sutherland approximation used by the repository
- Conductivity*: the repository's modified Eucken approximation
- Prandtl*: `mu* Cp* / k*`

The asterisks indicate **model estimates**, not direct RCS measurements, experimental fits, or NASA CEA/RPA results. A generic fallback is **not a proxy for validated equilibrium-reacting mixture transport**, especially where high-temperature species and condensed products matter. The model code itself warns of substantial conductivity and Prandtl uncertainty.

No confidence interval is assigned: the dominant uncertainty is model discrepancy / unknown real chemistry, not a small statistical error. Do not silently insert these values as verified transport into a `.ric` file.

## Display and ingestion rules

- Display `*` alongside **every** estimated field and any output derived from such fields. UI tooltip: `* Estimated; source chemistry and/or transport model is unverified. See provenance.`
- Keep source IDs and source-specific dates in machine-readable metadata. Do not apply a single confidence label to an entire propellant when some fields are reported and others estimated.
- Distinguish `manufacturer_reported`, `calculated_estimated`, `unverified_historical`, `unreported`. `manufacturer_reported` means provenance is verified, **not** that an uncertainty bound has been measured.
- Never infer identity from similar names alone (`Red` ≠ proven `Redline`, `White` ≠ proven `White Lightning`, `Blue` ≠ proven `Blue Thunder`).
- Do not coerce a single approximate conductivity/Cp pair into separately calculated frozen and equilibrium variants.
- Keep an explicit user confirmation boundary before promoting *-flagged data from a reference register into solver inputs.

## Next research/implementation boundary

RCS publishes individual propellant data sheets. Additional verified chamber-state fields can be imported from each **specific** sheet and reconciled with its product code. Where chemical composition remains proprietary, treat calculations as scenarios; do not label them as commercial product validation. NASA CEA 3 can provide transport only for the **specified hypothetical reactant set** and thermodynamic conditions. Even a perfectly converged CEA run cannot establish that a third-party forum recipe is the genuine product formulation.

**No `.ric` values or solver defaults are changed by this documentation-only contribution.**
