# Black-Hole and AGN baseline module

This document defines the baseline `black_hole_agn` model used by CosmoSim for the first full-physics BH/AGN path.

## Scope and model choices

- This is a **single explicit baseline model** (no alternative accretion or feedback closures mixed into this path).
- The module is controlled by `physics.enable_black_hole_agn` and defaults to **off** for lightweight modes.
- All calibration-sensitive parameters are explicit in config normalization and module sidecar metadata.
- The baseline uses a Newtonian conservative mass ledger: `mdot_acc_code` is gas rest mass transferred into BH dynamical/subgrid mass. `epsilon_r` controls the unresolved radiative-energy yield and the Eddington cap; it does not silently delete tracked gravitating mass.

## Equations and frame conventions

The module uses:

- Bondi-like rate:
  \[
  \dot M_\mathrm{Bondi} = \alpha\,4\pi G^2 M_\mathrm{BH}^2\rho_\mathrm{phys}\left(c_s^2 + v_\mathrm{rel}^2\right)^{-3/2}
  \]
- Eddington rate:
  \[
  \dot M_\mathrm{Edd} = \frac{4\pi G M_\mathrm{BH} m_p}{\epsilon_r\sigma_T c}
  \]
- requested accretion cap:
  \[
  \dot M_\mathrm{req} = \min(\dot M_\mathrm{Bondi},\dot M_\mathrm{Edd})
  \]
- actual transfer is additionally capped by available host gas mass, leaving the positive mass floor;
- feedback power from the **actual** transferred mass rate:
  \[
  \dot E_\mathrm{fb} = \epsilon_f\epsilon_r\dot M_\mathrm{acc}c^2.
  \]

Cosmological hydro density is stored in the repository's comoving density convention. At a source-stage evaluation scale factor `a`, the Bondi closure therefore receives

\[
\rho_\mathrm{phys}=\rho_\mathrm{stored}/a^3.
\]

Gas and BH velocity lanes are physical peculiar velocities, so `v_rel` is computed directly from their vector difference without an additional power of `a`.

## Conservative transaction

For each accepted accretion transfer `Delta M`:

- host gas mass decreases by `Delta M`;
- BH subgrid and dynamical particle masses increase by exactly `Delta M`;
- host density is reduced by the same mass fraction at fixed cell volume;
- gas metal mass is removed proportionally with swallowed gas and accounted in module counters;
- BH velocity is mass-weighted with the swallowed gas velocity so the local gas+BH linear-momentum ledger is conserved by the transfer;
- generated feedback energy is tracked separately from deposited feedback energy;
- `gas_cells.internal_energy_code` is **specific** internal energy, so deposited total feedback energy is divided by the remaining host gas mass before updating that lane.

BH seeding is also a gas-to-BH transaction. An eligible host must contain more than the seed mass plus the mass floor. Seed mass is removed from the host gas, density/metal mass are reduced consistently, and the newborn BH inherits the gas peculiar velocity. Production SourceRuntime obtains seed particle IDs through the same exact distributed particle-ID precommit registry used for other births and schedules the new particle for the next legal tick.

## Restart/schema/provenance implications

- Restart `state/black_holes` persists:
  - `host_cell_index`
  - `eddington_ratio`
  - `cumulative_accreted_mass_code`
  - `cumulative_feedback_energy_code`
  - `duty_cycle_active_time_code`
  - `duty_cycle_total_time_code`
- `physics.*` BH config keys are included in normalized config text and hash.
- Per-step module metadata sidecar key: `black_hole_agn`, including accreted/removed mass, swallowed metal mass, seeded mass, generated/deposited feedback and duty-cycle counters.

## Validation and remaining model limits

The implementation now closes the local mass/frame/energy-unit defects identified by the 2026-08 scientific workflow audit. This is **implementation closure, not empirical calibration**. Quantitative BH/AGN science still requires Bondi/Eddington unit tests in the configured code-unit system, convergence/calibration studies, and multi-rank ownership tests for host-cell migration and seeding.

The current baseline is host-cell local; it does not yet implement kernel-weighted accretion from a resolved gas neighbourhood, BH repositioning, BH mergers, or a separately tracked swallowed-metal reservoir. Those are explicit future model extensions rather than hidden behaviour.

## Reference engine (CHUI-RT-BH-AGN-REF-v1)

This legacy module is preserved unchanged as `physics.bh_engine_model =
legacy_host_cell` (default). The reference engine, disk reservoir, ledgers,
mergers, seeding and AGN channels are documented in `docs/rt_bh_agn.md`; the
two never run for the same source interval.
