# Radiation, black holes and AGN (CHUI-RT-BH-AGN-REF-v1, pass 01)

**Status:** source implementation, UNVALIDATED. No build, test, benchmark or
simulation has been run for this module (user-reserved source-only pass).
Nothing here is verified, calibrated or production-qualified. Detailed
coverage: `docs/architecture/chui_rt_bh_agn_implementation_pass01_status.md`;
decisions: `docs/architecture/chui_rt_bh_agn_compatibility_amendment_v01.md`.

## Model selection (`.param.txt`, `[physics]`)

| Key | Values (default first) | Meaning |
|---|---|---|
| `bh_engine_model` | `legacy_host_cell`, `reference_v1` | legacy host-cell Bondi/Eddington/thermal callback (unchanged) or the reference engine |
| `bh_supply_model` | `bondi_torque`, `bondi`, `bondi_angular_suppressed`, `torque`, `resolved_flux` (P1, refused) | requested supply at the disk boundary |
| `bh_disk_model` | `scalar_reservoir`, `instantaneous` | scalar viscous reservoir or the legacy-equivalent instantaneous comparator |
| `bh_mass_ledger` | `legacy_newtonian`, `relativistic_defect` | `dM_bh = dm0` or `dM_bh = (1-eps_acc) dm0 - E_spin/c^2` |
| `bh_eddington_cap` | `reference_rate`, `none` | inner inflow cap at `mdot_Edd,ref = L_Edd/(eps_ref c^2)`; must agree with `bh_use_eddington_cap` |
| `bh_dynamics_model` | `none`, `unresolved_chandrasekhar` | Tremmel-range Chandrasekhar correction |
| `bh_merger_model` | `none`, `bound_pair` | bound, converging/unresolved pairs only |
| `bh_seed_model` | `none`, `halo_fof_serial` | single-rank FoF halo seeding (distributed: refused) |
| `bh_spin_model`, `bh_accretion_state_model` | `none` (P1 values refused) | Kerr/BZ and accretion-state sources exist but are not admitted |
| `agn_feedback_model` | `legacy_thermal_host`, `channels_v1` | channels consume only the engine's pending budgets |
| `agn_*_efficiency` | `>= 0`, sum + `bh_epsilon_r` < 1 | fractions of `dm0 c^2` per channel |
| `agn_radiation_coupling` | `none`, `rt_transport`, `local_empirical` | radiative budget escapes, is transported once, or a declared fraction heats locally (mutually exclusive) |
| `rt_enabled` | `false`, `true` | on-the-fly multigroup M1 |
| `rt_spectral_profile` | `ionizing3`, `grey`, `broad_uv_ir` (refused), `agn_hard` (refused) | group edges (physical eV) |
| `rt_opacity_model` | `photoionization` (refused, see amendment A5), `constant_reference` | |
| `rt_integrator` | `imex_pd_ars`, `explicit_reference` | PD-ARS (He et al. 2024) or split comparator |
| `rt_c_fraction` | `(0, 1]` (1) | globally constant reduced speed; `< 1` refused in cosmological modes |

Units: kernel radii are PHYSICAL kpc (`*_kpc_phys`), masses code units
(`*_code`), times years, velocities km/s, opacities cm^2/g.

## Runtime state map

| Truth | Owner | Lifetime / persistence |
|---|---|---|
| Horizon mass, disk mass/composition/J, pending budgets, cumulative ledgers | `chui_bh_engine_v1` module sidecar (per BH particle ID) | restart (generic module block), migrates with the particle |
| Radiation moments `E, F` per group | `SimulationState::radiation` (`RadiationMomentSoa`) | restart v27 `/state/radiation` |
| Gas mass/momentum/energy/composition | hydro lanes (unchanged authority); mutated only through `SourceTransactionBatch` | — |
| Magnetic face flux | MHD (never touched by BH/AGN/RT) | — |
| Transported photo-rates | derived cache (`RadiationPhotoRateCache`) | not persistent |

## Chronology (one source interval)

At the SourceRuntime discrete-event epoch (coarse interval end, hierarchy
synchronized): star formation → stellar returns/feedback → BH engine
(4 collective support rounds: radius bins, kernel statistics, counter-rotation,
withdrawal prepare/vote/commit) → disk/horizon/budgets → DF → mergers →
seeding → AGN injection (statistics round, commit round) → RT interval
(AGN radiative budget emitted exactly once; PD-ARS subcycles; staged commit).
The legacy BH callback does not run when `reference_v1` is selected.

## Known limitations (pass 01)

See the gap ledger in the status document. In particular: photoionization RT,
stellar RT sources, multi-rank RT, RT with regrid, distributed seeding, P1
spin/state/jet-power models, nozzle and magnetic injection are refused at
admission; mergers and DF are rank-local; no science-snapshot export of the
new BH lanes; no measured accuracy of anything.
