# AMR Temporal Ghost Interpolation

Date: 2026-06-18

## Implemented contract

CHUÍ records coarse patch conserved state at both ends of an open local coarse interval and uses that
history for coarse-to-fine ghosts during local AMR hydro subcycling. The ghost-fill request carries a target
state time, requested ghost-fill time, current source-state time, temporal-history pointer, and explicit
coarse-to-fine temporal enablement. The requested time is the beginning of the current fine hydro update.
That convention is deliberate: the production `HydroCoreSolver` consumes ghost state before the
MUSCL-Hancock reconstruction/predictor path, so this layer must not guess a separate half-step external
boundary time.

Interpolation is performed on conserved state only:

```text
U(t_fill) = (1 - alpha) U_start + alpha U_end
alpha = (t_fill - t_start) / (t_end - t_start)
```

Endpoints are exact. A tiny endpoint tolerance handles floating-point representation only; genuine
extrapolation is rejected. The interpolated conserved state is checked for finite positive density and then
converted through the existing ideal-gas primitive recovery. A non-finite or non-positive recovered state is
a hard error.

## Stable identity and lifecycle

A temporal history record contains patch ID, patch level, geometry fingerprint, gas identity generation,
interval, completion marker, and one record per source patch cell. Per-cell records contain stable
`gas_cell_id`, patch-local cell index, and conserved start/end state. The fill path resolves both stable
identity and patch-local location. It does not use a dense storage row as physical identity.

A history becomes invalid when geometry or identity generation changes. The local implementation blocks
refine/derefine while histories are active. It has no MPI migration/remapping path. Same-level ghosts
require equal source/target time. Fine-to-coarse ghost reads require synchronization time and otherwise
reject with a named temporal-misalignment status.

## 2026-11-04 source repair: publication and retirement

Start capture, including per-level capture with parent intervals already open,
builds one admitted candidate containing all new and retained records. End
capture also completes a candidate before publication. Allocation, admission,
geometry, identity or physical-state validation failure leaves all canonical
records and their storage unchanged. Candidates count outer and nested cell
capacities while the old store remains live; geometry scratch is separately
admitted for start capture. Mathematical interpolation and endpoint values are
unchanged. This is SOURCE ONLY; fault/conservation tests are NOT RUN — USER-RESERVED.

Both retirement functions now take an optional trailing `MemoryGovernor*`.
Pass the capture governor to reconcile freed nested capacity. Record compaction
uses nonthrowing moves/destruction and retains the charged outer capacity. Calls
without a governor remain supported, with accounting reconciliation owned by
their caller. `applyCompletePendingFluxRegistersToSimulationState` similarly
accepts a trailing defaulted governor for its candidate workspace. Existing
source calls remain valid; dependents must be rebuilt. No restart schema changes.

These are individual-operation guarantees. A failure after earlier hydro
substeps committed still requires the coordinator's whole-step rollback; general
distributed recursive integration remains open. See
[`../repair/chui_conformance_closure_20261104.md`](../repair/chui_conformance_closure_20261104.md).

## Restart

Schema `cosmosim_restart_v19` serializes active histories in
`/state/amr_temporal_boundary_history`. A restart written while a local interval is active must contain a
completed coarse start/end history before remaining fine substeps can resume. The HDF5 test
`integration_restart_equivalence_amr_temporal_ghosts` exercises this exact continuation path.

## Non-claims

The implementation is not distributed temporal ghost exchange, not arbitrary-depth Berger-Colella time
subcycling, and not a temporal interpolation model for remote fine patches. It does not permit history to
survive AMR regridding or patch migration by remapping; it fails safely instead.
