# CHUI discrete-event chronology: physical event time vs publication epoch

Status: source contract (patches 51–53). Not numerically qualified; nothing
in this document was executed.

## Definitions

| Term | Meaning | Source authority |
|---|---|---|
| physical event time | Model-defined instant of a discrete stochastic SN: the stellar age `a_k = inf{a : Λ_c(a) ≥ Λ_k}` at which the star's birth-mass-scaled cumulative expected count of channel `c` (CC or Ia) first reaches the pre-drawn threshold `Λ_k` of ordinal `k` | `physics::stellarEventAgeYears`, `physics::realizedEventChronology` (`src/physics/stellar_population.cpp`) |
| publication epoch | The synchronized end of the accepted coarse (or hierarchical-block) source interval: every AMR level and particle bin is synchronized, no flux register, temporal history or coarse/fine split is open | `SourceRuntimeImpl::executeDiscreteSourceEvents` epoch guard |
| publication lag | `age_end − a_k ≥ 0` (years) | `DiscreteEventChronologyLedger::accumulate` |
| forecast | Earliest physical time of the next not-yet-realized event over the local stars, from the persistent population records and star ages | `SourceRuntime::discreteEventStepEndLimitTimeCode` |

Deposition is **not** performed at the physical event time. Energy, momentum
and ejecta of every event are deposited at the publication epoch; the
physical time is recorded (earliest/latest event time, lag statistics) so the
chronology is explicit and auditable, not hidden.

## Event families

| Family | Physical-time meaning | Publication | Lag policy |
|---|---|---|---|
| Star birth (S01) | interval end by model (birth probability integrated over the interval) | interval end, age 0 (F5) | none (identical instants) |
| Continuous SSP returns, table winds/pre-SN | integrated over the absolute age interval `[age_b, age_e]`; no intra-step instant exists in the model | interval end | first-order operator splitting |
| Stochastic CC/Ia SN (`ssp_stochastic_sn_events = true`) | exact threshold crossing `a_k` | interval end | recorded; bounded when `physics.ssp_max_event_publication_lag_years > 0` |
| Expectation-mode SN | fractional expected counts, interval-integrated | interval end | none (no discrete time) |
| BH seeding / accretion / feedback | seeding on the synchronized end state; rates × interval | interval end | none |

## `physics.ssp_max_event_publication_lag_years`

* `0` (default) means **the extra bound is disabled**: the lag is bounded only
  by the accepted coarse interval and is reported every interval
  (`event_lag_policy = bounded_by_accepted_interval_only`). It does **not**
  mean "zero lag required". No universal numeric limit is imposed; choosing a
  positive value is a scientific decision of the run owner.
* `L > 0`: before each global step every rank evaluates its forecast; the
  time coordinator reduces it with MIN over ranks (guarded phase + collective
  vote, then `allreduceMinDouble`) and clips the step end to
  `t_next_event + L` (`time.source_event_publication_clip`, ordered with
  snapshot clips; the snapshot clip keeps its restart-resume dt bookkeeping).
  Every event realized in that step therefore has lag ≤ L, and the stage
  refuses (fail closed) an interval whose forecast-covered event exceeds
  `L·(1 + 10⁻⁶)`. The criterion never shortens a step below `L`.
* Requires a Part 2 `stellar_population_model` with stochastic SN events.

## Exceptions, rollback and restart

* **Legacy catch-up.** A population record without an absolute cursor (stars
  from older checkpoints or initial conditions) adopts the current interval as
  its first cursor. Events whose threshold was already crossed at `age_b` are
  counted as `catch_up_discrete_events` and are not timed. Timed events of a
  record whose cursor was adopted in this interval are counted in
  `unforecast_timed_discrete_events` and are exempt from the bound (the
  step-start forecast could not see them).
* **Rollback.** The forecast is a pure function of the persistent records and
  star ages in `SimulationState`; the per-interval ledger is published only
  with the stage commit and is restored by the step rollback checkpoint. A
  rejected or retried step re-derives the identical chronology.
* **Restart.** No extra restart state exists: thresholds, ordinals and cursors
  are in the population sidecar. A checkpoint written before the key existed
  resumes only with the default `0`.
* **Refinement / ranks.** `a_k` depends only on (birth key, channel, ordinal,
  birth mass, birth metallicity, SSP model): it is invariant under AMR
  refinement, rank count and migration.

## Validation still required (not source gaps)

Event-age subdivision invariance (`testEventPhysicalTimeIsSubdivisionInvariant`,
written, not run); multi-rank forecast/clip agreement; restart equivalence with
`L > 0`; scientific assessment of the deposition-lag impact for the selected
feedback model and resolution.
