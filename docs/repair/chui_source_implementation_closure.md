# CHUÍ — patch 49 source-completion record (supersedes the patch-48 self-assessment)

Date: 2026-11-05. Scope: static source edits only.
**No compilation, tests, benchmarks, simulations, or scientific validation were performed in this campaign.**

**Status: NOT SOURCE-COMPLETE.** The patch-48 section below is retained as
history; its matrix was an unverified self-declaration and is superseded by the
matrix in this section (independent audit
`CHUI_patch48_auditoria_independente_fonte_2026-11-04.md`).

## Contradictions closed (P0 + additional findings)

| ID | Defect | Resolution (symbol) |
|---|---|---|
| P0-01 | stale `source_evaluation_scale_factor` / undefined `smooth_source_dt_code` | `SourceStagePlan` + `makeSourceStagePlan()`; every operator reads its clock from the plan |
| P0-02 | `executeThermochemistry` returned a report without value; caller arity | returns `void`; observations only accumulate (`accumulateSourceIntervalObservation`), publication in `closeSourceInterval()` |
| P0-03/04 | `const amr_options` mutated; undefined adapter/context types | mutable options; `RecursiveAmrLevelSourceAdapter` (hydro_amr_runtime.cpp); legacy callback API removed |
| P0-05 | `setMultilevelSourceRuntime ... override` with no base | removed; `bindSourceRuntime()` (already called by `reference_runtime_composition.cpp`) binds and calls `setAmrLevelSourceCouplingEnabled()` |
| P0-06 | CTU ghost telemetry fields undeclared | declared in `HydroProfileEvent` and `ProductionAmrHydroDiagnostics`; patch/level aggregation and global rank sums |
| P0-07 | private `invalidateOccupancyCache()` called from rollback | public `DistributedParticleIdRegistry::invalidateForRollback()` |
| new | `persist_v02_regrid_state` called, never defined; applied regrid never published V02 | defined; AMR half of V02 published on no-op and applied passes |
| new | `ProductionAmrRegridPlan::optional_admitted_*`, `optional_min_admitted_benefit` undeclared | declared |
| new | predictive AMR block: `catch` body swallowed the local forecast; `predictive_enabled` undeclared | restructured, admitted before allocation |
| new | `admitCompositionSparseActivation` called at three sites, never defined | defined (retained-capacity transaction, canonical class) |

## Production reachability (traced call flow)

* `.param.txt` → `SimulationConfig` → `reference_runtime_composition.cpp`
  builds hydro then source owners and calls `HydroAmrRuntime::bindSourceRuntime()`
  → coupling enabled iff `numerics.amr_time_stepping = recursive_subcycling`.
* `StepOrchestrator` dispatches `kSourceTerms` twice (`kStrangPreHalf`,
  `kStrangPostHalf`) in both the single-step and hierarchical-block paths →
  `SourceRuntimeImpl::execute()` → plan → smooth halves → (post only)
  `executeDiscreteSourceEvents()` → `closeSourceInterval()`.
* `runProductionAmrHydroPath()` → recursive driver
  (`advanceProductionAmrHydroRecursive[Impl]`, serial and distributed) → per level
  `RecursiveAmrLevelSourceOperator::apply(kPreHydroHalf)` before history start
  and `apply(kPostHydroHalf)` after sweep/register exchange →
  `executeAmrLevelSmoothSources()` → `executeThermochemistryInterval()`.
* Regrid cooling-length criterion: accepted chemistry history, else
  `SourceRuntime::fillAmrCoolingTimes()` (same rate set), else legacy provider.

## H01–V02 source-coverage matrix (patch 49)

| ID | Status | Evidence / exact remaining gap |
|---|---|---|
| H01 | Implemented in source, not tested | CTU ghost-slab increments + counted fallback; telemetry now declared/aggregated. Formal order/convergence: later validation |
| H02 | Implemented in source, not tested | dual-energy, floors, comoving factors (unchanged) |
| C01 | External data dependency | hybrid network + table loaders; production UVB/metal tables not supplied |
| C02 | Implemented in source, not tested + **External data dependency** | `metal_ion_network.*`: C/O/Si/Fe stage inventory, validated `chui_metal_ion_rates_v1` loader, ladder solver, layout/transport/diffusion/snapshot/restart-signature/config. No licensed, source-verified rate dataset exists in the repository; charge exchange with H/He not represented (schema v1) |
| C03 | Implemented in source, not tested | residual/error history, hysteresis/dwell (physical years), sparse mapping; ion ladders initialized from equilibrium on promotion (counted) |
| C04 | Implemented in source, not tested | shadow non-committing; controller published once per Strang interval |
| I01 | Implemented in source, not tested | regime split unchanged |
| I02 | Implemented in source, not tested | SGS shear transfer now receives the defined half duration |
| S01 | Implemented in source, not tested | births only in post half at `t_end`; global synchronized (not per level) |
| S02 | External data dependency | cursor/chronology once per interval; production yield tables absent |
| F01 | Implemented in source, not tested | feedback once per interval (post half) |
| F02 | External data dependency | winds/pre-SN once per interval; production wind dataset absent |
| D01 | Implemented in source, not tested | capability authorizes `core_elements`, refuses `total_only` on >1 rank; ion lanes diffuse with elements |
| A01 | Implemented in source, not tested | mandatory/optional budgets, transactional publication (patch 48, preserved) |
| A02 | Implemented in source, not tested | recursive distributed reflux (patch 48, preserved); level Strang ordering w.r.t. histories now explicit |
| A03 | Implemented in source, not tested | cooling-length fallback through `fillAmrCoolingTimes` |
| A04 | Implemented in source, not tested (research) | predictive block repaired |
| G01 | Conditional on later qualification | TreePM/KDK retained; no multilevel correction without evidence |
| T01 | Implemented in source, not tested | `SourceStagePlan`, once-only events, single publication, rollback of open interval |
| T02 | **Partial, exact missing code** | smooth thermochemistry/ISM/dust/SGS per level clock (pre+post). Diffusion, SGS shear transfer, SF, SSP, winds, SN, BH remain globally synchronized at coarse halves/endpoint. Missing: level-local cadence for element diffusion/SGS transfer with coarse/fine flux registers for composition-diffusion fluxes, and per-level event ownership for SF/feedback |
| M01 | **Partial — not exhaustively traced** | sparse species/ion storage, transport sparse activation (`admitCompositionSparseActivation`), element-diffusion workspace (`elementDiffusionWorkspaceBytes` over all lanes, including ion lanes) and predictive-AMR points are admitted before allocation. Not traced end-to-end in this campaign: that a single admission decision accounts for the coexisting rollback copy + PM + HDF5 + MPI buffer peak, and that `parallel/node_memory_pool` node-aggregate admission covers every phase-resident owner |
| M02 | Implemented in source, not tested | CPU serial/OpenMP/MPI unchanged; GPU optional |
| V01 | Implemented in source, not tested | ledgers/restart; four modes untouched |
| V02 | Implemented in source, not tested (research) | chemistry + AMR halves of the joint controller both publish; local estimators, not global error |
| X01 | Out of scope | no MHD/RT/new AGN |

## Remaining validation gates (owner-reserved)

Build all presets; unit tests (`test_amr_recursive_subcycling` updated for the
operator seam); Strang temporal convergence with mixed sources; exactly-once
event counts per interval under retry; recursive vs synchronized AMR source
equivalence; metal-ion ladder conservation/equilibrium against a qualified
dataset; restart round-trip with ion lanes; MPI rank-count independence.

---

# Patch-48 history (unverified self-assessment, superseded)

# CHUÍ — fechamento estático de implementação H01–V02

Data: 2026-11-04. Escopo: inspeção e edição exclusivamente estáticas do patch
47. Conforme determinação do proprietário, **não houve compilação, teste,
benchmark, simulação, sanitizer ou validação numérica** nesta campanha.

## Resultado executivo

Esta campanha fechou em fonte as contradições concretas do patch 47 e os
blocos subsequentes de integração:

1. **A02:** o refluxo recursivo distribuído deixou de ser proibido pelo sweep.
   As contribuições remotas entram no `PendingFluxRegisterStore` do proprietário
   grosso com intervalo/subpasso explícito, validação de duplicidade por origem,
   área espacial somada entre ranks no mesmo subpasso e aplicação somente após
   cobertura temporal completa. O store já é migrado e serializado no restart.
2. **T01:** a ordem canônica passou a
   `S(dt/2) -> H(dt) -> S(dt/2)`. Termoquímica, dust, SGS e difusão usam os
   intervalos físicos de meia etapa; nascimento, evolução/eventos, feedback,
   ventos e BH executam uma única vez no post-step no epoch absoluto final.
3. **D01:** a capability agora reflete o operador MPI elementar realmente
   existente e aceita apenas layout `core_elements` em build MPI. O passo de
   difusão recebe `dt/2` em cada metade Strang. O layout escalar `total_only`
   continua explicitamente recusado em múltiplos ranks.
4. **A01:** refinamentos mandatórios não consomem o orçamento opcional. O
   runtime prepara a hierarquia inteira numa cópia admitida pelo governador,
   vota coletivamente e publica por move no-throw somente após sucesso global.
5. **H01/I02:** CTU transversal em halos AMR/MPI e transferência SGS por
   cisalhamento foram ligados aos caminhos conservativos, scratch governado,
   reflux/regrid e telemetria existentes.
6. **M01:** as lanes químicas detalhadas adaptativas passaram a armazenamento
   esparso por `gas_cell_id`, com referência fixa densa, capacidade governada,
   transporte, migração, regrid, rollback e restart schema 3.
7. **A03/A04/V02:** a política força/resolução agora toma ação explícita; o
   forecast de eventos cruza ranks; e um controlador persistente selecionável
   coordena shadow químico, benefício AMR e custos tipados de cell-step,
   comunicação, churn e RAM sem substituir as referências fixas.

A inspeção pós-edição não encontrou chamadas desconectadas nas novas
assinaturas, gates de capability contraditórios, whitespace inválido ou a
antiga guarda A02. Ela encontrou e corrigiu ainda o teste de ordem que assumia
Lie splitting e a interseção geométrica que confundia contato de aresta/canto
com área positiva de face.

## Matriz normativa (24 requisitos)

`IMPLEMENTADO EM FONTE, NÃO TESTADO` significa apenas cobertura estática do
caminho normal de produção; não é qualificação científica. `PARCIAL` significa
que permanece código obrigatório ausente.

| ID | Status nesta árvore | Evidência e limite preciso |
|---|---|---|
| H01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | MUSCL–Hancock/PPM/CTU/HLLC-HLLE permanecem em `src/hydro/hydro_core_solver.cpp`. O CTU agora reconstrói, com o mesmo PCM/PLM/PPM e Riemann do avanço, os dois fluxos de cada eixo tangencial do ghost a partir do halo imutável já preenchido por AMR/MPI; o incremento usa o volume-alvo e não cria autoridade persistente. Apenas a região de aresta/canto cujo stencil cruza uma segunda interface usa fallback normal conservativo explícito, contado em `ctu_transverse_degraded_states`; correções/axes remotos têm telemetria própria. Scratch do halo foi incluído em `HydroScratchBuffers::workspaceBytes`, na admissão AMR e nos diagnósticos runtime. Prova de ordem e estabilidade remota pertencem à campanha numérica posterior. |
| H02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Dual energy, autoridade de energia interna, recuperação de primitivas, floors e conversões comoventes: `src/hydro/hydro_core_solver.cpp`, `src/amr/amr_ghost_fill.cpp`, `src/hydro/gas_thermodynamic_refresh.cpp`. Falta qualificação numérica, não uma rota conhecida de fonte. |
| C01 | **BLOQUEADO POR DADOS EXTERNOS** | Operador, loader SHA-256/provenance, H/He não-equilíbrio, heating/shielding e cooling elementar estão em `src/physics/thermochemistry*.cpp` e `metal_cooling_table.cpp`. As tabelas físicas de produção exigidas pelos perfis não acompanham o pacote. |
| C02 | **PARCIAL / BLOQUEADO POR DADOS EXTERNOS** | Backends analítico/equilíbrio existem; `fixed_detailed` em `src/physics/chemical_network.cpp` permanece rede reduzida primordial + C+/CO. Faltam inventário multiiônico aprovado e taxas/tabelas verificáveis; sem esses dados não é lícito inventar reações, mas a ligação do inventário definitivo também permanece código futuro. |
| C03 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Temperatura, densidade, H2/CO, cooling time, erro/resíduos, histerese, dwell e transição conservativa estão em `ThermochemistryOperator::decideChemistryFidelity`/`BaryonicSourceStage`. `GasCompositionSidecar` schema 3 preserva o histórico e a extensão esparsa por identidade através de commit, transporte MPI, migração, regrid, rollback/restart e HDF5. `SourceRuntimeImpl::updateGlobalChemistryControl` reduz o shadow global e, em modo V02 conjunto, grava bias, amostras e promoções multinível no controlador canônico. A qualificação dos limiares continua pendente de validação, não de código conhecido. |
| C04 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Shadow detalhado amostrado, determinístico, limitado e sem commit: `ThermochemistryOperator::advanceRow`; redução global em `SourceRuntimeImpl::updateGlobalChemistryControl`. O controlador agora avança uma vez por intervalo Strang, não entre as duas metades. |
| I01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Regime efetivo versus frio resolvido, ownership térmico único e ausência de double cooling: `effective_multiphase_ism.cpp`, `baryonic_source_stage.cpp`, `gas_thermodynamic_refresh.cpp`. |
| I02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `prepareSgsShearTransfer` em `ism_dust_sgs.cpp` usa `sgs_eddy_viscosity_coefficient` no stress traceless `nu_t=C_v Delta^2 |S*|`: impulsos e trabalho são antissimétricos por face, a produção sai da energia resolvida e entra na lane SGS no mesmo commit, e o limite viscoso pede retry numérico. `SourceRuntimeImpl::executeDiffusionAndSgsShear` reutiliza o grafo físico local/coarse-fine/periódico/MPI, valida interfaces remotas coletivamente, governa staging e refresca caches. Advecção/reflux/regrid/restart, feedback e dissipação permanecem na autoridade SGS existente. O modelo continua EXPERIMENTAL/UNQUALIFIED e requer validação numérica futura, não código obrigatório conhecido. |
| S01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Elegibilidade Jeans/fluxo/boundedness, vizinhos AMR/MPI, nascimento conservativo, ID determinístico e birth epoch final: `star_formation.cpp`, `star_formation_gradients.hpp`, `SourceRuntimeImpl::execute`. |
| S02 | **BLOQUEADO POR DADOS EXTERNOS** | SSP, idade/metallicidade, cronologia cumulativa, associações e loader de yields existem em `stellar_population.cpp` e `source_runtime.cpp`. O schema 2 da sidecar populacional preserva também o cursor absoluto do último intervalo aceito, recusa gaps/replays e migra schema 1; ordinais e cursor só avançam no mesmo commit dos depósitos. Faltam tabelas licenciadas de produção; nenhuma equivalência FIRE-3/COLIBRE é declarada. O depósito em barreiras AMR multinível continua limitado por T02. |
| F01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Ramos mecânico resolved/unresolved, orçamento exclusivo térmico/SGS, depósito conservativo e refresh termodinâmico: `supernova_feedback.cpp` e `SourceRuntimeImpl::executePart2StellarPopulation`. |
| F02 | **PARCIAL / BLOQUEADO POR DADOS EXTERNOS** | Vento e pre-SN, loader e orçamento não aditivo existem em `stellar_wind.cpp`/`source_runtime.cpp`; falta dataset físico de vento de produção e a cronologia por nível exigida por T02. |
| D01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Transporte conservativo de elementos/dust/moléculas/SGS por `CompositionFluxAccumulator`; difusão MPI por `DistributedElementDiffusionCoupling`; capability coerente em `buildRuntimeCapabilityReport`; integração Strang usa `source_dt_code`. |
| A01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `planProductionAmrRegrid` separa `admitted_optional_operations`; `runAutomaticAmrRegrid` usa candidato integral governado + voto coletivo + publicação no-throw. Prolong/restrict de composição e identidades continuam nos operadores existentes. |
| A02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `advanceProductionAmrHydroRecursiveImpl`, `mergeFluxRegistersIntoPendingStore`, `DistributedRecursiveHooks::exchangeLevelFluxRegisters`, histories start/end e restart de `PendingFluxRegisterRecord`. Intervalos incompatíveis e duplicatas por origem falham fechados; área fina repartida é recomposta no subpasso. |
| A03 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `patchDecision` mantém massa/Jeans/cooling length como mandatórios e jump, metal/species/vorticidade/eventos como opcionais tipados. Todo backend Part 2 publica `chemistry_last_cooling_time_seconds`; `runAutomaticAmrRegrid` o converte ao tempo-code autoritativo. `AmrForceResolutionPolicy` oferece a referência histórica contabilizada e a política produtiva que bloqueia refinamento opcional abaixo do softening e falha coletivamente se requisito mandatório violar softening, `amr_max_level` ou a hierarquia 2:1. G01 continua condicional à validação. |
| A04 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `planProductionAmrRegrid` calcula severidade normalizada por critério e custos separados de cell-steps incrementais, faces MPI, churn e RAM; o modo conjunto ordena benefício/cell-step e usa MPI/churn como desempate, aplica piso com histerese e orçamento, enquanto mandatórios bypassam os tetos. `runAutomaticAmrRegrid` troca posições atual/driftada por `allgatherBytesBounded`, resolve o host somente no owner e propaga falhas/admissão coletivamente. O modelo permanece RESEARCH/UNQUALIFIED até validação numérica. |
| G01 | **CONDICIONAL À VALIDAÇÃO** | TreePM/KDK e auditoria AMR-softening existem. A especificação condiciona uma correção gravitacional multinível adicional ao resultado de qualificação futura; esta rodada não executou essa qualificação e não inventou novo solver. |
| T01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `SourceSplitPhase`, ordem de nove dispatches em `time_integration.cpp`, inversão determinística de tempo cosmológico, `makeSourceInterval` por metade e eventos somente post. Rollback global restaura estado, schedulers, force endpoint e checkpoint privado da fonte. |
| T02 | **PARCIAL** | `advanceProductionAmrHydroRecursiveImpl` agora invoca `SourceRuntimeImpl::executeAmrLevelSmoothSources` em ambas as metades de cada subintervalo, inclusive nos ranks sem linhas locais; termoquímica/dust/SGS usam o relógio `dt_L`, entram nos históricos temporais e shadow, promoções adaptativas e controlador são agregados uma vez por passo com rollback privado. A autoridade de eventos permanece O(1) na sidecar populacional e ganhou cursor absoluto persistente com rejeição de gap/replay e migração de schema. SF/SSP/feedback/ventos/difusão ainda usam a barreira global enquanto faltam depósitos cross-level e difusão no grafo multinível. |
| M01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `GasCompositionSidecar` mantém apenas as lanes básicas densas e armazena C+/CO adaptativos em vetores contíguos ordenados por `gas_cell_id`, com projeção implícita conservativa e retenção do histórico após rebaixamento. Crescimento em química, transporte Godunov MPI, difusão e demais commits usa `RetainedCapacityTransaction`/governador; regrid e migração dimensionam exatamente o candidato já admitido. Pack schema 3, restart com migração schema 1/2, HDF5, AMR prolong/restrict, ownership, rollback e depósitos de feedback/poeira preservam a mesma autoridade. O modo detalhado fixo continua denso como referência determinística. Qualificação de pico MPI fica para a campanha de testes reservada. |
| M02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | CPU scalar/SIMD/OpenMP/MPI, decomposição/migração e balanceamento ponderado permanecem nos runtimes existentes. GPU segue opcional e não foi transformada em dependência. |
| V01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Contratos de configuração/provenance, ledgers, snapshots HDF5, restart versionado e quatro decks científicos permanecem conectados em `config.cpp`, `snapshot_hdf5.cpp`, `restart_checkpoint.cpp` e `configs/chui_part3`. Aceitação end-to-end foi deliberadamente adiada. |
| V02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `V02ControlMode` seleciona `fixed_reference` ou `joint_error_resource_v1`. O estado versionado em `adaptive_resource_control.hpp` reside na sidecar de `SimulationState`, portanto acompanha rollback/restart sem segunda autoridade. `SourceRuntimeImpl` publica shadow, bias e demanda química multinível; `runAutomaticAmrRegrid` reserva crescimento esparso antes de expor headroom opcional, persiste decisões admitidas ou no-op e usa estimadores AMR tipados sem inventar erro físico universal. O modo fixo preserva bias químico e ordem AMR por `patch_id`; física mandatória nunca é reduzida pelo controlador. O modo conjunto permanece RESEARCH/UNQUALIFIED até qualificação. |

## Contratos alterados nesta campanha

- `StepContext::source_split_phase` é transitório e não cria nova autoridade de
  tempo. `IntegratorState` e `CosmologicalTimeline` continuam autoritativos.
- `PendingFluxRegisterRecord` continua a única autoridade persistente do reflux
  incompleto; não foi criado store paralelo. A semântica de
  `fine_area_accumulated_comov` passou a ser a cobertura espacial do subpasso
  corrente/completo, e não o máximo de uma contribuição por rank.
- O regrid automático conserva uma única publicação autoritativa de
  `SimulationState`; o candidato é scratch transacional governado e destruído
  após commit/rejeição.
- Difusão distribuída somente é elegível com composição elementar. A capability
  não promove o operador escalar legado nem oculta a dependência de MPI.
- `chui_v02_joint_error_resource_control_v1` é o único estado persistente do
  controlador conjunto. Química continua dona de composição/fidelidade e AMR
  continua dono da hierarquia; o sidecar guarda apenas histórico de decisão e
  reservas, participa da cópia de rollback/restart e não duplica campos físicos.

## Pendências obrigatórias verdadeiras

Não é correto declarar fechamento integral H01–V02 nesta árvore. Permanecem
obrigatórias a implementação T02 multinível bariônica e a ligação do
inventário definitivo de C02. C02 requer primeiro a especificação/dados verificáveis do
inventário detalhado; S02/F02 requerem tabelas físicas de produção. G01 é
deliberadamente condicional à campanha de validação prevista pela norma.

X01 permanece fora de escopo: nenhuma implementação de MHD, moving mesh, SPH,
CHIMES, Grackle, RT multifrequência ou nova família AGN foi introduzida.
