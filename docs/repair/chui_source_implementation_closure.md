# CHUÍ — fechamento estático de implementação H01–V02

Data: 2026-11-04. Escopo: inspeção e edição exclusivamente estáticas do patch
47. Conforme determinação do proprietário, **não houve compilação, teste,
benchmark, simulação, sanitizer ou validação numérica** nesta campanha.

## Resultado executivo

Esta campanha fechou em fonte quatro contradições concretas do patch 47:

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
| H01 | **PARCIAL** | MUSCL–Hancock/PPM/CTU/HLLC-HLLE estão em `src/hydro/hydro_core_solver.cpp`; `transverseSourceCell()` ainda não transporta incrementos CTU de ghosts importados AMR/MPI e contabiliza degradação. A correção de área positiva foi feita em `src/amr/amr_hydro_orchestrator.cpp::intervalsOverlap`, mas a ordem multidimensional remota não está fechada. |
| H02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Dual energy, autoridade de energia interna, recuperação de primitivas, floors e conversões comoventes: `src/hydro/hydro_core_solver.cpp`, `src/amr/amr_ghost_fill.cpp`, `src/hydro/gas_thermodynamic_refresh.cpp`. Falta qualificação numérica, não uma rota conhecida de fonte. |
| C01 | **BLOQUEADO POR DADOS EXTERNOS** | Operador, loader SHA-256/provenance, H/He não-equilíbrio, heating/shielding e cooling elementar estão em `src/physics/thermochemistry*.cpp` e `metal_cooling_table.cpp`. As tabelas físicas de produção exigidas pelos perfis não acompanham o pacote. |
| C02 | **PARCIAL / BLOQUEADO POR DADOS EXTERNOS** | Backends analítico/equilíbrio existem; `fixed_detailed` em `src/physics/chemical_network.cpp` permanece rede reduzida primordial + C+/CO. Faltam inventário multiiônico aprovado e taxas/tabelas verificáveis; sem esses dados não é lícito inventar reações, mas a ligação do inventário definitivo também permanece código futuro. |
| C03 | **PARCIAL** | Temperatura, densidade, H2/CO, cooling time, histerese, dwell, transição conservativa e bias shadow persistente existem em `thermochemistry_operator.cpp`/`source_runtime.cpp`. A `GasCompositionSidecar` schema 2 agora retém temperatura/cooling time, erro normalizado, resíduos e ledgers radiativo/heating por célula; esses campos participam de commit atômico, migração, regrid conservativo, rollback/restart e export HDF5. Os novos limiares opcionais de erro/resíduo alimentam promoção/demora de rebaixamento com histerese. Ainda falta a decisão conjunta com custo/MPI/RAM de V02. |
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
| A03 | **PARCIAL** | Massa, Jeans, jump/shock, cooling length, metal/species/vorticity, forced event rows e telemetria softening existem em `amr_production_regrid.cpp`/`runAutomaticAmrRegrid`. Falta política operacional de erro de força/resolução que coordene esses sinais além da contagem de regime não qualificado. |
| A04 | **PARCIAL** | Previsão de evento/posição e forced refinement estão conectados; orçamento opcional e memória governada também. Falta função concreta benefício/cell-step/MPI/churn com histerese e decisão conjunta V02. Previsão de host driftado fora do rank ainda é limitação documentada. |
| G01 | **CONDICIONAL À VALIDAÇÃO** | TreePM/KDK e auditoria AMR-softening existem. A especificação condiciona uma correção gravitacional multinível adicional ao resultado de qualificação futura; esta rodada não executou essa qualificação e não inventou novo solver. |
| T01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | `SourceSplitPhase`, ordem de nove dispatches em `time_integration.cpp`, inversão determinística de tempo cosmológico, `makeSourceInterval` por metade e eventos somente post. Rollback global restaura estado, schedulers, force endpoint e checkpoint privado da fonte. |
| T02 | **PARCIAL** | `advanceProductionAmrHydroRecursiveImpl` agora invoca `SourceRuntimeImpl::executeAmrLevelSmoothSources` em ambas as metades de cada subintervalo, inclusive nos ranks sem linhas locais; termoquímica/dust/SGS usam o relógio `dt_L`, entram nos históricos temporais e o shadow/controller é agregado uma vez por passo com rollback privado. A autoridade de eventos permanece O(1) na sidecar populacional e ganhou cursor absoluto persistente com rejeição de gap/replay e migração de schema. SF/SSP/feedback/ventos/difusão ainda usam a barreira global enquanto faltam depósitos cross-level e difusão no grafo multinível. |
| M01 | **PARCIAL** | Governador, pool por nó MPI, reservas de pico, candidatos AMR e rollback governado existem. `GasCompositionSidecar` continua denso para todas as linhas; não há mapa esparso/tiled da química detalhada nem reserva coletiva completa por nó para todas as extensões adaptativas. |
| M02 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | CPU scalar/SIMD/OpenMP/MPI, decomposição/migração e balanceamento ponderado permanecem nos runtimes existentes. GPU segue opcional e não foi transformada em dependência. |
| V01 | **IMPLEMENTADO EM FONTE, NÃO TESTADO** | Contratos de configuração/provenance, ledgers, snapshots HDF5, restart versionado e quatro decks científicos permanecem conectados em `config.cpp`, `snapshot_hdf5.cpp`, `restart_checkpoint.cpp` e `configs/chui_part3`. Aceitação end-to-end foi deliberadamente adiada. |
| V02 | **PARCIAL** | Bias químico por shadow e orçamento AMR por headroom existem, mas são controladores separados. Faltam estimadores tipados persistentes e coordenador conjunto de erro local, benefício, custo, MPI/churn e RAM, com referência fixa e rollback. |

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

## Pendências obrigatórias verdadeiras

Não é correto declarar fechamento integral H01–V02 nesta árvore. Permanecem
obrigatórias as implementações H01 remoto, C03 histórico/erro, I02 produção
conservativa, T02 multinível bariônico, M01 química esparsa e o coordenador
A03/A04/V02. C02 requer primeiro a especificação/dados verificáveis do
inventário detalhado; S02/F02 requerem tabelas físicas de produção. G01 é
deliberadamente condicional à campanha de validação prevista pela norma.

X01 permanece fora de escopo: nenhuma implementação de MHD, moving mesh, SPH,
CHIMES, Grackle, RT multifrequência ou nova família AGN foi introduzida.
