# CHUÍ — C02: rede multiiônica C, O, Si e Fe

**Estado:** dataset *completo em estrutura*, **NÃO QUALIFICADO cientificamente** (`verification unverified`). Gerado em 8 de novembro de 2026 a partir de coeficientes e tabelas **publicados**, sem parâmetros físicos inventados. A correção do parser e a validade física integral **não foram demonstradas** por compilação/testes do CHUÍ.

## Artefato principal

`chui_metal_ion_rates_v1_C_O_Si_Fe_unverified.txt`

Esquema: `chui_metal_ion_rates_v1` conforme parser `src/physics/metal_ion_network.cpp` do CHUÍ patch 52. O documento contém **58 estágios** (C:7; O:9; Si:15; Fe:27), 54 transições ionizáveis, 54 potenciais de ionização e **189 temperaturas log10(K)** de 4 a 8. Cada estágio contém as três séries exigidas pelo parser: `collisional_ionization`, `recombination` e `cooling`. Os coeficientes são log10(cm³/s), log10(cm³/s), e log10(erg·cm³/s). `-inf` significa taxa nula, utilizada para os limites físicos (núcleo completamente ionizado, átomo neutro) e para `u>80` do método publicado de Voronov.

### Fontes científicas e procedimentos

1. **Resfriamento por estágio iônico:** O. Gnat & G. J. Ferland, *Ion-by-Ion Cooling Efficiencies*, ApJS **199**, 20 (2012), DOI `10.1088/0067-0049/199/1/20`. Dados em 4 arquivos do espelho GitHub de `pikarpov/CoolingCurve`: [`C.txt`](https://github.com/pikarpov/CoolingCurve/blob/master/IonByIon/C.txt), [`O.txt`](https://github.com/pikarpov/CoolingCurve/blob/master/IonByIon/O.txt), [`Si.txt`](https://github.com/pikarpov/CoolingCurve/blob/master/IonByIon/Si.txt), [`Fe.txt`](https://github.com/pikarpov/CoolingCurve/blob/master/IonByIon/Fe.txt). Arquivos contêm 189 temperaturas comuns e CIE (coluna final). **Foi utilizada cada curva ion-by-ion, NÃO a curva CIE.** Os valores de resfriamento foram compactados à resolução de 0,001 dex (erro de quantização de até 0,0005 dex, aproximadamente 0,12%). O arquivo `cool_lzw.b64` contém as 58 séries reconstituíveis usando `decode_cooling.py`.
2. **Ionização colisional:** G. S. Voronov, *A Practical Fit Formula for Ionization Rate Coefficients of Atoms and Ions by Electron Impact*, Atomic Data and Nuclear Data Tables **65**, 1 (1997), [dados e FORTRAN originais](https://www.pa.uky.edu/~verner/col.html). Equação: `u=E/(kB*T)`, `C=A(1+P sqrt(u))u**K exp(-u)/(X+u)`, `u>80 => 0`; `kB=8.617333262145e-5 eV/K`. Os 54 registros transcritos estão em `voronov_selected.tsv`.
3. **Recombinação radiativa:** D. A. Verner & G. J. Ferland, *Atomic Data for Astrophysics I*, ApJS **103**, 467 (1996), [table1.dat](https://www.pa.uky.edu/~verner/rec.html), aplicado aos íons de sequências H, He, Li e Na quando presentes. Fórmula `A/[sqrt(T/T0)*(1+sqrt(T/T0))**(1-B)*(1+sqrt(T/T1))**(1+B)]`. Demais íons: parametrização publicada `A*(T/1e4)**(-B)` de Aldrovandi & Pequignot/Shull & van Steenberg/Arnaud & Rothenflug, conforme `pl.dat` (idem). Para Fe com 2 a 13 elétrons recombinados, uso do ajuste Arnaud & Raymond (1992) `ar1.dat`, salvo exceções com fit mais amplo de Verner/Ferland. **As parametrizações são selecionadas, NÃO somadas uma à outra.**
4. **Recombinação dieletrônica:** C/O/Si via compêndio `hdr.dat`, Arnaud & Rothenflug e outros, `A*T**(-1.5)*exp(-T0/T)*(1+B*exp(-T1/T))` (idem). Para Fe via M. Arnaud & J. Raymond, ApJ **398**, 394 (1992), `ar2.dat`, `T**(-1.5)*sum_i c_i*exp(-E_i/kBT)`. Os valores são somados a **uma** contribuição de recombinação radiativa, pois são canais físicos distintos. Não existe canal dieletrônico para a recombinação ao íon hidrogenoide de um núcleo nu.
5. **Potenciais de ionização:** os 54 limiares em eV de **Voronov 1997** (`dE` em `cfit.dat`). Valores publicados, mas **arredondados para a precisão do ajuste** — **não** são os potenciais criticamente avaliados do NIST ASD. A atualização para NIST é desejável antes de declarar o contrato energético cientificamente qualificado.

### Limitações importantes (NÃO ignorar)

- **Faltam contribuições específicas de recombinação dieletrônica de baixa temperatura**, p.ex. ressonâncias de Nussbaumer/Storey/Badnell/CHIANTI. A série `recombination` está **preenchida em todos os estágios**, mas **não significa completude de todos os canais físicos**.
- **Não contém fotoionização ou fotoaquecimento de metais.** As seis curvas HI/HeI/HeII do TREECOOL FG20 não permitem calcular diretamente taxas de fotoionização para os 54 níveis metálicos. É necessário um UVB espectral em energia/frequência, seção de choque por íon e integração da taxa e aquecimento. O CHUÍ aceita a omissão dos campos `photo`, mas isso limita a física do UVB metálico.
- **A validade de cada ajuste não cobre necessariamente todo o intervalo global de temperatura.** A grade C02 foi escolhida pela cobertura comum das curvas Gnat/Ferland; os intervalos recomendados dos demais ajustes e a dinâmica fora do equilíbrio precisam ser examinados. A interpretação científica abaixo de `10^4 K` **não está coberta** pelo resfriamento C02 e requer outro regime de física.
- Gnat/Ferland assumem plasma de baixa densidade opticamente fino. Suas curvas são eficiências **íon por íon**, não uma tabela de metal cooling de equilíbrio; podem envolver trabalho térmico de ionização, resfriamento de recombinação, linhas e bremsstrahlung. Ao acoplar ao `ionPotentialDensity` do CHUÍ, conferir cuidadosamente **ausência de dupla contagem** entre `cooling`, `metal_ion_potential_delta_erg` e a energia térmica evoluída. A expressão `energy_convention collisional_ionization_work_in_cooling` é declarada como exigida pelo parser, mas não substitui uma qualificação científica.
- Os **direitos de redistribuição de todos os dados publicados, especialmente dos suplementos Gnat/Ferland**, não foram confirmados. O campo `license` sinaliza a pendência; **não considerar autorização para republicar em repositório público**. Verifique os termos dos fornecedores científicos.
- `verification unverified` significa que o loader do CHUÍ exige opt-in explícito `physics.thermochemistry_allow_unqualified_rates=true` para habilitar estas taxas no caminho de produção. **Não desativar a proteção por conveniência**. O arquivo não foi submetido ao parser C++ ou a experimentos, de acordo com a campanha de implementação por fonte.

### Regeneração e auditoria somente de dados

A pasta contém `build_c02_from_sources.py`, `decode_cooling.py`, `cool_lzw.b64` e as tabelas de parâmetros transcritas do site de D. Verner. Para regenerar sem qualquer compilação ou execução do CHUÍ:

```bash
python3 build_c02_from_sources.py
sha256sum chui_metal_ion_rates_v1_C_O_Si_Fe_unverified.txt
```

O script verifica a gramática dos registros, limites e contagem de 58×3 séries, 54 potenciais e 189 pontos, além de registrar checksums no `C02_manifest.json`. Isto é apenas **verificação estrutural de dados**, não validação científica.

**Status apropriado:** `C02_DATA_ASSEMBLED_UNQUALIFIED` — não `C02_SCIENCE_VERIFIED`. O código CHUÍ e seus testes/simulações continuam **não executados**.
