# Spec: avaluació de reassortment de 25CR6244 (H5 LPAIV) en context europeu H5 2025

**Data:** 2026-11-05
**Mostra d'interès:** `25CR6244` (Barcelona, 2025-11-26), notificada com a H5N1 LPAIV, no-Gs/GD.

## 1. Pregunta científica

25CR6244 és un H5 LPAIV que **no pertany al linatge Gs/GD**, però durant la tardor
els H5 LPAIV i els H5 HPAIV del clade 2.3.4.4b co-circulen en poblacions d'aus
silvestres que se solapen (gavians, *Larus* spp., i anàtids). El reassortment pot
afectar **qualsevol segment**. Cal determinar, **segment a segment**, si algun
segment de 25CR6244 s'agrupa amb virus HPAIV Gs/GD en lloc del seu context LPAIV.

Un segment de 25CR6244 que caigui dins d'un clade HPAIV 2.3.4.4b amb suport de
bootstrap és evidència de reassortment amb virus HP Gs/GD. Aquest és el resultat
que l'anàlisi ha de poder afirmar o descartar per a cada segment.

## 2. Dades d'entrada (ja al repositori)

| Fitxer | Contingut |
|---|---|
| `docs/fastas/cresaid.fasta` | 21 mostres pròpies, 156 registres. 25CR6244 té **6 segments**: HA, PB2, PB1, NP, MP, NS |
| `docs/fastas/AVIAN_JUNY_NOV_2025.fasta` | 10.829 registres GISAID, capçalera `EPI_ISL|SEG|DATE|STRAIN|LAB` |
| `docs/AVIAN_JUNY_NOV_2025.xls` | Metadades GISAID, 1.369 aïllats, 2025-06-04 → 2025-11-30 |
| `docs/metadata.csv` | Metadades pròpies en format FluTyper: `ID, DATE, LOCATION` |

Composició del conjunt GISAID (camp `Pathogenicity`):

- **HPAI:** 1.327, tots clade 2.3.4.4b. UK 491, DE 117, NL 104, IT 98, CZ 87, ES 67, BE 62, SE 52, PL 34, DK 32, NO 31, AT 27, FI 23, PT 23, IE 19, …
- **LPAI:** 33 — H5N3 ×19, H5N4 ×8, H5N2 ×3, H5 ×2, H5N1 ×1. Només Suècia (24), Itàlia (7), Països Baixos (2). Hostes: *Mallard* 23, *Duck* 5, *Teal* 2, altres anàtids 3. **Cap *Larus* spp.**
- **Patogenicitat en blanc:** 9 — s'han de classificar amb la crida pròpia de FluTyper.

**Limitació acceptada:** 25CR6244 no té PA ni NA en les dades disponibles. L'avaluació
de reassortment cobreix 6 dels 8 segments. Sense NA no es pot confirmar el subtipus
N1 des de la seqüència; el subtipus notificat es documenta com a procedent del
laboratori, no de l'anàlisi.

## 3. Disseny dels datasets

Dos datasets de context, **per segment**, cadascun amb 25CR6244 sempre inclosa
com a mostra de referència.

### 3.1 Context LPAIV
- Tots els H5 LPAIV europeus del conjunt (juny–novembre 2025), inclosos els 9
  de patogenicitat desconeguda que la crida pròpia determini com a LPAI.
- Tots els hostes. Es registra l'hoste de cada fulla per poder destacar
  *Larus* spp. (*L. argentatus*, *L. fuscus*, *L. ridibundus*, *L. michahellis*) i
  anàtids (ànecs, oques, cignes); cap d'aquests dos grups es filtra.
- Mida esperada: ~34–42 fulles per segment.

### 3.2 Context HPAIV
- **Tots** els aïllats amb clade `2.3.4.4b` (1.327), sense submostreig.
- Això cobreix per construcció el nord, centre i oest d'Europa i tots els genotips
  circulants de la temporada; el camp `Genotype` del GISAID (GenoFLU / Genin / ggFLU)
  s'arrossega com a anotació de cada fulla per poder-los identificar a l'arbre.

### 3.3 Període
Juny–novembre, que cobreix el període post-reproductor i la migració de tardor —
quan es barregen poblacions d'aus del nord, centre i oest d'Europa — i coincideix
amb l'augment de deteccions d'H5.

## 4. Sortides requerides

1. **Cribratge:** taula per aïllat amb subtipus inferit, patotip (LPAI/HPAI), clade
   i genotip, amb el patotip propi de FluTyper **comparat** amb el camp
   `Pathogenicity` del GISAID, i els desacords llistats.
2. **Arbres per segment:** HA, PB2, PB1, NP, MP, NS, per a cada dataset, amb suport
   de bootstrap i 25CR6244 identificable.
3. **Exportació Auspice:** un dataset Nextstrain per arbre, amb mapa d'Europa i
   coloracions per patotip, clade, genotip, hoste i país.

## 5. Restriccions

- El patotip no s'infereix mai del nom de la soca ni del camp del GISAID quan la
  seqüència permet calcular-lo; els dos valors es reporten per separat.
- Cap taxa de rellotge molecular, data o coordenada s'inventa. Una taxa sense font
  citable fa que l'arbre s'exporti només en divergència, amb avís — el comportament
  que el mòdul ja té per a un subtipus sense fila a `auspice_clock_rates.tsv`.
- Els 9 aïllats de patogenicitat desconeguda es documenten explícitament al resultat,
  classificats o no.

## 6. Estat actual del codi (el que cal estendre)

| Lloc | Estat | Cal |
|---|---|---|
| `Phylogenetics.nf:622` | `--phylogenetics is only available for the HUMAN protocol` | permetre `AVIAN` |
| `Phylogenetics.nf:807` `phyloKnownSubtypes()` | `['H3N2','H1N1pdm09']` | afegir **un sol** subtipus `H5` |
| `Phylogenetics.nf:826` `phyloSubtypeTag()` | mapa de 2 entrades | afegir `H5`, que fa *match* **només pel tag H** |
| `Phylogenetics.nf:63,80` | filtre `"${h_tag}${n_tag}" == subtype_tag` | un tag sense component N compara només `h_tag` |

El subtipus és `H5`, no una llista d'`H5Nx`: agrupar per H+N partiria el context LPAIV en
quatre arbres (H5N1, H5N2, H5N3, H5N4) i deixaria fora els 2 aïllats que al GISAID són
`A / H5` perquè el N no s'ha pogut cridar. Un cribratge de reassortment vol tots els H5 del
període al mateix arbre per segment, i el subtipus N de cada fulla queda com a anotació.
| `Phylogenetics.nf:813` `phyloSubtypeDefaults()` | refs humanes | referència/outgroup H5 |
| `Phylogenetics.nf:38` | llegeix `${protocol}_references.fasta` | el fitxer AVIAN es diu `Avian_references.fasta` |
| `PhyloPreflight.nf:183` | accepta prefixos CDS `{subtype, ref_subtype}` | AVIAN usa `H5_…` i `N1_…` |
| `main.nf:425-470` | blocs `if ('H3N2' in …)` / `if ('H1N1pdm09' in …)` hardcodejats | instanciació genèrica per subtipus |
| `MergeReports.nf:26` `PHYLO_KNOWN_SUBTYPES` | `["H1N1pdm09","H3N2"]` | mantenir sincronitzat |
| `RESOURCES/auspice_clock_rates.tsv` | només H3N2 i H1N1pdm09 | files H5N1 citades, o divergència-only |
| `RESOURCES/auspice_config.json` | coloracions humanes (age_group, sex) | patotip, clade, genotip, hoste |
