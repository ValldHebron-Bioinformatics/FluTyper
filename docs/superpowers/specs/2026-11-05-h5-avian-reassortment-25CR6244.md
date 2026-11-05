# Spec: avaluació de reassortment de 25CR6244 (H5 LPAIV) en context europeu H5 2025

**Data:** 2026-11-05
**Mostra d'interès:** `25CR6244` (Barcelona, 2025-11-26), notificada com a H5N1 LPAIV, no-Gs/GD.

## 1. Pregunta científica

**Objectiu (reformulat el 2026-11-05):** caracteritzar **d'on ve 25CR6244** — origen i
parentiu dins la diversitat H5 europea de juny–novembre 2025. Per a cada segment avaluable:
on cau 25CR6244 a l'arbre i amb quin suport de bootstrap; els seus 5–10 parents més propers
(país, hoste, data, genotip); si el clade on cau és geogràficament coherent (els LPAI del
conjunt són només de Suècia, Itàlia i Països Baixos) o cau en un lloc inesperat; i si HA, NS i
NP apunten al mateix origen geogràfic i d'hoste. Si convergeixen, la lectura d'origen és
sòlida; si divergeixen, es reporta com a observació.

**Reassortment: NO contestable amb aquesta mostra.** La pregunta original — si algun segment
s'agrupa amb virus HPAIV Gs/GD 2.3.4.4b — exigeix comparar segment a segment tot el genoma.
Amb la cobertura real de §2, 25CR6244 només té material interpretable en 3 dels 8 segments
(HA, NS i, com a exploratori, NP). Amb 3 segments de 8 **no es pot afirmar ni descartar
reassortment**, i l'informe ho ha de dir així; qualsevol lectura que ho suggereixi seria
excedir-se. Els segments que falten (PB2, PB1, PA, NA, MP) són el que caldria resequenciar.

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

**Cobertura real de 25CR6244** (definició de `PhyloCoverageQC`: posicions de la referència
cobertes per la mostra, no N ni forat, sobre la longitud de la referència; referència H5
`A/American_wigeon/South_Dakota/23-033792-006-original/2023`). Els registres estan farcits
amb N i guions fins a 2341 bp, cosa que amaga la mancança real:

| Segment | Cobertes / long. ref. | Cobertura | Estatus |
|---|---|---|---|
| NS | 838 / 838 | 100,00% | primari (llindar 0,95) |
| HA | 1695 / 1704 | 99,47% | primari (llindar 0,95) |
| NP | 1407 / 1497 | 93,99% | **secundari / exploratori** (llindar 0,93 declarat a part; per sota del primari) |
| PB1 | 258 / 2274 | 11,35% | no avaluable |
| PB2 | 238 / 2280 | 10,44% | no avaluable |
| MP | 65 / 982 | 6,62% | no avaluable |
| PA, NA | — | absents | no avaluable |

El llindar primari **no es mou** per encabir NP: triar-lo mirant el resultat de la mostra
d'interès no és defensable. NP s'executa a part, etiquetat com a exploratori; una col·locació
interessant demanaria confirmació independent abans d'afirmar res. No s'executen PB1, PB2 ni MP.

Sense NA no es pot confirmar el subtipus N1 des de la seqüència; el N1 ve de PCR i es
documenta com a procedent del laboratori, no de l'anàlisi.

**Pathotype:** el de 25CR6244 és LPAI llegit del motiu de tall d'HA (`RNVPQRETR/GLF`, 2
residus bàsics en P1–P5), no heretat del sufix de la referència del nextclade.

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
