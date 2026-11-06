# Screening of 1,390 H5 samples (2026-11-05): HPAI / LPAI from the HA cleavage site

Run `deadly_banach` (Nextflow, started 2026-11-05 14:07:22, `--protocol AVIAN`, input 10,978 records / 1,390 samples = 1,369 GISAID isolates from Europe 2025-06-01..11-30 + the 21 own samples of `docs/fastas/cresaid.fasta`, 25CR6244 de-duplicated). Raw tables in this folder: `2026-11-05-screening_comparison.csv` (one row per sample) and `2026-11-05-ha_cleavage.csv` (cleavage read with reason). Dataset construction: `2026-11-05-dataset_report.txt`.

## Which run this is

The recovered pair `screening_comparison2.csv` + `ha_cleavage2.csv` was chosen by content, not by suffix. Both comparison files have 1,390 rows, so the row count does not discriminate. What does: the unsuffixed pair is the earlier run (`nostalgic_magritte`, started 13:43:16), after the `sampleIdFromHeader` fix (commit 2d06288, 13:43:07) but **before** the `seqkit grep -f` fix (ce46124, 14:07:06). In that run 139 samples lost their HA1 (commas/quotes in headers made `seqkit grep -p` drop segments): `ha_cleavage.csv` has 1,385 rows, 147 indeterminate (139 of them "no HA1"), and 139 `only_gisaid` rows. The "2" pair comes from `deadly_banach`, started 14:07:22, 16 seconds after the second fix: 1,390 rows, 14 indeterminate, 1 `only_gisaid`. The 138 samples whose call differs between the two files are the ones that run lost. File mtimes are not usable evidence (the recovered copies were re-written on 2026-11-06). The original `RESULTS/screening*` folders no longer exist, so these CSVs are the only record.

## Counts

Cleavage-site read (`pathotype_cleavage`), 1390 samples: **HPAI 1335**, **LPAI 41**, **indeterminate 14**.

Rule: HPAI = at least 3 basic residues in P1-P5 or a furin motif; LPAI = a readable monobasic site; anything unreadable is `indeterminate`, never LPAI.

## The three pathotype columns against `pathotype_gisaid`

`pathotype_cleavage` is read from the HA protein of the sample. `pathotype_nextclade` is only the HPAI/LPAI suffix of the nearest reference found by `nextclade sort` (a nearest-neighbour label, H5/H7 only), not a cleavage read. `pathotype_gisaid` is the submitter's field; it is blank for the 21 own samples and for 9 GISAID isolates.

| pathotype_cleavage \ pathotype_gisaid | HPAI | LPAI | (blank) | Total |
|---|---|---|---|---|
| HPAI | 1326 | 0 | 9 | 1335 |
| LPAI | 0 | 33 | 8 | 41 |
| indeterminate | 1 | 0 | 13 | 14 |
| **Total** | 1327 | 33 | 30 | 1390 |

| pathotype_nextclade \ pathotype_gisaid | HPAI | LPAI | (blank) | Total |
|---|---|---|---|---|
| (blank) | 0 | 0 | 13 | 13 |
| HPAI | 1327 | 0 | 13 | 1340 |
| LPAI | 0 | 33 | 4 | 37 |
| **Total** | 1327 | 33 | 30 | 1390 |

`pathotype_gisaid` itself: HPAI 1,327, LPAI 33, blank 30 (21 own samples + 9 GISAID isolates).

**Agreement** (cleavage vs GISAID): agree 1359, **disagree 0** (no sample is HPAI by one and LPAI by the other), only_cleavage 17 (GISAID blank, cleavage read gives a call), only_gisaid 1 (cleavage indeterminate, GISAID has a call), neither 13.

## Every row that is not a plain agreement (31)

| ID | inferred_subtype | cleavage_motif | cleavage | nextclade | GISAID | clade | agreement |
|---|---|---|---|---|---|---|---|
| 24CR5067 | H5N2 | RNVPQRETR/GLF | LPAI | LPAI | - | - | only_cleavage |
| 24CR5076 | H13Nx | RNIPAISNR/GLF | LPAI | - | - | - | only_cleavage |
| 24CR5124 | H13N2 | - | indeterminate | - | - | - | neither |
| 24CR5197 | H13N2 | - | indeterminate | - | - | - | neither |
| 24CR5236 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 24CR5247 | H13N2 | RNVPAISNR/GLF | LPAI | - | - | - | only_cleavage |
| 25CR5363 | H4N6 | RNIPEKASR/GLF | LPAI | - | - | - | only_cleavage |
| 25CR5413 | H6N8 | RNVPQIENR/GLF | LPAI | - | - | - | only_cleavage |
| 25CR5659 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR5713 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR5714 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR5715 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR5824 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR5865 | H7N7 | NVPEIPKGR/GLF | LPAI | LPAI | - | - | only_cleavage |
| 25CR5922 | H13N6 | RNVPAISNR/GLF | LPAI | - | - | - | only_cleavage |
| 25CR5923 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR6062 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR6084 | H5N1 | PLREKRRKR/GLF | HPAI | HPAI | - | - | only_cleavage |
| 25CR6193 | H9Nx | - | indeterminate | - | - | - | neither |
| 25CR6229 | H6Nx | - | indeterminate | - | - | - | neither |
| 25CR6244 | H5Nx | RNVPQRETR/GLF | LPAI | LPAI | - | - | only_cleavage |
| EPI_ISL_20212724 | H5N2 | - | indeterminate | LPAI | - | 2.2.1 | neither |
| EPI_ISL_20246482 | H5N1 | - | indeterminate | HPAI | - | 2.3.4.4b | neither |
| EPI_ISL_20282939 | H5N1 | - | indeterminate | HPAI | - | 2.3.4.4b | neither |
| EPI_ISL_20289024 | HxNx | - | indeterminate | - | - | - | neither |
| EPI_ISL_20289027 | HxNx | - | indeterminate | - | - | - | neither |
| EPI_ISL_20290821 | H5N1 | CEFRL*LKT/GLF | indeterminate | HPAI | HPAI | 2.3.4.4b | only_gisaid |
| EPI_ISL_20349763 | HxN1 | - | indeterminate | - | - | - | neither |
| EPI_ISL_20354011 | HxNx | - | indeterminate | - | - | - | neither |
| EPI_ISL_20505315 | H5N1 | - | indeterminate | HPAI | - | 2.3.4.4b | neither |
| EPI_ISL_20505324 | H5Nx | - | indeterminate | HPAI | - | 2.3.4.4b | neither |

## The 14 indeterminate samples and why

| ID | motif | reason | GISAID |
|---|---|---|---|
| 24CR5124 | - | HA2 does not start with the conserved GLF anchor | - |
| 24CR5197 | - | HA2 does not start with the conserved GLF anchor | - |
| 25CR6193 | - | HA2 does not start with the conserved GLF anchor | - |
| 25CR6229 | - | no HA1 protein | - |
| EPI_ISL_20212724 | - | HA2 does not start with the conserved GLF anchor | - |
| EPI_ISL_20246482 | - | HA2 does not start with the conserved GLF anchor | - |
| EPI_ISL_20282939 | - | HA2 does not start with the conserved GLF anchor | - |
| EPI_ISL_20289024 | - | no HA1 protein | - |
| EPI_ISL_20289027 | - | no HA1 protein | - |
| EPI_ISL_20290821 | CEFRL*LKT/GLF | unreadable residue(s) * in the cleavage motif | HPAI |
| EPI_ISL_20349763 | - | no HA1 protein | - |
| EPI_ISL_20354011 | - | no HA1 protein | - |
| EPI_ISL_20505315 | - | HA2 does not start with the conserved GLF anchor | - |
| EPI_ISL_20505324 | - | HA2 does not start with the conserved GLF anchor | - |

By reason: HA2 does not start with the conserved GLF anchor: 8; no HA1 protein: 5; unreadable residue * in the cleavage motif: 1. Four are own samples (24CR5124, 24CR5197, 25CR6193, 25CR6229); ten are GISAID isolates, of which nine have a blank GISAID pathogenicity and stay UNKNOWN (listed in the dataset report), and EPI_ISL_20290821 keeps GISAID's HPAI call because its motif has a stop codon (`*`).

## 25CR6244

HA cleavage motif **`RNVPQRETR/GLF`**, 2 basic residues in P1-P5, read as **LPAI** (monobasic cleavage site (2 basic residues in P1-P5)). `pathotype_nextclade` also says LPAI, but that is a nearest-reference label and is not independent evidence. Sequence-derived subtype: **H5Nx**. The N1 in the lab's "H5N1" comes from PCR; the NA is not sequenced.

**Known defect in the metadata of this sample.** The note "N1 by PCR, NA not sequenced" was written into the `GENOTYPE` column of 25CR6244's row (`screening_meta.csv` and the three dataset metadata files), so the screening table above shows it under `genotype`, and any tree or Auspice colouring by genotype would show it as a genotype category. It is a note, not a genotype. Not corrected here; it needs its own column (or the note removed) before the trees are re-exported.

## Why Task 3 Step 7 was skipped

Step 7 of Task 3 (write this file) and Step 8 (commit it) were never done. What the ledger shows:
- I ruled early that Steps 5-7 were long-running I/O and would run in the background while Tasks 4-9 were implemented, and that Task 3's completion line would be written only once its own verification steps passed against real output.
- When the clean screening finished, I recorded the numbers in the ledger as `Task 3: SCREENING COMPLETE (RESULTS/screening2 ...)`. That entry carried the counts, so I treated the ledger as the record and the file as a formality. The file, the commit and the `task-done` completion line for Task 3 were never produced.
- Right after that, every next step came from a new instruction (stop before building datasets, the cleavage reader, the LPAI run, the change of objective, the HPAI/combined/NP runs), and nothing in the process went back to a task with no completion line.
- The consequence was concrete: the only copies of the screening tables lived in the session scratchpad and in `RESULTS/screening*`, and the latter were later gone. The CSVs survived only because they were copied out by hand.

What prevents a repeat: a task whose plan names an output file is not closed by a ledger note; `task-done` is the only thing that writes the completion line, so the check before the final review is to list every task without a `Task N: complete` line (Tasks 3, 9, 10, 11 and 12 had none when the branch was pushed) and to confirm that each file named in a plan step exists in git. This file and the raw tables close Task 3's missing output.
