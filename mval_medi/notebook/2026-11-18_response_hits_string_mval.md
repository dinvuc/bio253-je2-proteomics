# 2026-11-18 – Fig 2 proteins: function, STRING network, kdpD in stk vs stp (group mval)

## Question for today
What are the proteins in `step06_responses/fig_response_vs_anc_mval.png`, why might they
change, and is the kdpD result the same in both mutants? Make a STRING network of them.

## What we did
- new branch `response-hits-biology-mval`; new script `scripts/06d_response_hits_string_mval.R`
  (plan in its header), output in `mval_medi/results/step06_responses/hits_biology/`
- annotation from `read_protein_metadata()` (all 7 orthologs "good")
- STRING v12.5 API, *S. aureus* NCTC 8325 (taxon 93061): partners (score ≥ 0.7, 6 each)
  + network among them; answers cached in the folder; nodes mapped back to JE2 via `up_kegg`
- kdpD: raw values of the whole kdp operon (`plot_protein()`), peptides per sample, and a
  peptide-matched SN − TSB comparison (same peptide, same strain)
- literature: PubMed (abstracts only), DOIs in `hits_biology_mval.md`

## Results (key numbers)
- STRING: 30 nodes, 64 edges; only one edge among the 7 proteins (arcA–argF, 0.995);
  5 separate modules + B7H15_06000 alone
- kdpD protein, SN response: `anc` +1.03, `stk` −4.09 (interaction FDR 1e-10, lower in SN in
  `stk`), `stp` not testable (1/3); TSB: `stk` +1.35 vs `anc` (whole kdp region higher in
  `stk`, not in `stp`)
- kdpD peptide evidence: 12–27 peptides/sample in TSB, 0–3 in SN in all strains; 5 % kept
  (rank 9/1,121; median protein 100 %); depth equal (TSB 32,439 / SN 30,788 peptides)
- peptide-matched SN − TSB: kdpD `stk` median −4.03 (n = 4); `anc` +1.27 (n = 4, two up,
  two down); kdpE `anc` −4.09, `stk` −3.06, `stp` −4.90
- Hld (Agr/RNAIII readout, ortholog "check") falls in SN in all strains: `anc` −4.52
  (FDR 3e-5), `stk` −3.54, `stp` −2.63
- partners: betA moves with B7H15_14515 in `stk` (delta −1.76, not a hit); second arc
  copy B7H15_00350 does NOT show the `stk` rise of arcA
- all estimates observed (our limma treat model), n = 3

## Problems / surprises
- the kdpD `stk`-specific interaction is most likely a summarisation artefact: the protein
  value in SN rests on 1–3 peptides, and the `anc` "rise" comes from 2 inconsistent peptides
- JE2 has a second arcA/argF/arcC copy (B7H15_00330–00350) mapping to the same STRING node;
  rule written in the script (argF was a tie)
- `timeout` + Rscript and `glimpse()` with redirected output crashed under Git Bash
  (segfault); avoided

## AI use
Claude Code: wrote and ran the script, searched PubMed, wrote `hits_biology_mval.md` and
the figure guide. The biology is interpretation (rules section 10). The numbers and the
cited papers were marked as checked by the group on 22 Nov 2026 (the AI itself read only abstracts).

## Part 2: STRING network from the course's own JE2 organism
- question: can the network come from the course JE2 organism on STRING instead of the
  NCTC 8325 reference strain?
- found: the course organism `STRG0B09NBK` ("JE2_Bio253_HS26") answers the STRING API and
  uses our `protein_key` as its protein id; the older `STRG0A35SHG` (source of the local
  files in `data/raw/resources/`) does not answer the API
- changed `scripts/06d_response_hits_string_mval.R` (branch `string-je2-organism-mval`):
  query by protein_key, join nodes directly (ortholog mapping and the arc-copy rule
  removed); layout now one circle per hit module, dashed lines between modules
- result: 37 proteins, 70 links (score ≥ 0.7), all 37 found in our metadata as CDS;
  6 modules, only arcA–argF linked directly (0.918), 6 links between modules via shared
  partners (14470, 11190, 12195/12200); new partner of kdpD: Y1819 (B7H15_09315,
  `stk` delta −1.12, not a hit); kdpD/peptide results unchanged
- the NCTC 8325 numbers in "Results" above (30 nodes, 64 edges, 5 modules) are replaced
  by these

## Next steps
- decide if the poster says "kdp switched off by PASN in all strains" instead of
  "kdpD stk-specific"
- check further Agr targets in the data (PSMs, other RNAIII-dependent proteins)
- Xue 2011 (KdpDE–Agr), Novick 1993 (hld in RNAIII), Makhlin 2007 (ArcR): marked as checked by the group on 22 Nov 2026
