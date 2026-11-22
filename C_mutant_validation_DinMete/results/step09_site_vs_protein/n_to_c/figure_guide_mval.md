# Step 9b: N-to-C expression plots of all 25 proteins with a phosphosite hit

Made by [`C_mutant_validation_DinMete/09b_n_to_c_plots_mval.R`](../../../09b_n_to_c_plots_mval.R) in the style
of Fig. 1 of Wolski, Dittmann, Panse, Kunz & Grossmann, Chapter 8 (ISBN 978-1-0716-5519-1).
Method details are in the [Step 9 figure guide](../figure_guide_mval.md#n-to-c-expression-plots-n_to_c-added-17-nov).
The refined, paper-style subset is in [`../n_to_c_refined/`](../n_to_c_refined/).

**Description:** one figure per protein, with one panel per contrast in which the protein
has a site hit (1–3 panels, shared y-scale). x = amino-acid position from N to C;
light-yellow rectangle = protein log2FC in the same contrast (total proteome); vertical
lines = phosphosite log2FC (STY search), colour = residue; red `*` / `**` = site hit
FDR < 0.05 / < 0.01; solid = localised, dashed = not localised; grey tick = site too
sparse to test.

**How to read:** a line ending inside the rectangle = the site changes no more than its
protein (abundance-driven). A line reaching clearly beyond it = the site changes beyond
its protein (a "usage" change). **Signs:** "SN-resp" = TSB − SN, so negative = higher in SN;
"stk vs anc" etc. = first minus second, measured in TSB.

**Key points:** of 25 proteins, 11 have at least one hit that passes the refinement (Step 9c:
leave-one-out 18/18, also a hit in STYH, no QC flag). 9 of these 11 show a robust change
beyond the protein, all in the PASN response, in central-metabolism, stress and signalling proteins
(fda, pgi, eno, pflB, glyA, qoxB, traP, B7H15_04690, B7H15_04940). Stk1 and Y2093
follow their protein. 14 proteins have no robust hit.

**Overall conclusion:** PASN changes the phosphorylation of a group of glycolytic and
metabolic enzymes beyond their protein amount, in all three strains. Between strains the
only robust phospho change is the Stk1 site in `stp`, and it follows the Stk1 protein loss.

| figure | protein | hits / refined | key points (site: contrast, site log2FC vs protein log2FC, category) | conclusion |
|---|---|---|---|---|
| `fig_n_to_c_B7H15_00775_deoB_mval.png` | deoB (phosphopentomutase) | 5 / 0 | T85 SN-resp stk -5.3 vs protein -0.6, beyond protein; T87 SN-resp stk -5.0 vs protein -0.6, beyond protein; T85 stk vs anc -3.8 vs protein +0.1, beyond protein; T85 stp vs stk +5.1 vs protein +0.0, beyond protein; T87 stp vs stk +4.7 vs protein +0.0, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_01205_pflB_mval.png` | pflB (formate acetyltransferase) | 2 / 1 | S742 SN-resp anc -4.9 vs protein -3.7, beyond protein, refined; S742 SN-resp stp -4.3 vs protein -3.9, follows protein | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_02665_rplY_mval.png` | rplY (50S ribosomal protein L25/general str...) | 1 / 0 | S67 SN-resp stk -1.8 vs protein -0.7, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_02995_tuf_mval.png` | tuf (elongation factor Tu) | 2 / 0 | T300 SN-resp stk -3.2 vs protein -0.7, beyond protein; T300 SN-resp stp -2.3 vs protein -0.6, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_04345_nogene_mval.png` | B7H15_04345 (aldehyde dehydrogenase) | 1 / 0 | S210 SN-resp stk -1.8 vs protein -1.2, follows protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_04365_eno_mval.png` | eno (enolase) | 1 / 1 | S373 SN-resp stk -2.6 vs protein -0.4, beyond protein, refined | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_04690_Y845_mval.png` | Y845 (CsbD family protein) | 10 / 4 | T52 SN-resp stk +4.4 vs protein +3.8, follows protein, refined; T17 SN-resp stk +4.0 vs protein +3.8, follows protein, refined; T22 SN-resp stk +3.2 vs protein +3.8, follows protein, refined; T22 SN-resp stp +3.8 vs protein +2.8, beyond protein, refined; T22 SN-resp anc +2.6 vs protein +2.2, follows protein; S5 SN-resp stk +3.7 vs protein +3.8, follows protein; T52 SN-resp stp +4.9 vs protein +2.8, beyond protein; T17 SN-resp stp +4.5 vs protein +2.8, beyond protein; S5 SN-resp stp +4.5 vs protein +2.8, beyond protein; S37 SN-resp stp +4.0 vs protein +2.8, unresolved | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_04940_nogene_mval.png` | B7H15_04940 (NAD-specific glutamate dehydrogenase) | 4 / 3 | S264 SN-resp anc +2.2 vs protein +1.3, follows protein, refined; S264 SN-resp stk +1.9 vs protein +1.0, follows protein, refined; S264 SN-resp stp +2.9 vs protein +1.3, beyond protein, refined; T269 SN-resp stp +10.8 vs protein +1.3, beyond protein | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_04970_pgi_mval.png` | pgi (glucose-6-phosphate isomerase) | 3 / 2 | T143 SN-resp anc -2.4 vs protein -0.5, beyond protein, refined; T143 SN-resp stk -2.5 vs protein -0.6, beyond protein, refined; T143 SN-resp stp -2.1 vs protein -0.7, beyond protein | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_05485_qoxB_mval.png` | qoxB (cytochrome ubiquinol oxidase subunit I) | 2 / 1 | T552 SN-resp anc -2.9 vs protein -0.5, beyond protein, refined; T552 SN-resp stk -2.0 vs protein -0.9, beyond protein | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_05615_nogene_mval.png` | B7H15_05615 (phosphocarrier protein HPr) | 2 / 0 | S52 SN-resp stk +3.5 vs protein +0.9, beyond protein; T20 SN-resp stp +2.2 vs protein +0.8, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_05950_trxA_mval.png` | trxA (thiol reductase thioredoxin) | 3 / 0 | S71 SN-resp stk +3.3 vs protein +1.8, beyond protein; S71 SN-resp stp +2.8 vs protein +1.6, unresolved; S62 SN-resp stp +2.7 vs protein +1.6, unresolved | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_06355_nogene_mval.png` | B7H15_06355 (serine/threonine protein kinase) | 2 / 2 | T164 stp vs anc -4.2 vs protein -5.4, unresolved, refined; T166 stp vs anc -4.1 vs protein -5.4, unresolved, refined | robust, but the site **follows its protein** (abundance-driven) |
| `fig_n_to_c_B7H15_07545_nogene_mval.png` | B7H15_07545 (scaffolding protein) | 1 / 0 | S7 SN-resp stk +2.4 vs protein +1.9, follows protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_08315_nogene_mval.png` | B7H15_08315 (phosphogluconate dehydrogenase (NADP(...) | 3 / 0 | S34 SN-resp anc -2.5 vs protein -0.2, beyond protein; S34 SN-resp stk -3.2 vs protein -0.4, beyond protein; S34 SN-resp stp -2.2 vs protein -0.4, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_08885_nogene_mval.png` | B7H15_08885 (hypothetical protein) | 1 / 0 | Y32 SN-resp stk +2.0 vs protein +2.8, follows protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_09165_thrS_mval.png` | thrS (threonine--tRNA ligase) | 2 / 0 | S30 SN-resp anc -2.4 vs protein -1.9, follows protein; S30 SN-resp stp -2.6 vs protein -1.6, follows protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_09250_pyk_mval.png` | pyk (pyruvate kinase) | 1 / 0 | T537 stp vs anc +2.0 vs protein +0.1, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_09315_Y1819_mval.png` | Y1819 (universal stress protein UspA) | 1 / 0 | T157 interaction stp vs stk -2.3 vs protein -1.2, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_10065_traP_mval.png` | traP (signal transduction protein TRAP) | 2 / 1 | S150 SN-resp anc -2.7 vs protein +0.2, beyond protein, refined; S150 SN-resp stp -2.6 vs protein -0.2, beyond protein | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_10475_Y2093_mval.png` | Y2093 (hypothetical protein) | 2 / 1 | S51 SN-resp stk +4.0 vs protein +3.1, follows protein, refined; S51 SN-resp stp +2.3 vs protein +2.0, follows protein | robust, but the site **follows its protein** (abundance-driven) |
| `fig_n_to_c_B7H15_11735_glyA_mval.png` | glyA (serine hydroxymethyltransferase) | 1 / 1 | Y51 SN-resp stk +4.6 vs protein +0.7, beyond protein, refined | **robust change beyond the protein** (usage change) |
| `fig_n_to_c_B7H15_11980_glmM_mval.png` | glmM (phosphoglucosamine mutase) | 1 / 0 | S102 SN-resp stk -2.2 vs protein -0.1, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_13850_pgcA_mval.png` | pgcA (phosphoglucomutase) | 2 / 0 | S193 SN-resp stk -4.8 vs protein -0.7, beyond protein; T191 SN-resp stk -3.8 vs protein -0.7, beyond protein | hits **not robust** (failed a refinement check): do not interpret alone |
| `fig_n_to_c_B7H15_14480_fda_mval.png` | fda (class I fructose-bisphosphate aldolase) | 7 / 3 | S238 SN-resp anc -3.7 vs protein -1.4, beyond protein, refined; S238 SN-resp stk -3.3 vs protein -1.0, beyond protein, refined; S238 SN-resp stp -2.9 vs protein -1.4, beyond protein, refined; S25 SN-resp anc -2.3 vs protein -1.4, follows protein; S22 SN-resp stk -2.2 vs protein -1.0, beyond protein; S25 SN-resp stk -2.0 vs protein -1.0, follows protein; S22 SN-resp stp -2.4 vs protein -1.4, beyond protein | **robust change beyond the protein** (usage change) |
