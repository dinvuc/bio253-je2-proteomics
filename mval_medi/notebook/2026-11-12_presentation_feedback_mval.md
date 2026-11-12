# 2026-11-12 – presentation feedback into the workflow (group mval)

## Question for today
How do we build the feedback from our project presentation (Jonas) into the
workflow, and can we get peptide-level evidence without the msstats peptide table?

## What we did
- branch `feedback-presentation-mval`
- mapped every feedback note to a workflow step (new section "Feedback from the
  project presentation" in `PROJECT_WORKFLOW_mval.md`)
- searched our existing files for peptide-level data (the msstats peptide table
  will not be provided)
- searched the search FASTA headers for transposon proteins

## Results (key numbers)
- `DE_total.xlsx`, sheet `normalized_abundances`: column `nr_children_protein_Id`
  = peptides quantified per protein per sample (out of `nrPeptides`). First look:
  the 20-peptide protein has 9–15 peptides in anc/stk and 0–2 in stp; the
  37-peptide protein has 27–32 in anc, 0–1 in stk and 1–2 in stp. Which is Stp1
  and which is Stk1 must be confirmed by protein_id in the Step 3 script.
- Stk1 phosphopeptides (STYH/STY site tables): sites at positions 159–315, all
  downstream of the stk insertion (codon 46). Most are absent from stk, but
  `HIQDSVPNVTtDVR…` (~231–242) is detected 3/3 per medium in stk. Not yet
  interpreted: uniqueness, localisation and intensity to be checked.
- Stp1 has no phosphosites in our tables.
- `JE2_proteins.fasta` (4,155 entries): no ErmB, no GFP, only the genome's own IS
  transposases, so transposon products cannot be pulled from our data.

## Problems / surprises
- Peptide sequences/positions exist only for phosphopeptides, not for the total
  proteome: we can say how many Stp1 peptides are seen, not which part.
- Two note words were unclear: "tsp" (read as TSB) and "emre" (read as ermB);
  "stk1 broadly involved" needs confirmation from Jonas.

## AI use
Claude Code: read the notes, searched the files, wrote the workflow changes.
The first-look numbers above come from a quick console check; Step 3 must
reproduce them in a script before they count.

## Next steps
See "Next actions" in `PROJECT_WORKFLOW_mval.md`: Step 3 script first
(peptide counts + Stk1 phosphopeptide map), then Steps 6–7, 10, poster.
