# reference/ – public reference data

Two independent things live here:

1. **Protein metadata for the whole search database** —
   `JE2_protein_metadata.tsv.gz` (4,155 proteins x 80 columns) and
   `JE2_protein_terms.tsv.gz`. Locus tag, product, genome position, **UniProt
   accession + BLAST e-value**, gene symbol, EC number, GO / KEGG / STRING
   annotation. Load with `read_protein_metadata()`. Full column documentation:
   [README_protein_metadata.md](README_protein_metadata.md); the script that
   builds it is `make_protein_metadata.R` (it runs in the lecturer's resource
   folder, where the FASTA, GFF and STRING exports live).
2. **A proposal to extend the search FASTA** with ErmB — the rest of this file.

## Extending the search database

Small, public reference sequences. Unlike `data/`, this folder **is** committed:
nothing here is unpublished course data.

## Why

Our two mutants carry the *bursa aurealis* transposon (Nebraska Transposon
Mutant Library). It is 3,237 bp and carries two genes of its own:

* **`ermB`** – erythromycin resistance, the marker the mutants are selected with
* **`gfp`**

If either is expressed, its protein is physically present in the `stk` and `stp`
samples and absent from the wild type. That would be a **positive** control for
the transposon, which is much nicer evidence than the absence of Stk1/Stp1.

But a database search can only report sequences that are in the FASTA it
searched. Our database
(`p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta`, 4,155 entries)
contains **no ErmB at all**, so the question has never been asked of the data.
Notebook [`analyses/07_transposon_mutants.qmd`](../analyses/07_transposon_mutants.qmd)
walks through this.

## What is here

| file | what |
|---|---|
| `ermB_P20173.fasta` | ErmB, 245 aa. UniProt P20173 (*E. faecalis*) and Q799N2 (*S. aureus*) have **identical** sequences, so the exact source of the pBursa copy barely matters. |

## How to use it

```bash
cat p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta \
    reference/ermB_P20173.fasta > JE2_plus_ermB.fasta
```

then re-run the search (FragPipe/MSFragger) against `JE2_plus_ermB.fasta` and
look at ErmB across the 27 samples. The prediction is explicit:

* **detected in the 18 mutant samples, absent from the 9 `anc` samples** →
  the transposon is present and expressed; a clean positive control.
* **not detected anywhere** → inconclusive, not evidence of absence. `ermB` is
  classically regulated by translational attenuation and is strongly induced by
  erythromycin. If these cultures were grown **without** erythromycin, ErmB may
  simply be below the detection limit.

Check the growth protocol before interpreting the result either way.

## What is deliberately *not* here

**Truncated Stp1 (1–188) and Stk1 (1–45).** Tempting, but they would not help:
every tryptic peptide of a truncated protein is also a peptide of the
full-length entry already in the database, so those peptides are *already*
being found and assigned. Adding the truncations would only create redundant
protein groups and make inference harder. The way to see a truncation is the
one used in notebook 07 – look at **where along the protein** the detected
peptides fall. No new FASTA needed.

**The junction peptide.** The one genuinely new sequence would be the chimeric
peptide that runs from the last native codon into the transposon's inverted
repeat. To build it you need the exact pBursa sequence and the reading frame at
the insertion point. Worth doing if you want direct proof of the insertion in
the protein data, but it is a single low-abundance peptide, so absence would
prove little.

**GFP.** Already in the database as a contaminant entry (`zz_Y-FGCZCont00258_`,
*A. victoria* GFP) and **not detected**. Weak evidence: the transposon carries a
codon-optimised *gfp* whose protein sequence may differ from that entry. If you
re-search anyway, add the actual pBursa GFP sequence rather than trusting the
contaminant entry.

## Sources

* Fey P. D. *et al.* (2013) *mBio* 4:e00537-12 — the NTML and the *bursa
  aurealis* construct. <https://doi.org/10.1128/mBio.00537-12>
* Bae T. *et al.* (2004) *PNAS* 101:12312–12317 — the original *bursa aurealis*
  transposon and pBursa. PMID 15304642
* NTML insertion map (JE2 coordinates) — <https://ntml.unmc.edu/tninsertionmap>
* Prust N. *et al.* (2021) *Mol Cell Proteomics* 20:100034 — used the same two
  mutants, NE217 and NE1919.
