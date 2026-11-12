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
2. **The proteins carried by the transposon** (`transposon_proteins.fasta`) —
   the rest of this file.

## The proteins on the transposon

`transposon_proteins.fasta` — 11 sequences, pulled straight from the GenBank
records of the two delivery plasmids by
[`make_transposon_fasta.R`](make_transposon_fasta.R):

* **pBursa — GenBank [AY672109](https://www.ncbi.nlm.nih.gov/nuccore/AY672109)**,
  7,383 bp. GenBank annotates the transposon as `mobile_element 4147..7383`
  (3,237 bp), and exactly **two** coding sequences lie inside it:

  | protein | gene | aa | tryptic peptides (6–30 aa) |
  |---|---|---|---|
  | rRNA methylase (erythromycin resistance) | `ermB` | 245 | 16 |
  | green-fluorescent protein | `gfp` | 238 | 12 |

* Everything else is **plasmid backbone** (pBursa `cop-6`, `repC`, `cat`) or
  **helper plasmid** pFA545 ([AY672108](https://www.ncbi.nlm.nih.gov/nuccore/AY672108):
  `repC`, `bla`, `xylR`, `tnp` mariner transposase, `tetD`, `tetB`). Both
  plasmids are cured at 45 °C, so these proteins should be **absent** — they are
  included as a control: if any of them turns up, a plasmid was not lost.

### Why this is worth re-searching

A search engine can only report sequences that are in its FASTA. `ermB` is not
in the JE2 database, so "ErmB not found" so far only means "never looked for".
Finding ErmB in the 18 mutant samples and not in the 9 wild-type ones would be
a **positive** control for the transposon — better evidence than the absence of
Stk1 and Stp1.

Two honest expectations:

* **ErmB is the real test.** It has 16 tryptic peptides and has never been
  searchable. Note that `ermB` is induced by erythromycin, so if these cultures
  grew without it the protein may sit below the detection limit; a negative
  result would then be uninformative. Check the growth protocol first.
* **GFP will probably stay negative.** The pBursa GFP differs from the
  *A. victoria* GFP already in the database (contaminant `zz_Y-FGCZCont00258_`)
  at only three residues (S65G, S72A, D155G), so **10 of its 12** tryptic
  peptides were already searchable — and none were detected.

### How to use it

```bash
cat p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta \
    transposon_proteins.fasta > JE2_db1_plus_transposon.fasta
```

4,155 + 11 = 4,166 entries, so the effect on the FDR estimate is negligible.

## Appendix: the original ErmB-only proposal


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
| `transposon_proteins.fasta` | the 11 plasmid-encoded proteins (see above). The ErmB in pBursa is byte-identical to UniProt P20173, which this file supersedes. |

## How to use it

```bash
cat p41881_db1_SA_JE2_modNmoreTranslations_20260630.fasta \
    reference/transposon_proteins.fasta > JE2_plus_transposon.fasta
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
