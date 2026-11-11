# constants.R -- the experimental design and a few proteins we will meet often.

# Factor 1: genotype. "anc" = JE2 wild type (the "ancestor").
GENOTYPES <- c("anc", "stk", "stp")

# Factor 2: medium. We order them by "how much / how long PASN":
#   TSB     = control medium
#   SNpulse = grown in TSB, then 15 min PASN (immediate response)
#   SN      = grown ~15 h in 70 % TSB + 30 % cell-free, sterile-filtered spent
#             medium of P. aeruginosa PAO1 ("PASN"; secreted molecules only,
#             no living Pseudomonas; partly nutrient-depleted) -- long-term
#   (SNpulse = the same supernatant, but only for 15 min after growth in TSB)
# In the raw file names SNpulse is written "pulseSN"; we always use "SNpulse",
# because that is how the contrasts in the DEA tables are named.
MEDIA <- c("TSB", "SNpulse", "SN")

# The 9 groups (genotype_medium), in plotting order.
GROUPS <- as.vector(t(outer(GENOTYPES, MEDIA, paste, sep = "_")))

# Colours for genotypes (colour-blind friendly, Okabe-Ito palette).
GENOTYPE_COLOURS <- c(anc = "#0072B2", stk = "#E69F00", stp = "#009E73")

# Proteins that appear in the exercises (JE2 = GenBank CP020619.1 locus tags).
KEY_PROTEINS <- tibble::tribble(
  ~short_name, ~locus_tag,     ~protein_id,  ~role,
  "Stp1",      "B7H15_06350",  "ARG45702.1", "Ser/Thr phosphatase (mutated in 'stp')",
  "Stk1",      "B7H15_06355",  "ARG45703.1", "Ser/Thr kinase PknB (mutated in 'stk'), downstream of stp1",
  "TcyA",      "B7H15_13375",  "ARG46976.1", "L-cystine ABC transporter, substrate-binding protein",
  "TcyB",      "B7H15_13370",  "ARG46975.1", "L-cystine ABC transporter, permease",
  "TcyC",      "B7H15_13365",  "ARG46974.1", "L-cystine ABC transporter, ATPase",
  "TcyP-like", "B7H15_02105",  "ARG44923.1", "L-cystine symporter (alternative cystine uptake)",
  "AhpC",      "B7H15_02095",  "ARG44921.1", "alkyl hydroperoxide reductase (ROS scavenger)",
  "NorA",      "B7H15_03880",  "ARG45237.1", "multidrug efflux pump"
)

# Default thresholds used by the prolfquapp report (FDR < 0.1 and |log2FC| > 1).
FDR_CUTOFF <- 0.1
LFC_CUTOFF <- 1

# ---- The two transposon mutants --------------------------------------------
# Both strains come from the Nebraska Transposon Mutant Library (NTML), which
# is built in JE2 with the mariner-based transposon "bursa aurealis":
#   3,237 bp, carries ermB (erythromycin resistance) and gfp, delivered from the
#   temperature-sensitive plasmid pBursa, transposase supplied from pFA545.
# Insertion coordinates are from the NTML Tn insertion map (JE2 numbering),
# https://ntml.unmc.edu/tninsertionmap ; the same two mutants were used by
# Prust et al. 2021 (Mol Cell Proteomics), who call them NE217 and NE1919.
#
# `insertion_pos` is the position in the JE2 genome (CP020619.1); cds_start /
# cds_end come from the JE2 annotation. Notebook 07 recomputes `codon` from
# these numbers rather than trusting the table, which is the point of the
# exercise.
TRANSPOSON_INSERTIONS <- tibble::tribble(
  ~short_name, ~ntml_id,  ~locus_tag,    ~cds_start, ~cds_end,  ~insertion_pos, ~orientation,
  "Stp1",      "NE1919",  "B7H15_06350",    1232295,  1233038,         1232861, "F",
  "Stk1",      "NE217",   "B7H15_06355",    1233035,  1235029,         1233171, "F"
)

# Size of the inserted element, for the polar-effect discussion in notebook 07.
TRANSPOSON_LENGTH_BP <- 3237
