source('preamble.R')


# ==================================
# ==================================

# BACKGROUND INFO

# ==================================
# ==================================

# __ Objects created __
#   (tbl_df)     --> geneLocs
#   (data.frame) --> chrSizes, tissueRegions
#   (numeric)    --> firstStart, PAR_boundary
#   (character)  --> bedFiles
# See `mus_tidying.R` for how these were derived
load('human.RData')
# Changing to 1-based indexing...
firstStart <- firstStart + 1

PAR_length <- chrSizes[1, 'chrX'] - PAR_boundary + 1

# Looked at chrX from hg19 on UCSC
firstStart <- 60000