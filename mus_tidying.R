# Tidying background data for `mus.R` and saving to `mus.RData`

library('dplyr')
library('stringr')
library('readr')

# Path to folder containing BED files. Change for your system.
BEDparent <- './bed'
# Paths to all BED files
bedFiles <- list.files(BEDparent) %>% 
  str_subset('.bed.gz') %>% 
  str_c(BEDparent, '/', .)

# Data frame from GTF file (for gene densities)
# Took a few minutes to process so I commented it to prevent running unnecessarily
# geneLocs <- read_tsv(gzfile('./bed/mm9_Ensembl.gtf.gz'), 
#                   col_names = c('chrom', 'source', 'feature', 'start', 'end',
#                                 'score', 'strand', 'frame', 'attribute')) %>%
#   filter(feature != 'exon') %>%
#   rowwise %>%
#   mutate(gene = (attribute %>% str_split('[;"]'))[[1]][2]) %>%
#   ungroup %>%
#   group_by(gene) %>% 
#   summarize(chrom = min(chrom), start = min(start), end = max(end)) %>%
#   arrange(chrom, start, gene) %>%
#   select(chrom, gene, start, end)


# Short and long versions of all tissue regions
tissueRegions <- 
  "olfa,olfactory bulb
bone,bone marrow
thym,thymus
test,testis
smal,small intestine
cere,cerebellum
brow,brown adipose tissue
sple,spleen
kidn,kidney
hear,heart
lung,lung
cort,cortical plate
plac,placenta
live,liver" %>% 
  read_csv(., col_names = c('short', 'long'))

tissueRegions <- (tissueRegions %>% 
                    t %>% as.data.frame(., stringsAsFactors = F) %>% 
                    setNames(tissueRegions$short))[2,]
row.names(tissueRegions) <- NULL



# Chromosome sizes
chrSizes <- 
  "chr1,197195432
chr2,181748087
chr3,159599783
chr4,155630120
chr5,152537259
chr6,149517037
chr7,152524553
chr8,131738871
chr9,124076172
chr10,129993255
chr11,121843856
chr12,121257530
chr13,120284312
chr14,125194864
chr15,103494974
chr16,98319150
chr17,95272651
chr18,90772031
chr19,61342430
chrX,166650296
chrY,15902555
chrM,16299" %>% 
  read_csv(., col_names = c('chrom', 'len'))

chrSizes <- (chrSizes %>% 
               t %>% as.data.frame(., stringsAsFactors = F) %>% 
               setNames(chrSizes$chrom))[2,] %>% 
  mutate_each(funs(as.numeric))
row.names(chrSizes) <- NULL

# PAR boundary location
PAR_boundary <- 166410000


# First possible beginning (3000000) determined from reading chrXb9 and chr19b9 
#   fasta files in python
firstStart <- 3000000



# save(geneLocs, chrSizes, tissueRegions, bedFiles, firstStart, PAR_boundary,
# file = 'mus.RData', compress = TRUE)
