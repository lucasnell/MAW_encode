# Tidying background data for `human.R` and saving to `human.RData`


library('dplyr')
library('readr')
library('tidyr')
library('parallel')


options(stringsAsFactors = FALSE)



# =================
# PATHS TO BED FILES
# =================

# Path to folder containing BED files.
BEDparent <- './bed_files/human'

# We use broadPeak or gappedPeak peak BED files, for histone marks with diffuse or 
#   compact enrichment patterns, respectively.
# I found this information on Anshul Kundaje's website
#   - Author site: https://sites.google.com/site/anshulkundaje/projects/encodehistonemods
#   - DOI for paper he refers to: 10.1073/pnas.1318948111
#   - I saved author site on Evernote in case it changes:
#       "(2014) ENCODE: Histone modification ChIP-seq uniform peak calls - Anshul Kundaje"
bedFiles <- c(
    paste0(BEDparent, '/E097-', c('H3K9me3', 'H3K27me3', 'H3K36me3'), '.', 'broad', 
           'Peak.gz'),
    paste0(BEDparent, '/E097-', c('H3K4me1', 'H3K4me3', 'H3K27ac'), '.', 'gapped', 
           'Peak.gz')
    )
rm(BEDparent)

# =================
# SIZES AND BOUNDARIES
# =================

chrXSize <- 155270560

# PAR and XTR boundaries
boundsDF <- data.frame(PAR1 =  c(60001, 2699520),
                       PAR2 = c(154931044, 155260560),
                       # I believe below is 0-indexed, bc that's default in Perl, and
                       # I got these points from Cotter et al. (2016)'s Perl scripts
                       XTR = c(88193855, 93193855) + c(1, 0))

# Looked at chrX from hg19 on UCSC (`+1` is to convert to 1-based indexing)
firstStart <- 60000 + 1

# Adding a 3rd row with region lengths
boundsDF <- rbind(boundsDF, 
                  sapply(colnames(boundsDF), 
                         function(x){boundsDF[2,x] - boundsDF[1,x] + 1}))



# save(bedFiles, firstStart, boundsDF, chrXSize, numCores,
#      file = 'human.RData', compress = TRUE)








# Not using gene densities (yet) so keeping this commented

# =================
# GENE LOCATIONS
# =================

# # Data frame from GTF file (for gene densities)
# geneLocs <- read_tsv('/Volumes/MW_18TB/Lucas_Nell/lan/encode/human/',
#                                    'Homo_sapiens.GRCh37.75.gtf.gz'),
#                      skip = 5,
#                      col_names = c('chrom', 'source', 'feature', 'start', 'end',
#                                    'score', 'strand', 'frame', 'attribute'),
#                      col_types = paste0('cccii', 'cccc', collapse = ''))
# 
# 
# 
# # Parses attributes string to a matrix where [,1] is the attribute type (e.g., gene_id),
# #   and [,2] is the value of that attribute type
# # Then returns value of matrix [,2] in row where `focalAttribute` == value in [,1]
# parseAttributes <- function(attrDF, focalAttr, cores, attrColName = 'attribute'){
#     # Inner function to get attribute(s) from a single row
#     # stringr and `%>%` not used to improve speed
#     oneLine <- function(attributeString, focalAttribute){
#         xmat <- matrix(
#             strsplit(
#                 gsub(pattern = '[";]', replacement = '', 
#                      x = strsplit(
#                          attributeString, split = ', ')[[1]]), split = ' ')[[1]],
#             ncol = 2, byrow = TRUE)
#         if (length(focalAttribute) > 1){
#             outString <- paste0(xmat[xmat[,1] %in% focalAttribute,2], collapse = ';')
#         } else {
#             outString <- xmat[xmat[,1] == focalAttribute,2]
#         }
#         return(outString)
#     }
#     
#     if (cores > 1){
#         outVector <- mclapply(attrDF[[attrColName]], oneLine, 
#                               focalAttribute = focalAttr, 
#                               mc.cores = cores) %>% unlist
#     } else {
#         outVector <- sapply(attrDF[[attrColName]], oneLine, 
#                             focalAttribute = focalAttr, 
#                             USE.NAMES = FALSE)
#     }
#     
#     return(outVector)
# }
# 
# 
# 
# parsedAttr <- parseAttributes(geneLocs, c('gene_id', 'gene_name', 'gene_biotype'),
#                               cores = numCores)
# #    user  system elapsed
# # 166.483  12.926 203.847
# 
# 
# geneLocs <- geneLocs %>% 
#     mutate(attrInfo = parsedAttr) %>% 
#     separate(attrInfo, c('id', 'name', 'biotype'), sep = ';')
# 
# # ~~~~~~
# #  Maybe eventually do something like below
# # ~~~~~~
# # geneLocs %>%
# #     filter(feature %in% c("gene", "transcript", "exon", 
# #                           "CDS", "start_codon", "stop_codon"),
# #            !grepl('pseudogene', biotype)) %>% 
# #     select(chrom, start, end, id, name, biotype) %>%
# #     group_by(id) %>%
# #     summarize(chrom = min(chrom), start = min(start), end = max(end),
# #               name = min(name), biotype = min(biotype))
# # 
# # #  user  system elapsed
# # # 7.769   0.240   8.140



# save(geneLocs, bedFiles, firstStart, boundsDF, chrXSize, numCores,
#      file = 'human.RData', compress = TRUE)

  
  