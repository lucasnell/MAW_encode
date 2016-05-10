library('dplyr')
library('ggplot2')
library('stringr')
library('readr')
library('parallel')



# ==================================
# ==================================

# BACKGROUND INFO

# ==================================
# ==================================

# Objects created
#   tbl_df: geneLocs
#   data.frame: chrSizes, tissueRegions
#   numeric: firstStart, PAR_boundary
# See `mus_tidying.R` for how these were derived
load('mus.RData')
# Changing to 1-based indexing...
firstStart <- firstStart + 1

# Number of cores to use for parallel processes
numCores <- 3

options(stringsAsFactors = FALSE)

PAR_length <- chrSizes[1, 'chrX'] - PAR_boundary + 1




# ==================================
# ==================================

# READING BED FILES

# ==================================
# ==================================

readBED <- function(bedFile){
    
    # Derive sample attributes based on file name: target, tissue, sex, ID
    attrDF <- str_split(bedFile, '/') %>% unlist %>% tail(., 1) %>%
        str_split('_') %>% unlist %>% t %>%
        as.data.frame(., stringsAsFactors = F) %>%
        rename(target = V1, tissue = V2, sex = V3, id = V4) %>%
        mutate(id = str_replace(id, '.bed.gz', ''))
    
    bedDF <- read_tsv(gzfile(bedFile),
                      col_names = c('chrom', 'start', 'end', 'name', 'score',
                                    'strand', 'signal', 'p', 'q')) %>%
        mutate(target = attrDF[1,'target'], tissue = attrDF[1,'tissue'],
               sex = attrDF[1,'sex'], id = attrDF[1,'id'],
               # Converting to 1-based, inclusive indexing
               start = start + 1) %>%
        select(target, tissue, sex, id, chrom, start, end, signal) %>%
        as.tbl
    return(bedDF)
}

# Elapsed time typically ~5-10 seconds
allSampsDF <- lapply(bedFiles, readBED) %>% bind_rows

# No longer need these...
rm(bedFiles, readBED)





# ==================================
# ==================================

# PERMUTATIONS

# ==================================
# ==================================

# =================
# Permutations via moving the chrX peak locations
# =================

permuteXpeaks <- function(df, ID, B = 1000, returnVector = FALSE){
    
    X_DF <- df %>% filter(id == ID, chrom == 'chrX') %>%
        mutate(len = end - start + 1)
    X_n <- nrow(X_DF)
    
    # What proportion of peaks overlap PAR?
    PAR_peaks <- X_DF %>% filter(end >= PAR_boundary) %>% nrow
    obs <- PAR_peaks / X_n
    
    possibleEnds <- c(firstStart + min(X_DF$len), 
                      chrSizes[1, 'chrX'])
    
    onePerm <- function(possRange, sampleSize, cutOff){
        ends <- ceiling(runif(n = sampleSize, min = possRange[1], max = possRange[2]))
        n_obs <- length(ends[ends >= cutOff])
        return(n_obs / sampleSize)
    }
    
    permOvers <- sapply(seq(B), function(i) onePerm(possibleEnds, X_n, PAR_boundary))
    
    # Looking for enrichment, so using `>=`
    pval <- length(permOvers[permOvers >= obs]) / B
    
    if (returnVector){
        permDF <- data.frame(perms = permOvers)
        colnames(permDF) <- ID
        return(list(df_p = data.frame(id = ID, p = pval), df_perm = permDF))
    }
    
    return(data.frame(id = ID, p = pval))
    
}

# 
# RNGkind("L'Ecuyer-CMRG")
# set.seed(569)
# # Takes ~20-30 seconds
# perm_peak_locs <- lapply(allSampsDF$id %>% unique,
#                     function(i){ permuteXpeaks(allSampsDF, i,
#                                            returnVector = TRUE) })
# # save(perm_peak_locs, file = 'permute_X_peaks.RData', compress = TRUE)
# load('permute_X_peaks.RData')
# PAR_permSumm <- simplify2array(perm_peak_locs)[1,] %>% bind_rows
# PAR_perms <- simplify2array(perm_peak_locs)[2,] %>% bind_cols
# sigSamps <- PAR_permSumm %>% filter(p <= 0.05) %>% as.data.frame
# 
# summDF <- allSampsDF %>% 
#     filter(id %in% sigSamps$id) %>% 
#     group_by(id) %>%
#     summarize(target = target[1],
#               tissue = tissueRegions[1,tissue[1]],
#               sex = sex[1]) %>%
#     arrange(target, tissue, sex) %>%
#     select(target, tissue, sex, id)
# 
# summDF %>% rowwise %>% 
#     mutate(p = sigSamps$p[sigSamps$id == id]) %>% ungroup
# 
# # allSampsDF %>% filter(tissue == 'test') %>% select(target) %>% unique









# =================
# Permutations via moving the PAR location
# =================

permutePAR <- function(df, ID, B = 1e3, justX = FALSE, returnVector = FALSE){
    
    onePerm <- function(i){
        chrom <- sample(colNames, 1, prob = chromWeights)
        # Using `runif` bc it doesn't appear to create the sequence object first, so
        # saves mucho time. Incorporating `ceiling` makes it equivalent to `sample`.
        start_i <- ceiling(runif(1, firstStart, chrSizes[1, chrom] - PAR_length + 1))
        end_i <- start_i + PAR_length - 1
        peaks <- nrow(id_DF[id_DF$chrom == chrom & id_DF$start >= start_i & 
                                id_DF$end <= end_i, ])
        return(peaks)
    }
    
    id_DF <- df %>% filter(id == ID) %>% select(chrom, start, end)
    if (justX){
        colNames <- c('chrX')
    } else {
        colNames <- c(paste0('chr', seq(19)), 'chrX')
    }
    
    id_DF <- id_DF %>% filter(chrom %in% colNames)
    
    # How many peaks are fully inside PAR?
    PAR_peaks <- id_DF %>% filter(chrom == 'chrX', start >= PAR_boundary) %>% nrow
    
    # Weighting chromosome-name sampling by amount of sequence per chromosome
    # (not considering 3 Mb "N"s at beginning)
    chromWeights <- (chrSizes[1, colNames] - firstStart + 1) / 
        sum(chrSizes[1, colNames] - firstStart + 1)
    
    # If `PAR_peaks==0`, then every perm. will be `>=PAR_peaks`, so to save time:
    if (PAR_peaks == 0){
        permIn <- rep(1, B)
    } else {
        permIn <- sapply(seq(B), onePerm)
    }
    
    # Looking for enrichment, so using `>=`
    pval <- length(permIn[permIn >= PAR_peaks]) / B
    
    if (returnVector){
        permDF <- data.frame(perms = permIn)
        colnames(permDF) <- ID
        return(list(df_p = data.frame(id = ID, p = pval), df_perm = permDF))
    }
    
    return(data.frame(id = ID, p = pval))
}

# RNGkind("L'Ecuyer-CMRG")
# set.seed(946)
# perm_PAR_loc_X <- mclapply(allSampsDF$id %>% unique,
#                            function(i){ permutePAR(allSampsDF, i, justX = TRUE,
#                                                    returnVector = TRUE, B = 1e4) },
#                            mc.cores = numCores)
# # user  system elapsed
# # 80.603   1.969  50.313

# set.seed(328)
# perm_PAR_loc <- mclapply(allSampsDF$id %>% unique,
#                          function(i){ permutePAR(allSampsDF, i, justX = FALSE,
#                                                    returnVector = TRUE) },
#                          mc.cores = numCores)
# #    user  system elapsed
# # 245.868  47.605 180.661
# save(perm_PAR_loc, perm_PAR_loc_X, file = 'permute_PAR_location.RData', compress = T)
load('permute_PAR_location.RData')

X_perm <- simplify2array(perm_PAR_loc_X)[1,] %>% bind_rows
X_permVector <- simplify2array(perm_PAR_loc_X)[2,] %>% bind_cols

perm <- simplify2array(perm_PAR_loc)[1,] %>% bind_rows
permVector <- simplify2array(perm_PAR_loc)[2,] %>% bind_cols


X_sigSamps <- X_perm %>% filter(p <= 0.05)

allSampsDF %>% filter()

X_summDF <- allSampsDF %>%
    filter(id %in% X_sigSamps$id) %>%
    group_by(id) %>%
    summarize(target = target[1],
              tissue = tissueRegions[1,tissue[1]],
              sex = sex[1]) %>% 
    rowwise %>% 
    mutate(p = X_sigSamps$p[X_sigSamps$id == id]) %>% 
    ungroup %>%
    arrange(target, tissue, sex) %>%
    select(target, tissue, sex, id, p)

X_summDF


allSampsDF$target %>% unique






# # =================
# # SUMMARIZED DATA FRAME
# # =================
# 
# 
# # Note: You need to re-do this part.
# 
# summBed <- function(bedFile, B = 1000, boundary = PAR_boundary){
#     
#     # Derive sample attributes based on file name: target, tissue, sex, ID
#     attrDF <- str_split(bedFile, '/') %>% unlist %>% tail(., 1) %>%
#         str_split('_') %>% unlist %>% t %>%
#         as.data.frame(., stringsAsFactors = F) %>%
#         rename(target = V1, tissue = V2, sex = V3, id = V4) %>%
#         mutate(id = str_replace(id, '.bed.gz', ''))
#     
#     bedDF <- readBED(bedFile)
#     
#     summDF <- attrDF %>%
#         mutate(pval = permutePAR(bedDF, B, boundary)) %>%
#         as.tbl
#     
#     return(summDF)
# }
# 
# allSamps <- lapply(bedFiles, function(bed) summBed(bed, B = 100)) %>% 
#     bind_rows
# 
# 
# 
# allSamps %>% filter(pval <= 0.1) %>% as.data.frame
# 


























# =================
# GENE DENSITY
# =================

# Do this after permutations, only in samples you see enrichment

# _____ For the PAR... _____
# Length in Mb
PAR_length <- (chrSizes[1, 'chrX'] - PAR_boundary + 1) / 1e6
# Number of overlapping genes
PAR_genes <- geneLocs %>% filter(chrom == 'chrX', end >= PAR_boundary) %>% nrow
# Overlapping genes per Mb
PAR_density <- PAR_genes / PAR_length

# Function to get a region's gene density (genes/Mb)
getDense <- function(focalChrom, limits = NULL){
    if (is.null(limits)){
        len <- (chrSizes[1,focalChrom] - firstStart + 1) / 1e6
        genes <- geneLocs %>% filter(chrom == focalChrom, start >= firstStart) %>% nrow
    } else {
        if (length(limits) != 2){
            stop('`limits` should have length == 2.')
        }
        len <- (max(limits) - min(limits) + 1) / 1e6
        genes <- geneLocs %>% 
            filter(chrom == focalChrom, start <= max(limits), end >= min(limits)) %>% 
            nrow
    }
    geneDensity <- genes / len
    return(geneDensity)
}

Xwin <- seq(chrSizes[1,'chrX'] - (PAR_length*1e6) + 1, firstStart, -(PAR_length*1e6))
Xdens <- mclapply(Xwin, 
                  function(x){getDense('chrX', c(x, x + (PAR_length*1e6)))},
                  mc.cores = numCores) %>%
    unlist



matches <- which(Xdens == PAR_density)

# Gene density along chrX:
ggplot(aes(x = Xwin, y = Xdens), data = NULL) + 
    geom_hline(yintercept = PAR_density, linetype = 3) +
    geom_point(alpha = 0.4) + 
    ylab('Gene density') +
    xlab('Position on chromosome (Mb)') +
    geom_smooth(method = 'loess', span = 0.2, se = FALSE) +
    geom_point(data = NULL, inherit.aes = FALSE, color = 'red',
               aes(x = Xwin[matches], y = Xdens[matches]))



ggplot(aes(x = Xdens[Xdens > 0]), data = NULL) + 
    geom_histogram(bins = 20)


# 




































# =================
# PLOTTING
# =================

longBedDF <- bedDF %>% 
    rowwise %>%
    mutate(mid = median(start, end)) %>%
    ungroup %>%
    gather(position_type, position, -chrom, -signal)


longBedDF %>%
    filter(chrom == 'chrX', position_type == 'mid') %>%
    ggplot(aes(x = position, y = signal)) + 
    geom_vline(xintercept = PAR_boundary, linetype=2) +
    geom_vline(xintercept = chrSizes[1,'chrX'], linetype = 1) +
    geom_point(alpha = 0.5, shape = 21, color = NA, fill = 'black')


