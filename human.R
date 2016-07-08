source('preamble.R')



# ==================================
# ==================================

# BACKGROUND INFO

# ==================================
# ==================================

# __ Objects created __
#   (data.frame) --> boundsDF
#   (numeric)    --> firstStart, chrXSize, numCores
#   (character)  --> bedFiles
# See `human_tidying.R` for how these were derived
load('human.RData')






# ==================================
# ==================================

# READING BED FILES

# ==================================
# ==================================


readBED <- function(bedFile){
    
    # Derive antibody target based on filename
    antiTarget <- str_split(bedFile, '[-/]') %>% unlist %>% tail(., 1) %>%
        str_split('[.]') %>% unlist %>% head(., 1)
    
    # Whether gappedPeak BED format is used
    gappedPeak <- str_detect(bedFile, 'gapped')
    
    if (gappedPeak){
        colNames <- c('chrom', 'start', 'end', 'name', 'score', 'strand', 
                      'thickStart', 'thickEnd', 'itemRgb', # <-- useless for gappedPeak
                      'blockCount', 'blockSizes', 'blockStarts', 
                      'signal', 'p', 'q')
        colTypes <- paste0('ciiccc', 
                           'ccc',
                           'icc',
                           'cnn',
                           collapse = '')
        outCols <- c('target', 'gapped', 'chrom', 'start', 'end', 'p', 
                     'blockCount', 'blockSizes', 'blockStarts')
    } else {
        colNames <- c('chrom', 'start', 'end', 'name', 'score', 'strand', 
                      'signal', 'p', 'q')
        colTypes <- paste0('ciiccc', 
                           'cnn',
                           collapse = '')
        outCols <- c('target', 'gapped', 'chrom', 'start', 'end', 'p')
    }
    
    bedDF <- read_tsv(gzfile(bedFile),
                      col_names = colNames, col_types = colTypes) %>%
        mutate(target = antiTarget,
               gapped = gappedPeak,
               start = start + 1 # <-- converting to 1-based indexing
        ) %>%
        select_(.dots = outCols) %>% 
        filter(p >= -log10(0.01)) %>% # <-- p <= 0.01  (was <= 0.10)
        arrange(chrom, start)
    
    if (gappedPeak){
        
        # For separating below
        maxBlocks <- bedDF$blockCount %>% na.omit %>% max
        
        bedDF <- bedDF %>% 
            separate(blockSizes, paste0('b',seq(maxBlocks),'_sizes'), sep = ',', 
                     fill = 'right', convert = TRUE) %>%
            separate(blockStarts, paste0('b',seq(maxBlocks),'_starts'), sep = ',', 
                     fill = 'right', convert = TRUE) %>% 
            gather_('bN_sizes', 'bSize', paste0('b',seq(maxBlocks),'_sizes'), 
                    na.rm = TRUE) %>%
            gather_('bN_starts', 'bStart', paste0('b',seq(maxBlocks),'_starts'), 
                    na.rm = TRUE) %>%
            mutate(
                bN_starts = bN_starts %>% str_replace_all('[_a-z]', '') %>% 
                    as.numeric,
                bN_sizes = bN_sizes %>% str_replace_all('[_a-z]', '') %>% 
                    as.numeric
                ) %>%
            filter(bN_sizes == bN_starts) %>%
            mutate(start = start + bStart, end = start + bSize - 1) %>%
            select_(.dots = outCols[1:6])
    }
    
    return(bedDF)
}


# Even using 3/4 cores, it takes ~40 seconds
allSampsDF <- mclapply(bedFiles, readBED, mc.cores = numCores) %>% 
    bind_rows

allSampsChrX <- allSampsDF %>% filter(chrom == 'chrX')

# No longer need these...
rm(bedFiles, readBED)



# ==================================
# ==================================

# PERMUTATIONS

# ==================================
# ==================================

# =================
# Permutations via moving the focal regions (PAR1/2 and XTR) boundaries
#   and calculating proportion of new regions covered in peaks
# =================

# Get proportion of sequence range covered by ChIP-seq peaks
# `filteredDF` must be filtered by target and chrom

getProp <- function(filteredDF, start, end){
    limitEnd <- max(as.numeric(end))
    limitStart <- min(as.numeric(start))
    limitLength <- limitEnd - limitStart + 1
    # Filter for peaks that overlap limits
    overlappingDF <- filter(filteredDF, end >= limitStart, start <= limitEnd)
    if (nrow(overlappingDF) == 0){
        return(0)
    }
    newStarts <- ifelse(overlappingDF$start > limitStart, overlappingDF$start, limitStart)
    newEnds <- ifelse(overlappingDF$end < limitEnd, overlappingDF$end, limitEnd)
    indivOverPeaks <- newEnds - newStarts + 1
    overPeakSum <- sum(indivOverPeaks)
    peakProp <- overPeakSum / limitLength
    return(peakProp)
}

# Vectorized version
getPropVec <- Vectorize(getProp, c('start', 'end'), USE.NAMES = FALSE)




# A single permutation of all 3 focal regions' boundaries, making sure no overlapping
# or extension beyond chromosome limits occurs
onePerm <- function(i, chromLims, focalRegBounds, strippedDF){
    
    # Using `runif` bc it doesn't appear to create the sequence object first, so
    # saves mucho time. Incorporating `ceiling` and `-1` on the start position
    # makes it equivalent to `sample`.
    permStarts <- ceiling(runif(3, min(chromLims) - 1, 
                                max(chromLims) - min(focalRegBounds[3,]) + 1))
    names(permStarts) <- colnames(focalRegBounds)
    permEnds <- permStarts + as.numeric(focalRegBounds[3,]) - 1
    
    # Making sure they don't go past end of chromosome or overlap
    # Either of these should be rare, which is why I just used a while loop
    pastEnd <- any(permEnds > max(chromLims))
    overlap <- any((tail(sort(permStarts), -1) - head(sort(permEnds), -1)) <= 0)
    while (pastEnd | overlap){
        permStarts <- ceiling(runif(3, min(chromLims) - 1, 
                                    max(chromLims) - min(focalRegBounds[3,]) + 1))
        names(permStarts) <- colnames(focalRegBounds)
        permEnds <- permStarts + as.numeric(focalRegBounds[3,]) - 1
        pastEnd <- any(permEnds > max(chromLims))
        overlap <- any((tail(sort(permStarts), -1) - head(sort(permEnds), -1)) <= 0)
    }
    
    posDF <- data.frame(permNum = i,
                        region = colnames(focalRegBounds), 
                        start = permStarts, 
                        end = permEnds)
    
    rownames(posDF) <- NULL
    
    return(posDF)
}






# `filteredDF` should be filtered by target and chrom, but should still include
#   `"target"` column
# It's assumed that bounds in `focalRegBounds` are on the focal chrom
permuteProps <- function(filteredDF, focalRegBounds = boundsDF, 
                         chromLims = c(firstStart, chrXSize), 
                         B = 1e3, returnVector = FALSE){
    
    strippedDF <- filteredDF %>% select(start, end)
    
    focalRegBounds <- as.data.frame(focalRegBounds)
    
    # Observed proportions of each region covered in peak(s)
    obsProps <- sapply(colnames(focalRegBounds), 
                       function(x){
                           lims <- focalRegBounds[1:2,x]
                           prop <- getProp(strippedDF, lims[1], lims[2])
                           return(prop)
                       }, USE.NAMES = TRUE)
    
    permPosDF <- bind_rows(lapply(seq(B), onePerm, 
                                      chromLims = chromLims, 
                                      focalRegBounds = focalRegBounds, 
                                      strippedDF = strippedDF))

    # Final permutation df with peak proportion column 'prop'
    permDF <- permPosDF %>% 
        mutate(target = filteredDF$target[1],
               prop = getPropVec(start, end, filteredDF = strippedDF)) %>%
        select(target, permNum, region, prop)
    
    pDF <- permDF %>%
        group_by(region) %>%
        summarize(target = min(target), 
                  pval = {length(prop[prop >= obsProps[[min(region)]]]) / 
                          length(prop)}) %>%
        select(target, region, pval)
    
    if (returnVector){
        return(list(df_p = pDF, df_perm = permDF))
    } else {
        return(pDF)
    }
}




# RNGkind("L'Ecuyer-CMRG")
# set.seed(823)
# permX <- mclapply(allSampsChrX$target %>% unique,
#                             function(i){
#                                 permuteProps(filter(allSampsChrX, target == i),
#                                              returnVector = TRUE, B = 1e4)
#                             },
#                            mc.cores = numCores)
# # user  system elapsed
# # 143.871   5.221 130.470
# 
# save(permX, file = './R_data/human/permuteX.RData', compress = TRUE)

load('./R_data/human/permuteX.RData')


permP <- simplify2array(permX)[1,] %>% bind_rows
permV <- simplify2array(permX)[2,] %>% bind_rows

# write_csv(permP, './R_data/human/perm_summary.csv')

