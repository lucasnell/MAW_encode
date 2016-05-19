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
load('mus.RData')
# Changing to 1-based indexing...
firstStart <- firstStart + 1

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
# Permutations via moving the PAR location
#   and calculating proportion of new PAR covered in peaks
# =================

# Get proportion of sequence range covered by ChIP-seq peaks
# `DF` must be filtered by id
getProp <- function(DF, focalChrom = 'chrX', limits = c(PAR_boundary,chrSizes[1,'chrX'])){
    limitEnd <- max(limits)
    limitStart <- min(limits)
    limitLength <- limitEnd - limitStart + 1
    filteredDF <- filter(DF, chrom == focalChrom, end >= limitStart, start <= limitEnd)
    if (nrow(filteredDF) == 0){
        return(0)
    }
    newStarts <- ifelse(filteredDF$start > limitStart, filteredDF$start, limitStart)
    newEnds <- ifelse(filteredDF$end < limitEnd, filteredDF$end, limitEnd)
    indivOverPeaks <- newEnds - newStarts + 1
    overPeakSum <- sum(indivOverPeaks)
    peakProp <- overPeakSum / limitLength
    return(peakProp)
}


permutePARprop <- function(ID, df = allSampsDF, B = 1e3, justX = FALSE, 
                           returnVector = FALSE){
    
    id_DF <- df %>% filter(id == ID) %>% select(chrom, start, end)
    if (justX){
        colNames <- c('chrX')
    } else {
        colNames <- c(paste0('chr', seq(19)), 'chrX')
    }
    
    id_DF <- id_DF %>% filter(chrom %in% colNames)
    
    # What proportion of PAR is peak?
    PAR_peakProp <- getProp(id_DF)
    
    # Weighting chromosome-name sampling by amount of sequence per chromosome
    # (not considering 3 Mb "N"s at beginning)
    chromWeights <- (chrSizes[1, colNames] - firstStart + 1) / 
        sum(chrSizes[1, colNames] - firstStart + 1)
    
    onePerm <- function(i){
        chrom <- sample(colNames, 1, prob = chromWeights)
        # Using `runif` bc it doesn't appear to create the sequence object first, so
        # saves mucho time. Incorporating `ceiling` makes it equivalent to `sample`.
        newPAR_start <- ceiling(runif(1, firstStart, chrSizes[1, chrom] - PAR_length + 1))
        newPAR_end <- newPAR_start + PAR_length - 1
        peakProp <- getProp(id_DF, chrom, c(newPAR_start, newPAR_end))
        return(peakProp)
    }
    
    # If `PAR_peakProp==0`, then every perm. will be `>=PAR_peakProp`, so to save time:
    if (PAR_peakProp == 0){
        permIn <- rep(1, B)
    } else {
        permIn <- sapply(seq(B), onePerm)
    }
    
    # Looking for enrichment, so using `>=`
    pval <- length(permIn[permIn >= PAR_peakProp]) / B
    
    if (returnVector){
        permDF <- data.frame(perms = permIn)
        colnames(permDF) <- ID
        return(list(df_p = data.frame(id = ID, p = pval), df_perm = permDF))
    }
    
    return(data.frame(id = ID, p = pval))
}

# ~~~~~~~~~~~~~~~~
#  Running permutations
# (commented below bc I saved output as RData file, since it takes ~3 min)
# ~~~~~~~~~~~~~~~~
# RNGkind("L'Ecuyer-CMRG")
# set.seed(333)
# perm_PAR_prop_X <- mclapply(allSampsDF$id %>% unique,
#                            function(i){ permutePARprop(i, justX = TRUE,
#                                                        returnVector = TRUE, B = 1e4) },
#                            mc.cores = numCores)
# #    user  system elapsed
# # 414.085  15.637 185.571
# set.seed(833)
# perm_PAR_prop <- mclapply(allSampsDF$id %>% unique,
#                          function(i){ permutePARprop(i, justX = FALSE,
#                                                      returnVector = TRUE) },
#                          mc.cores = numCores)
# #    user  system elapsed
# # 216.568  23.381 133.877
# # save(perm_PAR_prop_X, perm_PAR_prop, file = './R_data/mus/permute_PAR_proportion.RData',
# # compress = T)

load('./R_data/mus/permute_PAR_proportion.RData')

X_perms <- simplify2array(perm_PAR_prop_X)[1,] %>% bind_rows
X_permVecs <- simplify2array(perm_PAR_prop_X)[2,] %>% bind_cols

ALL_perms <- simplify2array(perm_PAR_prop)[1,] %>% bind_rows
ALL_permVecs <- simplify2array(perm_PAR_prop)[2,] %>% bind_cols





# ==================================
# ==================================

# SUMMARY DATA FRAME

# ==================================
# ==================================

summDF <- allSampsDF %>%
    group_by(id) %>%
    summarize(target = target[1],
              tissue = tissueRegions[1,tissue[1]],
              sex = sex[1],
              PAR_peaks = length(chrom[chrom == 'chrX' & start >= PAR_boundary]),
              X_peaks = length(chrom[chrom == 'chrX'])) %>%
    rowwise %>%
    mutate(X_p = X_perms$p[X_perms$id == id],
           ALL_p = ALL_perms$p[ALL_perms$id == id]) %>%
    ungroup %>%
    arrange(target, tissue, sex, id) %>%
    select(target, tissue, sex, id, PAR_peaks, X_peaks, ALL_p, X_p)

# write_csv(summDF, './R_data/mus/perm_summary.csv')


sigSumm <- summDF %>% filter(X_p <= 0.05)





# ==================================
# ==================================

# GENE DENSITY

# ==================================
# ==================================



# _____ For the PAR... _____
# Number of overlapping genes
PAR_genes <- geneLocs %>% filter(chrom == 'chrX', end >= PAR_boundary) %>% nrow
# Overlapping genes per Mb
PAR_density <- PAR_genes / (PAR_length / 1e6)

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

# Getting gene densities of non-overlapping, PAR-sized windows across all chromosomes
getDenseDF <- function(focalChrom){
    wins <- seq(chrSizes[1,focalChrom] - PAR_length + 1, firstStart, -PAR_length)
    dens <- mclapply(wins, 
                     function(x){ getDense(focalChrom, c(x, x + PAR_length - 1)) },
                     mc.cores = numCores) %>%
        unlist
    return(data.frame(chrom = focalChrom, start = wins, dens = dens) %>% as.tbl)
}



denseDF <- lapply(c(paste0('chr', seq(19)), 'chrX'), getDenseDF) %>% bind_rows
# user  system elapsed
# 72.205  17.711  32.224

parDensDF <- denseDF %>% 
    filter(dens == PAR_density) %>% 
    mutate(end = start + PAR_length - 1) %>%
    select(chrom, start, end) %>%
    mutate(par = ifelse(chrom == 'chrX' & start == PAR_boundary, TRUE, FALSE))


getPARdensProps <- function(ID){
    DF <- allSampsDF %>% filter(id == ID)
    
    props <- apply(parDensDF, 1, function(row){
        getProp(DF, focalChrom = row['chrom'],
                limits = as.numeric(c(row['start'], row['end'])))
    })
    outDF <- data.frame(x = props)
    colnames(outDF) <- ID
    return(outDF)
}

PARdensProps <- lapply(sigSumm$id, getPARdensProps) %>% bind_cols


densProps <- list(parDensDF, PARdensProps) %>% 
    bind_cols %>%
    gather(sample, prop, -chrom, -start, -end, -par) %>% 
    rowwise %>% 
    mutate(target = sigSumm$target[sigSumm$id == sample],
           tissue = sigSumm$tissue[sigSumm$id == sample],
           sex = sigSumm$sex[sigSumm$id == sample]) %>% 
    ungroup %>%
    mutate(sample = factor(sample),
           target = factor(target),
           tissue = factor(tissue),
           sex = factor(sex)) %>%
    select(target, tissue, sex, chrom, start, end, par, sample, prop)



lan_theme <- function(base_size = 10, base_family = 'Helvetica') {
    theme_minimal(base_size = base_size, base_family = base_family) %+replace%
        theme(
            strip.text = element_text(face = 'bold'),
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "gray50"),
            axis.ticks = element_line(color = "gray50"),
            axis.ticks.length = unit(2, 'points'),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_text(face = 'plain', lineheight = 0.6)
        )
}



# Bootstrapped CI of median peak proportion for non-PAR windows of the same length and 
# gene density as the PAR
bootNonPar <- function(props, B = 1e3, fun = median){
    bootVec <- replicate(B, sample(props, replace = TRUE) %>% fun)
    ci <- quantile(bootVec, probs = c(0.025, 0.975)) %>% as.numeric
    return(paste(c(fun(props), ci), collapse = ':'))
}

set.seed(673)
densBoots <- densProps %>% 
    filter(!par) %>% 
    group_by(target, tissue) %>%
    summarize(ci = bootNonPar(prop)) %>%
    separate(ci, c('prop', 'lo', 'hi'), sep = ':', convert = TRUE)

set.seed(433)
densMeanBoots <- densProps %>% 
    filter(!par) %>% 
    group_by(target, tissue) %>%
    summarize(ci = bootNonPar(prop, fun = mean)) %>%
    separate(ci, c('prop', 'lo', 'hi'), sep = ':', convert = TRUE)

densPlot <- densProps %>%
    filter(!par) %>%
    ggplot(aes(x = 1, y = log2(prop))) +
    lan_theme() +
    ylab(expression("Peak proportion: " ~ log[2]*"[ "*Mb[peak] ~ Mb[total]^-1*" ]")) +
    geom_point(alpha = 0.4, position = position_jitter(height = 0, width = 0.2), 
               shape = 21, color = NA, fill = 'black', size = 1) + 
    geom_hline(data = densProps %>% filter(par), aes(yintercept = log2(prop)),
               color = 'dodgerblue', size = 1) +
    geom_errorbar(data = densBoots, aes(x = 1.0, ymax = log2(hi), ymin = log2(lo)), 
                  width = 0.02, size = 0.75, color = 'black', alpha = 0.7) +
    geom_point(data = densBoots, aes(x = 1.0), alpha = 0.7, size = 3, 
               color = 'black', shape = 18) +
    facet_grid(. ~ target + tissue, labeller = label_both) +
    geom_text(data = densProps %>% filter(par, target == 'H3K27me3', tissue == 'kidney'),
              aes(x = 1.1, y = log2(prop)), hjust = 1, vjust = -0.5, family = 'mono',
              label = 'PAR', color = 'dodgerblue', fontface = 'bold') +
    geom_text(data = densProps %>% filter(par, target == 'H3K27me3', tissue == 'kidney'),
              aes(x = 1.0, y = log2(prop) + 1.0), hjust = 0.5, vjust = 0.5,
              label = 'non-PAR', color = 'black', alpha = 0.4, fontface = 'bold',
              family = 'mono')


densMeanPlot <- densProps %>%
    filter(!par) %>%
    ggplot(aes(x = 1, y = log2(prop))) +
    lan_theme() +
    ylab(expression("Peak proportion: " ~ log[2]*"[ "*Mb[peak] ~ Mb[total]^-1*" ]")) +
    geom_point(alpha = 0.4, position = position_jitter(height = 0, width = 0.2), 
               shape = 21, color = NA, fill = 'black', size = 1) + 
    geom_hline(data = densProps %>% filter(par), aes(yintercept = log2(prop)),
               color = 'dodgerblue', size = 1) +
    geom_errorbar(data = densMeanBoots, aes(x = 1, ymax = log2(hi), ymin = log2(lo)), 
                  width = 0.02, size = 0.75, color = 'black', alpha = 0.7) +
    geom_point(data = densMeanBoots, aes(x = 1), alpha = 0.7, size = 3, 
               color = 'black', shape = 18) +
    facet_grid(. ~ target + tissue, labeller = label_both) +
    geom_text(data = densProps %>% filter(par, target == 'H3K27me3', tissue == 'kidney'),
              aes(x = 1.1, y = log2(prop)), hjust = 1, vjust = -0.5, family = 'mono',
              label = 'PAR', color = 'dodgerblue', fontface = 'bold') +
    geom_text(data = densProps %>% filter(par, target == 'H3K27me3', tissue == 'kidney'),
              aes(x = 1.0, y = log2(prop) + 1.0), hjust = 0.5, vjust = 0.5,
              label = 'non-PAR', color = 'black', alpha = 0.6, fontface = 'bold',
              family = 'mono')


densPlot
densMeanPlot





# 
# # Gene density along chrX:
# ggplot(aes(x = Xwin, y = Xdens), data = NULL) + 
#     geom_hline(yintercept = PAR_density, linetype = 3) +
#     geom_point(alpha = 0.4) + 
#     ylab('Gene density') +
#     xlab('Position on chromosome (Mb)') +
#     geom_smooth(method = 'loess', span = 0.2, se = FALSE) +
#     geom_point(data = NULL, inherit.aes = FALSE, color = 'red',
#                aes(x = Xwin[matches], y = Xdens[matches]))
# 
# 
# 
# ggplot(aes(x = Xdens[Xdens > 0]), data = NULL) + 
#     geom_histogram(bins = 20)
# 

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


