library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(parallel)

library(GenomicRanges)
library(plyranges)


# i/o functions -------------------------------------------------------------------------------

read_metadata <- function() {
  raw_info <- read_tsv("data/neo.impute.1000g.sampleInfo_clusterInfo.txt")

  info <-
    raw_info %>%
    filter(region %in% c("SouthernEurope", "WesternEurope", "WesternAsia", "NorthernEurope", "CentralEasternEurope"),
           country != "Greenland") %>%
    filter(groupAge != "Archaic") %>%
    mutate(ageAverage = ifelse(groupAge == "Modern", 0, ageAverage)) %>%
    mutate(coverage = ifelse(groupAge == "Modern", Inf, coverage))

  info
}

read_tracts <- function(set = NULL) {
  metadata <- read_metadata()
  
  if(!is.null(set)){
    info <- filter(metadata, groupAge == set) 
  }

  raw_tracts <- read_tsv(here::here("data/Vindija33.19_raw_eurasian_wModern"))

  tracts <- raw_tracts %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    filter(ID %in% unique(info$sampleId)) %>%
    select(ID, chrom, start, end) %>%
    mutate(length = end - start, set = set)

  tracts
}

# archaic deserts -----------------------------------------------------------------------------

generate_windows <- function(gaps_gr, window_size, step_size) {
  autosomes_gr <- GenomeInfoDb::getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  seqinfo(autosomes_gr) <- seqinfo(tracts_gr)

  windows_grl <- slidingWindows(autosomes_gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
    windows_grl[[i]]$gap <- FALSE
    to_remove <- queryHits(findOverlaps(windows_grl[[i]], gaps_gr))
    if (length(to_remove) > 0)
      windows_grl[[i]][to_remove]$gap <- TRUE
  }

  unlist(windows_grl)
}

generate_windows_sim <- function(n_chr=1, len_chr=100e6, window_size, step_size) {
  
  seq_info <- Seqinfo(paste0("chr", 1:n_chr), seqlengths = len_chr, isCircular = FALSE, genome = "sim")
  gr <- GRanges(paste0("chr", 1:n_chr), IRanges(rep(1, n_chr), len_chr), seqinfo = seq_info)
  
  windows_grl <- slidingWindows(gr, width = window_size, step = step_size)
  for (i in seq_along(windows_grl)) {
    windows_grl[[i]]$midpoint <- (start(windows_grl[[i]]) + end(windows_grl[[i]])) / 2
  }
  
  unlist(windows_grl)
}

compute_ancestry <- function(tracts_gr, windows_gr, keep_gaps = FALSE) {
  autosomes <- getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()

  # first compute coverage...
  cov <- coverage(tracts_gr)
  # ... then convert that to proportions
  cov <- lapply(cov, function(x) x / length(unique(tracts_gr$sampleId)))

  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
    chrom_coverage <- cov[[chrom]]
    chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]

    # count overlaps between windows and tracts
    average_coverage_per_window <- sapply(
      seq_along(chrom_gr),
      function(i) {
        start_idx <- start(chrom_gr[i])
        end_idx <- end(chrom_gr[i])
        mean(chrom_coverage[start_idx:end_idx])
    })

    mcols(chrom_gr)$coverage <- average_coverage_per_window
    if (!keep_gaps)
      mcols(chrom_gr)$coverage[mcols(chrom_gr)$gap] <- NA

    chrom_gr
  }, mc.cores = detectCores()-2)

  ancestry_grl <- GRangesList(ancestry_list)


  unlist(ancestry_grl)
}

plot_desert_ancestry <- function(ancestry_gr, deserts_gr, chrom, full = FALSE) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, ancient, modern, gap, midpoint)

  if (!full) {
    ancestry_df <- ancestry_df %>%
      filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1))
  }

  p <- ancestry_df %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.05) +

        geom_line(aes(midpoint, ancient), color = "blue") +

        geom_line(aes(midpoint, modern), color = "orange") +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$ancient, na.rm = TRUE), color = "blue", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern, na.rm = TRUE), color = "orange", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("blue", "orange")) +
        guides(color = guide_legend("", override.aes = list(size = 5)),
               linetype = guide_legend("")) +
        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        coord_cartesian(ylim = c(0, 0.5)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom", text = element_text(size = 13)) +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }
    if (!full) {
      p <- p +
        geom_point(data = filter(ancestry_df, ancient > 0), aes(midpoint, ancient, color = "ancient individuals"), size = 0.8) +
        geom_point(data = filter(ancestry_df, modern > 0), aes(midpoint, modern, color = "present-day individuals"), size = 0.8)
    }
    p
}

plot_desert_correlation <- function(ancestry_gr, chrom) {
  ancestry_df <- as_tibble(ancestry_gr) %>% filter(seqnames == chrom)

  rho <- ancestry_gr %>% filter(within_desert, seqnames == chrom) %>% { cor(.$modern, .$ancient) }

  ggplot() +
    geom_smooth(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert),
                formula = y ~ x, color = "red", fill = "black", method = "lm", linetype = "dashed", linewidth = 0.8, alpha = 0.35) +

    geom_point(data = filter(ancestry_df, !within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = desert, shape = "outside desert"),
               color = "lightgray", alpha = 0.5) +
    geom_point(data = filter(ancestry_df, within_desert, ancient > 0, modern > 0), aes(modern, ancient, color = within_desert, shape = "within desert"), color = "black") +

    geom_abline(slope = 1, linetype = "dashed") +

    geom_vline(aes(color = "modern", xintercept = mean(ancestry_df$modern, na.rm = TRUE)), linetype = "dashed", color = "blue") +
    geom_hline(aes(color = "ancient", yintercept = mean(ancestry_df$ancient, na.rm = TRUE)), linetype = "dashed", color = "orange") +

    scale_x_log10(breaks = c(0.0001, mean(ancestry_df$modern, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    scale_y_log10(breaks = c(0.0001, mean(ancestry_df$ancient, na.rm = TRUE), 1), labels = scales::percent_format(accuracy = 0.01), limits = c(0.00001, 1)) +
    labs(x = "Neanderthal ancestry proportion\nin present-day Eurasians [log scale]",
         y = "Neanderthal ancestry proportion\nin ancient Eurasians [log scale]") +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "bottom", text = element_text(size = 13)) +
    guides(shape = guide_legend("window", override.aes = list(alpha = 1, size = 3)),
           linetype = "none") +
    scale_shape_manual(values = c(4, 20)) +
    ggtitle("", subtitle = paste0("Pearson correlation within desert = ",
                                  formatC(rho, format = "f", digits = 3)))
}

plot_desert_ancestry2 <- function(ancestry_gr, deserts_gr, chrom) {
  desert_df <- deserts_gr %>% filter(seqnames == chrom) %>% as_tibble

  ancestry_df <- as_tibble(ancestry_gr) %>%
    filter(seqnames == chrom) %>%
    select(chrom = seqnames, start, end, modern, chen, gap, midpoint)

  ancestry_df %>%
    filter(start >= (desert_df$start * 0.9) & end <= (desert_df$end * 1.1)) %>%
    {
      ggplot(data = .) +
        geom_rect(data = desert_df, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "red", alpha = 0.1) +

        geom_line(aes(midpoint, chen), color = "orange") +
        geom_point(data = filter(., chen > 0), aes(midpoint, chen, color = "Chen et al."), size = 0.8) +

        geom_line(aes(midpoint, modern), color = "blue") +
        geom_point(data = filter(., modern > 0), aes(midpoint, modern, color = "Alba et al."), size = 0.8) +

        geom_vline(data = desert_df, aes(xintercept = start, linetype = "desert boundary"), color = "red") +
        geom_vline(data = desert_df, aes(xintercept = end, linetype = "desert boundary"), color = "red") +

        geom_hline(yintercept = mean(.$chen), color = "orange", linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = mean(.$modern), color = "blue", linetype = "dashed", alpha = 0.5) +

        scale_color_manual(values = c("Chen et al." = "orange", "Alba et al." = "blue")) +
        guides(color = guide_legend("", override.aes = list(size = 5)), linetype = guide_legend("")) +

        labs(x = "genomic coordinate [bp]", y = "proportion of Neanderthal ancestry") +
        scale_x_continuous(labels = scales::comma) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
        # coord_cartesian(ylim = c(0, 0.1)) +
        scale_linetype_manual(values = "dashed") +
        theme_minimal() +
        theme(legend.position = "bottom") +
        ggtitle(paste("Archaic ancestry desert on chromosome", gsub("chr", "", .$chrom[1])))
    }

}

# LD-based admixture dating -------------------------------------------------------------------

# Define positions of archaic ancestry informative sites
generate_info_sites <- function(tracts_gr, interval) {
  sites_grl <- lapply(seqlevels(tracts_gr), function(chrom) {
    positions <- seq(from = 1, to = seqlengths(tracts_gr)[chrom], by = interval)

    gr <- GRanges(seqnames = chrom, ranges = IRanges(start = positions, end = positions))
    mcols(gr)$index <- seq_len(length(gr))

    gr
  }) %>% GRangesList()
  seqlevels(sites_grl) <- seqlevels(tracts_gr)
  seqlengths(sites_grl) <- seqlengths(tracts_gr)

  sites_grl
}

# Define list of pairs of sites at given distances
# (one element of the list for each distance bin)
collect_pairs <- function(sites_grl, distances, recmap = NULL, ncores = parallel::detectCores()/2) {

  if (is.null(recmap)) {
    chrom_lengths <- seqlengths(sites_grl)
  } else {
    chrom_lengths <- split(recmap, recmap$chrom) %>% sapply(function(recmap_chr) max(recmap_chr$posg))
  }

  chroms <- sapply(sites_grl, function(x) as.character(unique(seqnames(x))))

  chr_pairs <- lapply(chroms, function(chrom) {

    sites_gr <- sites_grl[chroms == chrom, ] %>% unlist
    if (is.null(recmap)) {
      site_pos <- start(sites_gr)
    } else {
      site_pos <- convert_genetic(recmap, as.data.table(sites_gr), "start", chrom = "seqnames")$start_gen
    }

    pairs <- parallel::mclapply(distances, function(distance) {

      pair1 <- c()
      pair2 <- c()

      # iterate through each site one by one...
      for (i in sites_gr$index) {
        index1 <- i
        # ... and find the index of the first site that is at a given distance
        index2 <- sites_gr[site_pos >= site_pos[i] + distance]$index[1]

        if (is.na(index2)) {
          if (chrom_lengths[chrom] < site_pos[i] + distance  + distance / 10)
            break
          else
            next
        }

        # otherwise record the indices of the pair of sites and proceed with searching
        # for the next pair
        pair1 <- c(pair1, index1)
        pair2 <- c(pair2, index2)
      }

      list(pair1 = pair1, pair2 = pair2)

    }, mc.cores = ncores)

    pairs

  })

  names(chr_pairs) <- chroms
  chr_pairs
}

# Compute covariances of allele states at pairs of sites
compute_tract_covariances <- function(tracts_gr, sites_grl, pairs) {
  lapply(seqlevels(sites_grl), function(chrom) {

    chrom_sites_gr <- sites_grl[seqlevels(sites_grl) == chrom, ] %>% unlist

    parallel::mclapply(unique(tracts_gr$name), function(name) {

      ind_tracts_gr <- tracts_gr %>% filter(name == !!name, seqnames == chrom)
      ind_sites_gr <- chrom_sites_gr

      # mark sites falling within an introgressed tract
      tract_overlaps <- queryHits(findOverlaps(ind_sites_gr, ind_tracts_gr))
      mcols(ind_sites_gr)$neand <- FALSE
      if (length(tract_overlaps) > 0)
        mcols(ind_sites_gr[tract_overlaps])$neand <- TRUE
      mcols(ind_sites_gr)$neand <- as.integer(mcols(ind_sites_gr)$neand)

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_sites_gr[pairs[[chrom]][[i]]$pair1]$neand
        sites2 <- ind_sites_gr[pairs[[chrom]][[i]]$pair2]$neand
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = unique(ind_tracts_gr$sample_age),
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}


# Compute covariances of allele states at pairs of sites
compute_match_covariances <- function(info_gt, pairs, metadata) {
  archaic_name <- "NEA_1"
  samples <- setdiff(colnames(info_gt), c("chrom", "pos", archaic_name))

  lapply(unique(info_gt$chrom), function(chrom) {

    chrom_info_gt <- info_gt[, .SD[chrom %in% ..chrom, ], .SDcols = !c("chrom", "pos")]

    parallel::mclapply(samples, function(name) {

      ind_matches <- chrom_info_gt[, .(match = .SD[, get(name) == .SD[, get(archaic_name)]])]

      covariances <- sapply(seq_along(distances), function(i) {
        sites1 <- ind_matches[pairs[[chrom]][[i]]$pair1]$match
        sites2 <- ind_matches[pairs[[chrom]][[i]]$pair2]$match
        cov(sites1, sites2)
      })

      tibble(
        chrom = chrom,
        name = name,
        sample_age = filter(metadata, name == gsub("_hap\\d", "", !!name))$sample_age,
        distance = distances,
        covariance = covariances
      )
    }, mc.cores = detectCores()) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}

fit_exponential <- function(cov_df, distance) {
  cov_df <- cov_df %>% mutate(name = gsub("_hap\\d", "", name))

  distance <- match.arg(distance, choices = c("physical", "genetic"))

  grid_df <- expand_grid(name = unique(cov_df$name))

  fit_df <- lapply(1:nrow(grid_df), function(i) {
    name <- grid_df[i, ]$name

    data_df <- filter(cov_df, name == !!name) %>%
      group_by(name, sample_age, distance) %>%
      summarise(covariance = mean(covariance), .groups = "keep")

    lambda <- tryCatch({
      nls_res <- nls(covariance ~ SSasymp(distance, Asym, R0, lrc), data = data_df)
      exp(coef(nls_res)["lrc"])
    }, error = function(e) NA
    )

    if (!is.na(lambda)) {
      df <- tibble(
        distance = data_df$distance,
        covariance = predict(nls_res, newdata = data_df[, "distance"])
      )
    } else
      df <- NULL

    tibble(
      name = name,
      sample_age = data_df$sample_age[1],
      lambda = lambda,
      t_gens_before = ifelse(distance == "physical", lambda / 1e-8, lambda * 100),
      t_admix = t_gens_before * gen_time + sample_age,
      fit = list(df)
    )
  }) %>%
    do.call(rbind, .) %>%
    unnest(fit)
}


# conversion of physical distances to genetic distances ---------------------------------------

read_recmap <- function(path) {
  if (!dir.exists(path))
    stop("Path '", path, "' does not exist", call. = FALSE)

  files <- list.files(path, "plink.chr\\d+.*.map", full.names = TRUE)
  lapply(files, function(f) readr::read_tsv(f, col_names = c("chrom", "_", "posg", "pos"),
                                            show_col_types = FALSE)) %>%
    do.call(rbind, .) %>%
    mutate(chrom = paste0("chr", chrom))
}

read_recmap2 <- function(path, cumulative = TRUE) {
  if (!dir.exists(path))
    stop("Path '", path, "' does not exist", call. = FALSE)

  files <- list.files(path, ".*_recombination_map_hg19_chr_\\d+.bed", full.names = TRUE)
  lapply(files, function(f) {
    df <- read_tsv(f, show_col_types = FALSE) %>% setNames(c("chrom", "start", "end", "rate"))
    if (cumulative) {
      df <- rbind(tibble(chrom = df[1, ]$chrom, start = 0, end = df[1, ]$start, rate = 0), df) %>%
        mutate(posg = cumsum((end - start) * rate * 100)) %>%
        select(chrom, posg, pos = end)
    } else
      df <- mutate(df, length = end - start) %>% select(chrom, start, end, length, rate)
    df
  }) %>% do.call(rbind, .)
}

convert_genetic <- function(recmap, df, cols, chrom = "chrom") {
  interpolators <- recmap %>%
    split(.$chrom) %>%
    lapply(function(chrom_map) approxfun(chrom_map$pos, chrom_map$posg, rule = 2))

  df %>%
    { split(., .[[chrom]]) } %>%
    lapply(function(chrom_df) {
      if (!nrow(chrom_df)) return(NULL)
      chrom <- as.integer(gsub("chr", "", chrom_df[[chrom]][1]))
      for (c in cols) {
        chrom_df[[paste0(c, "_gen")]] <- interpolators[[!!chrom]](chrom_df[[c]])
      }
      chrom_df
    }) %>%
    do.call(rbind, .)
}

# Functions to compute tract frequency ------------------------------------

# This function computes the frequency of Neanderthal tract across windows
# NOTE: it returns them as a vector, ideally to be added at the windows_gr object as an metadata column

compute_tract_freq <- function(tracts_gr, windows_gr, age_group = NULL) {
  if(!is.null(age_group)){
    tracts_gr <- tracts_gr[tracts_gr$age_group==age_group]
  }
  autosomes <- getChromInfoFromUCSC("hg19") %>%
    filter(grepl("chr\\d+$", chrom)) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = size) %>%
    makeGRangesFromDataFrame()
  
  # first compute coverage...
  cov <- coverage(tracts_gr)
  # ... then convert that to proportions
  cov <- lapply(cov, function(x) x / length(unique(tracts_gr$sampleId)))
  
  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
                              chrom_coverage <- cov[[chrom]]
                              chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                              
                              # count overlaps between windows and tracts
                              # check Views object
                              cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                              average_coverage_per_window <- viewMeans(cov_view)
                              average_coverage_per_window
                            }, mc.cores = 1)
  
  # ancestry_grl <- GRangesList(ancestry_list)
  unlist(ancestry_list)
}

# Taken from Signac
StringToGRanges <- function(regions, sep = c("-", "-"), ...) {
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  return(granges)
}


# Plot frequencies over time ----------------------------------------------

# Take in inout a GR object

plot_frequencies <- function(gr_obj){
  # Get frequencies as df
  meta_df <- as.matrix(mcols(gr_obj))
  rownames(meta_df) <- paste0(seqnames(gr_obj), ":",ranges(gr_obj))
  meta_df <- meta_df %>% as_tibble(rownames="bin")
  # Long format for plotting
  meta_long <- meta_df %>% pivot_longer(cols = colnames(meta_df)[-c(1,ncol(meta_df))]) %>% transform(name=factor(name, levels=rev(c("present-day", "(0,2e+03]", "(2e+03,5e+03]", "(5e+03,1e+04]", "(1e+04,1.2e+04]" ,"(1.2e+04,3e+04]", "(3e+04,Inf]"))))
  p <- ggplot()+geom_line(data = meta_long, aes(x = name, y = as.numeric(value), group=bin, color=gene))+theme_minimal()
  print(p)
}


# Topic modeling ----------------------------------------------------------

compute_bin_vec <- function(tracts_gr, windows_gr, sample_id = NULL) {
  if(is.null(sample_id)){
    stop("Forgot to specify sampleID.")
  }
  
  # Get sample
  sample_gr <- tracts_gr %>% filter(sampleId==sample_id) 
  seqinfo(sample_gr) <- seqinfo(windows_gr)
  # first compute coverage...
  cov <- coverage(sample_gr)
  
  ancestry_list <- mclapply(as.character(unique(seqnames(windows_gr))),
                            function(chrom) {
                              chrom_coverage <- cov[[chrom]]
                              chrom_gr <- windows_gr[seqnames(windows_gr) == chrom]
                              
                              # count overlaps between windows and tracts
                              # check Views object
                              cov_view <- Views(chrom_coverage, ranges(chrom_gr))
                              average_coverage_per_window <- viewMeans(cov_view)
                              average_coverage_per_window
                            }, mc.cores = 1)
  
  # ancestry_grl <- GRangesList(ancestry_list)
  unlist(ancestry_list)
}


# Tract frame visualization -----------------------------------------------

plot_st_tract <- function(metadata, id_tract = NULL, dir_name = NULL) {
  if(is.null(id_tract)) stop("Provide a tract id.")
  if(is.null(dir_name)) stop("Provide a name for the directory.")
  if(!dir.exists(dir_name)) dir.create(dir_name)
  pdf(paste0(dir_name, "/st_", id_tract, ".pdf"))
  for (age in seq(27000, 1000, by=-1000)) {
    p <- metadata %>% filter(ageAverage > age) %>%
      filter(!is.na(latitude) & !is.na(longitude)) %>%
      st_as_sf(coords = c("longitude", "latitude")) %>%
      st_set_crs(4326) %>% st_jitter(amount = 1.4) %>%
      ggplot() +
      geom_sf(data = western_eurasia) +
      geom_sf(aes(
        alpha = ageAverage,
        color = as.factor(get(id_tract)),
        shape = as.factor(get(id_tract))
      )) +
      coord_sf(crs = 3035) + scale_alpha_binned(n.breaks = 15) + ggtitle(paste0(">", age))+ scale_color_manual(values=c("coral", "blue"))
    print(p)
  }
  dev.off()
}

# raster_map --------------------------------------------------------------

# stack
# https://gis.stackexchange.com/questions/375345/dividing-polygon-into-parts-which-have-equal-area-using-r
split_poly <- function(sf_poly, n_areas) {
  # Create random points
  set.seed(1234)
  points_rnd <- st_sample(sf_poly, size = 10000)
  # k-means clustering
  points <- do.call(rbind, st_geometry(points_rnd)) %>%
    as_tibble() %>% setNames(c("lon","lat"))
  k_means <- kmeans(points, centers = n_areas)
  # Create voronoi polygons
  voronoi_polys <- dismo::voronoi(k_means$centers, ext = sf_poly)
  # Clip to sf_poly
  crs(voronoi_polys) <- crs(sf_poly)
  voronoi_sf <- st_as_sf(voronoi_polys)
  equal_areas <- st_intersection(voronoi_sf, sf_poly)
  equal_areas$area <- st_area(equal_areas)
  return(equal_areas)
}

eucDist <- function(v1, v2) sqrt(rowSums(sweep(v1, 2, v2)^2))


# krigingST ---------------------------------------------------------------


CreateSpatioTemporalGrid <- function(data.UTM, sp.grid.UTM, numtimepoints, oldesttime = -Inf, youngesttime = 0){
  
  # Years as seconds
  dataTM <- as.POSIXlt(-data.UTM$TIME, origin="1970-01-01")
  oldestdata <- max( oldesttime, min(-data.UTM$TIME) )
  oldestdataTM <- as.POSIXlt(oldestdata, origin="1970-01-01")
  youngestdata <- min( youngesttime, max(-data.UTM$TIME) )
  youngestdataTM <- as.POSIXlt(youngestdata, origin="1970-01-01")
  
  # Create temporal grid
  if(!is.na(numtimepoints)){
    rawtimegrid <- round(seq(oldestdata,youngestdata,length.out=numtimepoints))
    tm.grid <- seq(as.POSIXct(oldestdataTM),as.POSIXct(youngestdataTM),length.out=numtimepoints)
  } else {
    rawtimegrid <- sort(unique(-data.UTM$TIME))
    tm.grid <- as.POSIXct(rawtimegrid, origin="1970-01-01")
  }
  
  # Create spatiotemporal grid
  grid.ST <- STF(sp.grid.UTM,tm.grid)
  
  return( list( dataTM, grid.ST, rawtimegrid, tm.grid))
  
}

BoundKriging <- function(KriggedData,min,max){
  leftbound <- which(KriggedData < min)
  rightbound <- which(KriggedData > max)
  KriggedData[leftbound,1] <- min
  KriggedData[rightbound,1] <- max
  return(KriggedData)
}

ComputeVariogram <- function(data.UTM, vartokrig, dataSP, dataTM, maxtlags=50, ncores = 3){
  dataDF <- data.frame(LABEL=data.UTM[[vartokrig]]) 
  timeDF <- STIDF(dataSP, dataTM, data=dataDF) 
  
  var <- variogramST(LABEL~1,data=timeDF,tunit="mins",assumeRegular=F,na.omit=T, tlags = 1:maxtlags, cores = ncores)
  # saveRDS(var, "variogramST.rds")
  
  anistart <- 10
  aniend <- 500
  anistep <- 10
  pars.l <- c(sill.s = 0, range.s = 0, nugget.s = 0,sill.t = 0, range.t = 0, nugget.t = 0,sill.st = 0, range.st = 0, nugget.st = 0, anis = 0)
  #pars.u <- c(sill.s = 1e7, range.s = 1e7, nugget.s = 1e7,sill.t = 1e7, range.t = 1e7, nugget.t = 1e7,sill.st = 1e7, range.st = 1e7, nugget.st = 1e7,anis = 1e7)
  finalVgmMSE <- Inf
  finalVgm <- NULL
  for( anisotropy in seq(anistart,aniend,anistep)){
    try( {
      metric <- vgmST("metric", joint = vgm(psill=1,"Exp", range=5e3, nugget=1e1), stAni=anisotropy)
      metric_Vgm <- fit.StVariogram(var, metric, method="L-BFGS-B",lower=pars.l,tunit="mins")
      mse <- attr(metric_Vgm,"MSE")
      #print(paste("Anisotropy: ",anisotropy,"; MSE: ",mse,sep=""))
      if(mse < finalVgmMSE){
        finalVgmMSE <- mse
        finalVgm <- metric_Vgm
      }
    }, silent = TRUE)
  }
  return(list(timeDF, var, finalVgm, finalVgmMSE))
}

