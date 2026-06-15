

# Required libraries for functions
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

subset_counts_for_lmm <- function(counts, 
                                   annotation, 
                                   subset.list){ 
  
  subset.counts <- counts
  subset.annotation <- annotation
  
  # Subset the object based on the given annotations
  for(column in names(subset.list)){ 
    
    subset.annotation <- subset.annotation %>% 
      filter(.[[column]] %in% subset.list[[column]])
    
    subset.IDs <- subset.annotation$sample_ID
    
    subset.columns <- c("gene", subset.IDs)
    
    subset.counts <- subset.counts %>% 
      dplyr::select(all_of(subset.columns))
    
    # Factor the columns with relevant annotations
    subset.annotation[[column]] <- factor(subset.annotation[[column]])
    
  }
  
  # Factor the slide column
  subset.annotation[["slide_name"]] <- factor(subset.annotation[["slide_name"]])
  
  # Create log2 counts
  subset.counts.log2 <-  subset.counts %>%
    mutate(across(where(is.numeric), log2))
  
  return(list("subset.counts" = subset.counts, 
              "subset.log.counts" = subset.counts.log2, 
              "subset.annotation" = subset.annotation))
  
}

subset_object_for_lmm <- function(object, 
                                  subset.group.1,
                                  subset.group.2){ 
  
  # Set up the object to subset
  subset.object <- object
  
  # Get the annotation
  annotation <- pData(object)
  
  
  # Get the AOIs for the first group
  subset.1.annotation <- annotation
  
  for(field in names(subset.group.1)){
    
    values <- subset.group.1[[field]]

    subset.1.annotation <- subset.1.annotation %>%
      filter(.data[[field]] %in% values)
    
  }
  
  # Final AOI list for the first group
  subset.1.AOIs <- rownames(subset.1.annotation)
  
  # Gather the AOIs for the second group
  subset.2.annotation <- annotation
  
  for(field in names(subset.group.2)){
    
    values <- subset.group.2[field]

    subset.2.annotation <- subset.2.annotation %>%
      filter(.data[[field]] %in% values)
    
  }
  
  # Final AOI list for the second group
  subset.2.AOIs <- rownames(subset.2.annotation)
  
  # Combine into a total AOI list
  subset.AOIs <- c(subset.1.AOIs, subset.2.AOIs)
  
  # Create the subset returns
  subset.object <- object[,subset.AOIs]
  
  subset.annotation <- pData(subset.object)
  subset.annotation$sampleID <- gsub(".dcc", "", rownames(subset.annotation))
  
  # Create log2 counts for normalized and raw
  assayDataElement(object = subset.object, elt = "log_q") <-
    assayDataApply(subset.object, 2, FUN = log, base = 2, elt = "q_norm")
  
  assayDataElement(object = subset.object, elt = "log_raw") <-
    assayDataApply(subset.object, 2, FUN = log, base = 2, elt = "exprs")
  
  return(list("subset.object" = subset.object, 
              "subset.annotation" = subset.annotation))
  
}


run_limma <- function(counts, 
                      annotation, 
                      contrast, 
                      contrast.levels){
  
  # Create the DGE object
  DGE.list <- DGEList(counts = counts, 
                      samples = annotation)
  
  # Create the LM model design
  design <- model.matrix(formula(paste0("~ 0 + ", contrast)), 
                         data = DGE.list$samples)
  
  # Gather slide correlation for using random effect
  corfit <- duplicateCorrelation(DGE.list$counts, 
                                 design, 
                                 block = DGE.list$samples$slide_name)
  
  # Fit model with blocking (adjust for slide variance)
  fit <- lmFit(DGE.list$counts, 
               design, 
               block = DGE.list$samples$slide_name, 
               correlation = corfit$consensus)
  
  # Set up the contrast
  contrast.level.ref <- paste0(contrast, contrast.levels[[1]])
  contrast.level.condition <- paste0(contrast, contrast.levels[[2]])
  
  
  contrast <- makeContrasts(paste0(contrast.level.condition, 
                                   " - ", 
                                   contrast.level.ref),
                            levels = colnames(coef(fit)))
  
  # Generate the estimate of the contrast
  contrast.estimate <- contrasts.fit(fit, contrast)
  
  # Run Empirical Bayes smoothing of standard errors 
  fit.eb <- eBayes(contrast.estimate, robust = TRUE)
  
  # Generate the results table
  results <- topTable(fit.eb, sort.by = "P", n=Inf)
  
  
  return(list("results" = results, 
              "fit" = fit.eb, 
              "design" = design))
}



run_lmm <- function(object, contrast, within.slide){
  
  if(within.slide == TRUE){
    
    # Run the linear model with random slope
    lmm.results <- mixedModelDE(object, 
                                elt = "log_q", 
                                modelFormula = formula(paste0("~ 1 + ", 
                                                              contrast, 
                                                              " + (1 + " , 
                                                              contrast, 
                                                              " | slide_name)")), 
                                groupVar = contrast, 
                                nCores = parallel::detectCores(), 
                                multiCore = TRUE)
  } else {
    
    lmm.results <- mixedModelDE(object, 
                                elt = "log_q", 
                                modelFormula = formula(paste0("~ 1 + ", 
                                                              contrast, 
                                                              " + (1 | slide_name)")), 
                                groupVar = contrast, 
                                nCores = parallel::detectCores(), 
                                multiCore = TRUE)
    
  }
  
  # Gather the results into an output table
  lmm.results.summary <- do.call(rbind, lmm.results["lsmeans", ])
  lmm.results.summary <- as.data.frame(lmm.results.summary)
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  lmm.results.summary$gene <- 
    unlist(lapply(colnames(lmm.results),
                  rep, nrow(lmm.results["lsmeans", ][[1]])))
  
  # Run multiple test correction
  lmm.results.summary$FDR <- p.adjust(lmm.results.summary$`Pr(>|t|)`,
                                      method = "fdr")
  
  
  # Rename columns
  lmm.results.summary$pval <- lmm.results.summary[["Pr(>|t|)"]]
  lmm.results.summary$adj_pval <- lmm.results.summary$FDR
  lmm.results.summary$logfc <- lmm.results.summary$Estimate
  
  # Format final summary data frame
  lmm.results.summary <- lmm.results.summary[, c("gene", "logfc", 
                                                 "pval", "adj_pval")]
  
  return(list("results" = lmm.results.summary, "lm.output" = lmm.results))
  
}

# Set up the MA plot table
make_MA <- function(contrast.field, 
                    condition.label, 
                    reference.label, 
                    results.df, 
                    log.counts, 
                    raw.log.counts, 
                    annotation){
  
  # Gather the sample IDs for condition and reference groups
  condition.samples <- rownames(annotation[annotation[[contrast.field]] == condition.label, ])
  reference.samples <- rownames(annotation[annotation[[contrast.field]] == reference.label, ])
  
  # Gather normalized and raw counts for both groups
  condition.counts <- as.data.frame(log.counts[, condition.samples])
  reference.counts <- as.data.frame(log.counts[, reference.samples])
  
  condition.raw.counts <- as.data.frame(raw.log.counts[, condition.samples])
  reference.raw.counts <- as.data.frame(raw.log.counts[, reference.samples])  
  
  # Get the mean log score for each gene for both 
  # normalized counts
  condition.row.order <- rownames(condition.counts)
  condition.counts <- as.data.frame(sapply(condition.counts, as.numeric))
  condition.counts$cond_mean <- rowMeans(condition.counts)
  condition.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.counts)
  reference.counts <- as.data.frame(sapply(reference.counts, as.numeric))
  reference.counts$ref_mean <- rowMeans(reference.counts)
  reference.counts$gene <- reference.row.order
  
  # raw counts
  condition.row.order <- rownames(condition.raw.counts)
  condition.raw.counts <- as.data.frame(sapply(condition.raw.counts, as.numeric))
  condition.raw.counts$cond_raw_mean <- rowMeans(condition.raw.counts)
  condition.raw.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.raw.counts)
  reference.raw.counts <- as.data.frame(sapply(reference.raw.counts, as.numeric))
  reference.raw.counts$ref_raw_mean <- rowMeans(reference.raw.counts)
  reference.raw.counts$gene <- reference.row.order
  
  
  # Create a new data frame of the gene and group means with M and A values
  normalized.counts <- merge(condition.counts, reference.counts, by = "gene") %>% 
    dplyr::select(gene, cond_mean, ref_mean) %>% 
    mutate(M.value = cond_mean - ref_mean) %>% 
    mutate(A.value = (cond_mean + ref_mean)/2)
  
  raw.counts <- merge(condition.raw.counts, reference.raw.counts, by = "gene") %>% 
    dplyr::select(gene, cond_raw_mean, ref_raw_mean) %>% 
    mutate(M.raw.value = cond_raw_mean - ref_raw_mean) %>% 
    mutate(A.raw.value = (cond_raw_mean + ref_raw_mean)/2)
  
  # Add the DE results and log counts together
  ma.plot.counts <- merge(normalized.counts, raw.counts, by = "gene")
  
  # Set the bounds for the y axix so that they are aligned
  min.y <- min(c(min(ma.plot.counts$M.value),min(ma.plot.counts$M.raw.value)))
  max.y <- max(c(max(ma.plot.counts$M.value),max(ma.plot.counts$M.raw.value)))
  
  ma.plot.norm <- ggplot(ma.plot.table, aes(x = A.value, y = M.value)) +
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=lm, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Post-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  ma.plot.raw <- ggplot(ma.plot.table, aes(x = A.raw.value, y = M.raw.value)) + 
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=lm, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Pre-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  combined.MA.plots <- arrangeGrob(ggplotGrob(ma.plot.raw), 
                                   ggplotGrob(ma.plot.norm), 
                                   nrow = 1, ncol = 2)
  
  return(combined.MA.plots)
  
  
  
  
}

run_GSEA <- function(){
  
  
  
}

make_heatmap <- function(normalized.log.counts.df = q3.norm.log.counts, 
                         de.results = NULL, 
                         top.degs = FALSE, 
                         top.variable = FALSE, 
                         logfc.column = NULL, 
                         logfc.cutoff = NULL, 
                         annotation.column, 
                         annotation.row = NULL, 
                         anno.colors, 
                         cluster.rows = TRUE, 
                         cluster.columns = TRUE, 
                         main.title, 
                         row.gaps = NULL, 
                         column.gaps = NULL, 
                         show.rownames = FALSE, 
                         show.colnames = FALSE, 
                         min.genes.to.display = 2, 
                         max.genes.to.display = 500, 
                         font.size.row = 4){
  
  
  
  if(top.degs == TRUE & top.variable == TRUE){ 
  
    stop("Set only one of top.degs or top.variable to TRUE, not both")  
    
  }
  
  #if (top.variable == TRUE){
    
  #}
  
  # Filter genes by top DEGs, if applicable
  if(top.degs == TRUE){ 
    
    # Arrange by adjusted p-value
    degs.df <- de.results %>% 
      filter(padj < 0.05) %>% 
      arrange(desc(padj))
    
    # Arrange by log FC
    degs.df <- degs.df %>% arrange(desc(logfc))
    
    if(!is.null(logfc.cutoff)){
      
      degs.df <- degs.df %>% 
        filter(.data[[logfc.column]] > logfc.cutoff | .data[[logfc.column]] < -(logfc.cutoff))
    }
    
    # Use DEGs meeting adj p-val and logfc cutoffs
    if(length(rownames(degs.df)) >= min.genes.to.display){
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [logfc (+-", logfc.cutoff, ")")
      
      main.title <- paste0(main.title, " adj p-val (<0.05)]")
      
    }    
    
    # Revert to adj p-val cutoff and no logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(padj < 0.05) %>% 
        arrange(desc(padj))
      
      print("Not enough DEGs with listed logFC cutoff, reverting to all DEGs with adj p-value < 0.05")
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [no logfc cutoff")
      
      main.title <- paste0(main.title, " adj p-val (<0.05)]")
      
    }
    
    # Revert to nonadj p-val cutoff and logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(pval < 0.05) %>% 
        arrange(desc(padj)) %>% 
        filter(.data[[logfc.column]] > logfc.cutoff | .data[[logfc.column]] < -(logfc.cutoff))
      
      print("Not enough DEGs with adj p-val < 0.05, reverting to all DEGs with p-value < 0.05")
      
      main.title <- paste0(main.title, " [logfc (+-", logfc.cutoff, ")")
      
      main.title <- paste0(main.title, " p-val (<0.05, no adj)]")
      
    }
    
    # Revert to only p-val if no DEGs with p-val cutoff and logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(pval < 0.05) %>% 
        arrange(desc(pval))
      
      print("Not enough DEGs with adj p-val < 0.05, reverting to NON-adjusted p-value < 0.05")
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [no logfc cutoff ")
      
      main.title <- paste0(main.title, " p-val (<0.05, no adj)]")
      
    }
    
    # If there are more then 500 DEGs, trim down to top 500
    if(length(rownames(degs.df)) > max.genes.to.display){
      degs.df <- degs.df %>% slice(1:max.genes.to.display)
    }
    
    # Grab the list of DEGs
    degs.list <- degs.df$gene
    
    # Subset the counts df for the DEGs and order based on the DEGs list
    counts <- normalized.log.counts.df[rownames(normalized.log.counts.df) %in% degs.list, ]
    counts <- counts[match(degs.list, rownames(counts)), ]
    
  } else {
    
    counts <- normalized.log.counts.df
    
  }
  
  # Arrange by annotations if no column clustering
  if(cluster.columns == FALSE){
    
    # First arrange the annotation by the annotation groups
    anno.col.names <- colnames(annotation.column)
    
    for(col in anno.col.names){
      
      #annotation.column[[col]] <- as.factor(annotation.column[[col]])
      
      annotation.column <- annotation.column %>% 
        arrange(.data[[col]])
      
    }
    
    # Next match the counts file to the row order of the annotation
    counts <- counts[, rownames(annotation.column), drop = FALSE]
    
  }
  
  heatmap.plot <- pheatmap(counts, 
                           main = main.title, 
                           show_rownames = show.rownames, 
                           scale = "row",   
                           show_colnames = show.colnames,
                           border_color = NA, 
                           cluster_rows = cluster.rows, 
                           cluster_cols = cluster.columns, 
                           clustering_method = "average", 
                           clustering_distance_rows = "correlation", 
                           clustering_distance_cols = "correlation", 
                           color = colorRampPalette(c("blue", "white", "red"))(120), 
                           annotation_row = annotation.row, 
                           annotation_col = annotation.column,  
                           annotation_colors = anno.colors, 
                           gaps_row = row.gaps, 
                           gaps_col = column.gaps, 
                           fontsize_row = font.size.row)
  
  
  
  return(heatmap.plot)
  
}

make_heatmap_gene_sets <- function(normalized.log.counts.df,
                         gene.sets,
                         annotation.column,
                         annotation.row = NULL,
                         anno.colors,
                         cluster.rows = FALSE,
                         cluster.columns = TRUE,
                         main.title,
                         row.gaps = NULL,
                         column.gaps = NULL,
                         show.rownames = FALSE,
                         show.colnames = FALSE,
                         min.genes.to.display = 2,
                         max.genes.to.display = 500,
                         font.size.row = 4,
                         keep.empty_sets = FALSE,
                         add.set.annotation.row = TRUE) {
  
  # Checks
  if (is.null(gene.sets) || !is.list(gene.sets)) {
    stop("gene_sets must be a named list of character vectors.")
  }
  if (is.null(names(gene.sets)) || any(names(gene.sets) == "")) {
    stop("gene_sets must be a *named* list (names are used for section labels).")
  }
  
  counts.all <- normalized.log.counts.df
  
  # Keep only genes present, preserve order within each set
  present.by.set <- lapply(gene.sets, function(gs) {
    gs <- unique(as.character(gs))
    gs[gs %in% rownames(counts_all)]
  })
  
  if (!keep.empty.sets) {
    present.by.set <- present.by.set[lengths(present.by.set) > 0]
  }
  
  genes.present <- unlist(present.by.set, use.names = FALSE)
  
  # Handle duplicates across sets: keep first occurrence only
  genes.present <- genes.present[!duplicated(genes.present)]
  
  if (length(genes.present) < min.genes.to.display) {
    stop(
      paste0(
        "Only ", length(genes_present),
        " genes from gene_sets were found in normalized.log.counts.df rownames. ",
        "Check gene symbols / rownames."
      )
    )
  }
  
  # Cap number of genes displayed (keeps ordering)
  if (length(genes_present) > max.genes.to.display) {
    genes.present <- genes.present[1:max.genes.to.display]
    present.by.set <- lapply(present.by.set, function(v) v[v %in% genes.present])
    present.by.set <- present.by.set[lengths(present.by.set) > 0]
  }
  
  # Subset & order counts
  counts <- counts.all[genes.present, , drop = FALSE]
  
  # Create row gaps between sets (if not provided)
  if (is.null(row.gaps)) {
    set.sizes <- lengths(present.by.set)
    # pheatmap gaps.row is indices after which to draw a gap
    row.gaps <- cumsum(set.sizes)[-length(set.sizes)]
    if (length(row.gaps) == 0) row.gaps <- NULL
  }
  
  # Add gene set annotation row (optional) 
  if (add.set.annotation.row) {
    set.name.vec <- rep(names(present.by.set), times = lengths(present.by.set))
    set.anno <- data.frame(GeneSet = set.name.vec, row.names = genes.present)
    
    if (is.null(annotation.row)) {
      annotation.row <- set.anno
    } else {
      # Ensure rownames align; subset/reorder if needed
      annotation.row <- annotation.row[rownames(annotation.row) %in% genes.present, , drop = FALSE]
      annotation.row <- annotation.row[match(genes.present, rownames(annotation.row)), , drop = FALSE]
      annotation.row <- cbind(annotation.row, set.anno)
    }
  } else if (!is.null(annotation.row)) {
    # Ensure annotation.row aligns with counts
    annotation.row <- annotation.row[rownames(annotation.row) %in% genes.present, , drop = FALSE]
    annotation.row <- annotation.row[match(genes.present, rownames(annotation.row)), , drop = FALSE]
  }
  
  # ---- Arrange by annotations if no column clustering ----
  if (cluster.columns == FALSE) {
    
    anno.col.names <- colnames(annotation.column)
    
    for (col in anno.col.names) {
      annotation.column <- annotation.column %>%
        dplyr::arrange(.data[[col]])
    }
    
    # Match the counts columns to the annotation row order
    counts <- counts[, rownames(annotation.column), drop = FALSE]
  }
  
  # ---- Plot ----
  heatmap.plot <- pheatmap::pheatmap(
    counts,
    main = main.title,
    show_rownames = show.rownames,
    scale = "row",
    show_colnames = show.colnames,
    border_color = NA,
    cluster_rows = cluster.rows,
    cluster_cols = cluster.columns,
    clustering_method = "average",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    color = colorRampPalette(c("blue", "white", "red"))(120),
    annotation_row = annotation.row,
    annotation_col = annotation.column,
    annotation_colors = anno.colors,
    gaps_row = row.gaps,
    gaps_col = column.gaps,
    fontsize_row = font.size.row
  )
  
  return(heatmap.plot)
}


calculate_signal2noise <- function(){
  
  
}


normalize_counts <- function() {}

gsea_preranked_list <- function(contrast.field, 
                                contrast.levels, 
                                annotation, 
                                log.counts){
  
  # Gather the signal to noise ratio for GSEA ranking
  # Default method for ranking genes from GSEA manual:
  # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
  
  # Contrast level A is the "condition" (positive when calculating fold change)
  contrast.A.annotation <- annotation %>% 
    filter(!!sym(contrast.field) == contrast.levels[1])
  
  contrast.A.sampleIDs <- rownames(contrast.A.annotation)
  
  contrast.A.counts <- as.data.frame(log.counts) %>% 
    dplyr::select(all_of(contrast.A.sampleIDs))
  
  contrast.A.counts$gene <- rownames(contrast.A.counts)
  
  # Contrast level B is the "reference" (negative when calculating fold change)
  
  contrast.B.annotation <- annotation %>% 
    filter(!!sym(contrast.field) == contrast.levels[2])
  
  contrast.B.sampleIDs <- rownames(contrast.B.annotation)
  
  contrast.B.counts <- as.data.frame(log.counts) %>% 
    dplyr::select(all_of(contrast.B.sampleIDs))
  
  contrast.B.counts$gene <- rownames(contrast.B.counts)
  
  # Add a column to each contrast level for the mean and standard deviation
  contrast.A.counts <- contrast.A.counts %>% 
    mutate(mean.A = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.A = apply(select_if(., is.numeric), 1, sd))
  
  contrast.B.counts <- contrast.B.counts %>% 
    mutate(mean.B = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.B = apply(select_if(., is.numeric), 1, sd))
  
  GSEA.preanked.df <- merge(contrast.A.counts, contrast.B.counts, by = "gene")
  
  GSEA.preanked.df <- GSEA.preanked.df %>% 
    mutate(signal2noise = (mean.A - mean.B)/(stdev.A + stdev.B)) %>% 
    arrange(desc(signal2noise)) %>% 
    dplyr::select(c(gene, mean.A, mean.B, stdev.A, stdev.B, signal2noise))
  
  return(GSEA.preanked.df)
  
}

make_volcano <- function(lmm.results, 
                         title, 
                         legend.title, 
                         x.axis.title, 
                         fc.limit = 1, 
                         pos.label.limit = 1, 
                         neg.label.limit = -1, 
                         custom.gene.labels = NULL){ 
  
  ## Make a volcano plot for the comparison
  
  # Create a column for direction of DEGs
  lmm.results$de_direction <- "NONE"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc > fc.limit] <- "UP"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc < -fc.limit] <- "DOWN"
  
  # Create a label for DEGs based on label limits
  lmm.results$deglabel <- ifelse((lmm.results$logfc > pos.label.limit | 
                                    lmm.results$logfc < neg.label.limit) & 
                                   lmm.results$padj < 0.05, 
                                 lmm.results$gene,
                                 NA
  )
  
  # Create a label for DEGs
  if(is.null(custom.gene.labels)){
    
    lmm.results$deglabel <- ifelse(lmm.results$de_direction == "NONE", 
                                   NA, 
                                   lmm.results$gene)
    
  } else {
    
    lmm.results$deglabel <- ifelse(lmm.results$gene %in% custom.gene.labels, 
                                   lmm.results$gene, 
                                   NA)
    
  }
  
  # Compute the scale for the volcano x-axis
  log2.scale <- max(abs(lmm.results$logfc))
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("steelblue4", "grey", "violetred4")
  names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
  
  # Make the volcano plot based on custom gene labels
  if(is.null(custom.gene.labels)){
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
    
    volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                                   y = -log10(padj), 
                                                   col = de_direction, 
                                                   label = deglabel)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors) + 
      geom_text_repel(max.overlaps = Inf, 
                      show.legend = FALSE) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    
    # Label the custom genes depending on significance
    lmm.results <- lmm.results %>% 
      mutate(custom.label = ifelse(!is.na(deglabel) & de_direction == "NONE", 
                                   "BLACK", 
                                   ifelse(!is.na(deglabel) & de_direction != "NONE", 
                                          de_direction, 
                                          "NONE")))
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4", "black")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP", "BLACK")
    
    lmm.results.labeled <- lmm.results %>%
      filter(custom.label != "NONE")
    
    lmm.results.unlabeled <- lmm.results %>% 
      filter(custom.label == "NONE")
    
    
    volcano.plot <- ggplot() + 
      geom_point(data = lmm.results.unlabeled, aes(x = logfc, 
                                                   y = -log10(padj), 
                                                   col = custom.label, 
                                                   alpha = 0.5)) + 
      geom_point(data = lmm.results.labeled, aes(x = logfc, 
                                                 y = -log10(padj), 
                                                 col = custom.label, 
                                                 alpha = 1)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors, 
                         breaks = c("DOWN", "UP")) + 
      geom_text_repel(data = lmm.results.labeled,
                      aes(x = logfc, 
                          y = -log10(padj), 
                          label = deglabel, 
                          col = custom.label), 
                      max.overlaps = Inf, 
                      size = 6, 
                      show.legend = FALSE) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_alpha_identity(guide = "none")
    
  }
  
  return(list("volcano.plot" = volcano.plot))
  
}



improved_make_volcano <- function(lmm.results, 
                         title, 
                         title.size = 16,
                         legend.title, 
                         fc.limit = 1, 
                         custom.gene.labels = NULL, 
                         remove.controls = FALSE, 
                         remove.genes = NULL, 
                         remove.all.gene.labels = FALSE,
                         legend.coordinates = c(.99, .01), 
                         x.lab = "Log2 Fold Change", 
                         y.lab = "-Log10 adjusted p-value", 
                         dotted.line.color = "gray", 
                         alpha = 1, 
                         nonDE.color = "gray", 
                         upDE.color = "violetred4", 
                         downDE.color = "steelblue4", 
                         label.size = NULL, 
                         label.color = "custom", 
                         axis.tick.label.size = 8, 
                         legend.text.size = 8, 
                         axis.title.size = 8){ 
  
  # Ensure that titles are characters
  legend.title <- as.character(legend.title)
  
  # Remove controls if applicable
  #if(remove.controls == TRUE){
#    
#    NEG.indices <- grep("NEG_", rownames(de.results))
#    POS.indices <- grep("POS_", rownames(de.results))
#    control.indices <- c(NEG.indices, POS.indices)
#    
#    # Remove the control probes
#    de.results <- de.results[-control.indices,]
#    
#  } 
  
  # Remove custom gene input
  if(!is.null(remove.genes)){
    
    lmm.results <- lmm.results[!rownames(lmm.results) %in% remove.genes, ]
    
  }
  
  
  
  # Create a column for direction of DEGs
  lmm.results$de_direction <- "NONE"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc > fc.limit] <- "UP"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc < -fc.limit] <- "DOWN"
  
  # Create a label for DEGs
  if(is.null(custom.gene.labels)){
    
    lmm.results$deglabel <- ifelse(lmm.results$de_direction == "NONE", 
                                  NA, 
                                  lmm.results$gene)
    
  } else {
    
    lmm.results$deglabel <- ifelse(lmm.results$gene %in% custom.gene.labels, 
                                  lmm.results$gene, 
                                  NA)
    
  }
  
  # Remove all gene labels, in the case of too many DEGs
  if(remove.all.gene.labels){
    
    lmm.results$deglabel <- NA
    
  }
  
  
  
  # Compute the scale for the volcano x-axis
  log2.scale <- ceiling(max(abs(lmm.results$logfc)))
  
  # Establish the color scheme for the volcano plot
  if(is.null(custom.gene.labels)){
    
    contrast.level.colors <- c(downDE.color, nonDE.color, upDE.color)
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
    
    
    
    if(is.null(legend.coordinates)){ 
      
      volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                                    y = -log10(padj), 
                                                    col = de_direction, 
                                                    label = deglabel)) +
        geom_vline(xintercept = c(-fc.limit, fc.limit), col = dotted.line.color, linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = dotted.line.color, linetype = 'dashed') + 
        labs(x = x.lab,
             y = y.lab, 
             title = title) + 
        geom_point(size = 2, alpha = alpha) +
        scale_color_manual(legend.title, 
                           values = contrast.level.colors, 
                           breaks = "UP", "NONE", "DOWN") + 
        geom_text_repel(max.overlaps = Inf, 
                        show.legend = FALSE, 
                        size = label.size, 
                        color = "black") + 
        xlim(-log2.scale-0.2, log2.scale+0.2) + 
        theme_classic() + 
        theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
              plot.title = element_text(hjust = 0.5, 
                                        face = "bold", 
                                        size = title.size), 
              axis.title = element_text(face = "bold", 
                                        size = axis.title.size), 
              axis.text = element_text(size = axis.tick.label.size), 
              legend.position="inside", 
              legend.position.inside = legend.coordinates, 
              legend.justification = c(1, 0), 
              legend.box.background = element_rect(
                colour = "black", 
                linewidth = 0.5,  
                fill = "white"), 
              legend.text = element_text(face = "bold", 
                                         size = legend.text.size), 
              legend.title = element_text(face = "bold", 
                                          size = legend.text.size)) + 
        scale_x_continuous(breaks = scales::breaks_pretty(n = 7),
                           limits = c(-log2.scale, log2.scale))
      
    } else {
      
      
      volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                                    y = -log10(padj), 
                                                    col = de_direction, 
                                                    label = deglabel)) +
        geom_vline(xintercept = c(-fc.limit, fc.limit), col = dotted.line.color, linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = dotted.line.color, linetype = 'dashed') + 
        labs(x = x.lab,
             y = y.lab, 
             title = title) + 
        geom_point(size = 2, alpha = alpha) +
        scale_color_manual(legend.title, 
                           values = contrast.level.colors, 
                           breaks = c("UP", "NONE", "DOWN")) + 
        geom_text_repel(max.overlaps = Inf, 
                        show.legend = FALSE, 
                        size = label.size, 
                        color = "black") + 
        theme_classic() + 
        theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
              plot.title = element_text(hjust = 0.5, 
                                        face = "bold", 
                                        size = title.size), 
              axis.title = element_text(face = "bold", 
                                        size = axis.title.size), 
              axis.text = element_text(size = axis.tick.label.size), 
              legend.position="inside", 
              legend.position.inside = legend.coordinates, 
              legend.justification = c(1, 0), 
              legend.box.background = element_rect(
                colour = "black", 
                linewidth = 0.5,  
                fill = "white"), 
              legend.text = element_text(face = "bold", 
                                         size = legend.text.size), 
              legend.title = element_text(face = "bold", 
                                          size = legend.text.size)) + 
        scale_x_continuous(breaks = scales::breaks_pretty(n = 7),
                           limits = c(-log2.scale, log2.scale))
      
      
    }
    
    
    
    
    
  } else {
    
    # Label the custom genes depending on significance
    lmm.results <- lmm.results %>% 
      mutate(custom.label = ifelse(!is.na(deglabel) & de_direction == "NONE", 
                                   "BLACK", 
                                   ifelse(!is.na(deglabel) & de_direction != "NONE", 
                                          de_direction, 
                                          "NONE")))
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4", "black")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP", "BLACK")
    
    lmm.results.labeled <- lmm.results %>%
      filter(custom.label != "NONE")
    
    lmm.results.unlabeled <- lmm.results %>% 
      filter(custom.label == "NONE")
    
    
    volcano.plot <- ggplot() + 
      geom_point(data = lmm.results.unlabeled, aes(x = logfc, 
                                                  y = -log10(padj), 
                                                  col = custom.label, 
                                                  alpha = 0.5)) + 
      geom_point(data = lmm.results.labeled, aes(x = logfc, 
                                                y = -log10(padj), 
                                                col = custom.label, 
                                                alpha = 1)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      labs(x = "log2 Fold Change",
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors, 
                         breaks = c("DOWN", "UP")) + 
      geom_text_repel(data = lmm.results.labeled,
                      aes(x = logfc, 
                          y = -log10(padj), 
                          label = deglabel, 
                          col = custom.label), 
                      max.overlaps = Inf, 
                      show.legend = FALSE) + 
      scale_alpha_identity(guide = "none") + 
      theme_classic() + 
      theme(plot.title = element_text(hjust = 0.5, 
                                      face = "bold", 
                                      size = title.size), 
            axis.title = element_text(face = "bold", 
                                      size = axis.title.size), 
            axis.text = element_text(size = axis.tick.label.size), 
            legend.position="inside", 
            legend.position.inside = legend.coordinates, 
            legend.justification = c(1, 0), 
            legend.box.background = element_rect(
              colour = "black", 
              linewidth = 0.5,  
              fill = "white"), 
            legend.text = element_text(face = "bold", 
                                       size = legend.text.size), 
            legend.title = element_text(face = "bold", 
                                        size = legend.text.size)) + 
      scale_x_continuous(breaks = scales::breaks_pretty(n = 7),
                         limits = c(-log2.scale, log2.scale))
    
  }
  
  
  return(volcano.plot)
  
} 

make_dual_volcano <- function(lmm.results.1, 
                              lmm.results.2,
                              lmm.results.1.label,
                              lmm.results.2.label, 
                              title, 
                              legend.title, 
                              x.axis.title, 
                              fc.limit = 1, 
                              pos.label.limit = 1, 
                              neg.label.limit = -1, 
                              custom.gene.labels = NULL, 
                              legend.coordinates = NULL){
  
  # A combined list of the two results
  lmm.results.list <- list()
  lmm.results.list[[lmm.results.1.label]] <- lmm.results.1
  lmm.results.list[[lmm.results.2.label]] <- lmm.results.2
  
  # Create a column for direction of DEGs
  for(lmm.results.label in names(lmm.results.list)){
    
    lmm.results <- lmm.results.list[[lmm.results.label]]
    
    lmm.results$de_direction <- "NONE"
    lmm.results$de_direction[lmm.results$padj < 0.05 & 
                               lmm.results$logfc > fc.limit] <- "UP"
    lmm.results$de_direction[lmm.results$padj < 0.05 & 
                               lmm.results$logfc < -fc.limit] <- "DOWN"
    
    #lmm.results$de_direction <- factor(lmm.results$de_direction, 
    #                                   levels = "UP", "NONE", "DOWN")
    
    # Create a label for DEGs based on label limits
    lmm.results$deglabel <- ifelse((lmm.results$logfc > pos.label.limit | 
                                      lmm.results$logfc < neg.label.limit) & 
                                     lmm.results$padj < 0.05, 
                                   lmm.results$gene,
                                   NA)
    
    # Convert to -log10 p-value
    lmm.results$neg.log10.pval <- -log10(lmm.results$padj)
    
    # Add an analysis label 
    lmm.results$analysis.label <- lmm.results.label
    
    lmm.results.list[[lmm.results.label]] <- lmm.results
    
  }
  
  # Convert -log10 p-value to neg for mirrored second analysis
  lmm.results.list[[lmm.results.2.label]]$neg.log10.pval <- -(lmm.results.list[[lmm.results.2.label]]$neg.log10.pval)
  
  # Combine the two analyses into a master df
  lmm.results.combine <- bind_rows(lmm.results.list[[lmm.results.1.label]], 
                                   lmm.results.list[[lmm.results.2.label]])
  
  # Establish the limits and breaks for x and y axes
  log2.scale <- ceiling(max(abs(lmm.results$logfc)))
  pval.scale <- ceiling(max(abs(lmm.results.combine$neg.log10.pval)))

  y.axis.breaks <- seq(-pval.scale, pval.scale, by = 1)
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("violetred4", "grey", "steelblue4")
  names(contrast.level.colors) <- c("UP", "NONE", "DOWN")
  
  
  
  # Make the plot
  if(is.null(legend.coordinates)){
    
    dual.volcano.plot <- ggplot(data = lmm.results.combine, 
                                aes(x = logfc, 
                                    y = neg.log10.pval, 
                                    col = de_direction, 
                                    label = deglabel)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), 
                 col = "darkgray", 
                 linetype = 'dashed') + 
      geom_vline(xintercept = 0, 
                 col = "black") + 
      geom_hline(yintercept = c(-log10(0.05), -(-log10(0.05))), 
                 col = "darkgray", 
                 linetype = 'dashed') + 
      geom_hline(yintercept = 0, 
                 col = "black") + 
      xlim(-log2.scale - 0.2, log2.scale + 0.2) + 
      scale_y_continuous(
        limits = c(-pval.scale, pval.scale),
        breaks = y.axis.breaks,
        labels = abs(y.axis.breaks)) + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors) + 
      geom_text_repel(max.overlaps = Inf, show.legend = FALSE) + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.title = element_text(face = "bold"), 
            axis.title = element_text(face = "bold")) + 
      theme_linedraw()
    
  } else {
    
    
    dual.volcano.plot <- ggplot(data = lmm.results.combine, 
                                aes(x = logfc, 
                                    y = neg.log10.pval, 
                                    col = de_direction, 
                                    label = deglabel)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), 
                 col = "darkgray", 
                 linetype = 'dashed') + 
      #geom_vline(xintercept = 0, 
      #           col = "black", 
      #           linetype = 'dashed') + 
      geom_hline(yintercept = c(-log10(0.05), -(-log10(0.05))), 
                 col = "darkgray", 
                 linetype = 'dashed') + 
      scale_y_continuous(
        limits = c(-pval.scale, pval.scale),
        breaks = y.axis.breaks,
        labels = abs(y.axis.breaks)) + 
      scale_x_continuous(limits = c(-log2.scale, log2.scale),
                         breaks = seq(-log2.scale, log2.scale, length.out = 7)) + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2, alpha = 0.7) +
      geom_hline(yintercept = 0, 
                 col = "black", 
                 linewidth = 0.5) + 
      scale_color_manual(legend.title, 
                         values = contrast.level.colors) + 
      geom_text_repel(max.overlaps = Inf, show.legend = FALSE) + 
      theme_linedraw() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position="inside", 
            legend.position.inside = legend.coordinates, 
            legend.justification = c(1, 1), 
            legend.box.background = element_rect(
              colour = "black", 
              linewidth = 0.5,  
              fill = "white"
            ), 
            legend.title = element_text(face = "bold"), 
            axis.title = element_text(face = "bold"), 
            panel.grid = element_blank())
      
    
  }
  
  
  
  #dual.volcano.labeled <- ggdraw() +
  #  draw_plot(dual.volcano.plot, x = 0.02, width = 0.95) +  # Move the plot slightly right
  #  draw_text(lmm.results.1.label, x = 0.01, y = 0.6, angle = 90, size = 14, hjust = 0) +
  #  draw_text(lmm.results.2.label, x = 0.01, y = 0.2, angle = 90, size = 14, hjust = 0)
  
  return(dual.volcano.plot)

}

default_if_null <- function(x, default) {
  if (is.null(x)) default else x
}

# Human-readable label for Quarto output.
de_method_label <- function(de.method) {
  if (de.method == "default") {
    "Default DSPWorkflow diffExpr"
  } else if (de.method == "standr") {
    "StandR limma-voom"
  } else {
    de.method
  }
}



# DE result formatting helpers


# Convert DSPWorkflow::diffExpr() output into the common format expected by the
# downstream volcano and DEG counting code:
#
#   gene, logfc, pval, padj

standardize_default_results <- function(results.df) {
  # diffExpr() creates method-specific column names, so identify them by pattern.
  logfc.column <- grep("logFC", colnames(results.df), value = TRUE)
  pval.column <- grep("_pval", colnames(results.df), value = TRUE)
  adj.pval.column <- grep("adjpval", colnames(results.df), value = TRUE)
  
  # These checks make the function fail clearly if diffExpr() output changes.
  if (length(logfc.column) != 1) {
    stop(
      "Expected exactly one logFC column from diffExpr, found: ",
      paste(logfc.column, collapse = ", ")
    )
  }
  
  if (length(pval.column) != 1) {
    stop(
      "Expected exactly one _pval column from diffExpr, found: ",
      paste(pval.column, collapse = ", ")
    )
  }
  
  if (length(adj.pval.column) != 1) {
    stop(
      "Expected exactly one adjpval column from diffExpr, found: ",
      paste(adj.pval.column, collapse = ", ")
    )
  }
  
  # Use Gene if available; otherwise fall back to gene or rownames.
  gene <- if ("Gene" %in% colnames(results.df)) {
    results.df$Gene
  } else if ("gene" %in% colnames(results.df)) {
    results.df$gene
  } else {
    rownames(results.df)
  }
  
  data.frame(
    gene = gene,
    logfc = results.df[[logfc.column]],
    pval = results.df[[pval.column]],
    padj = results.df[[adj.pval.column]],
    stringsAsFactors = FALSE
  )
}


# Convert limma::topTable() output from the StandR workflow into the same common
# format used by the rest of the pipeline:
#
#   gene, logfc, pval, padj
standardize_standr_results <- function(results.df) {
  if (!"logFC" %in% colnames(results.df)) {
    stop("StandR/limma results did not contain a logFC column.")
  }
  
  if (!"P.Value" %in% colnames(results.df)) {
    stop("StandR/limma results did not contain a P.Value column.")
  }
  
  if (!"adj.P.Val" %in% colnames(results.df)) {
    stop("StandR/limma results did not contain an adj.P.Val column.")
  }
  
  # topTable() often stores genes as rownames.
  gene <- if ("gene" %in% colnames(results.df)) {
    results.df$gene
  } else if ("Gene" %in% colnames(results.df)) {
    results.df$Gene
  } else {
    rownames(results.df)
  }
  
  standardized.df <- data.frame(
    gene = gene,
    logfc = results.df$logFC,
    pval = results.df$P.Value,
    padj = results.df$adj.P.Val,
    stringsAsFactors = FALSE
  )
  
  # Preserve average expression if present.
  mean.expr.col <- intersect(
    c("mean_expr", "AveExpr", "Amean"),
    colnames(results.df)
  )
  
  if (length(mean.expr.col) > 0) {
    standardized.df$mean_expr <- results.df[[mean.expr.col[1]]]
  }
  
  standardized.df
}


# When a DE result CSV already exists, read it back in and make sure it has the
# columns needed by downstream plotting/export code.
validate_cached_results <- function(results.df, result.file) {
  required.cols <- c("gene", "logfc", "pval", "padj")
  missing.cols <- setdiff(required.cols, colnames(results.df))
  
  if (length(missing.cols) > 0) {
    stop(
      "Cached DE result file does not contain the expected columns: ",
      result.file,
      "\nMissing columns: ",
      paste(missing.cols, collapse = ", ")
    )
  }
  
  # Keep required columns first, then preserve any extra columns.
  # Drop the common accidental CSV row-number column X.
  extra.cols <- setdiff(colnames(results.df), c(required.cols, "X"))
  
  results.df <- results.df[
    ,
    c(required.cols, extra.cols),
    drop = FALSE
  ]
  
  # Ensure numeric fields are numeric after reading from CSV.
  results.df$logfc <- as.numeric(results.df$logfc)
  results.df$pval <- as.numeric(results.df$pval)
  results.df$padj <- as.numeric(results.df$padj)
  
  results.df
}



# Volcano plot helpers


# Identify a limited number of up- and down-regulated genes to label on volcano
# plots. This prevents overcrowding when many genes are significant.
get_top_gene_labels <- function(results.df,
                                fc.limit = 1.5,
                                padj.limit = 0.05,
                                n = 10) {
  down <- results.df %>%
    dplyr::filter(!is.na(logfc), !is.na(padj)) %>%
    dplyr::filter(logfc < -fc.limit, padj < padj.limit) %>%
    dplyr::arrange(padj, logfc)
  
  up <- results.df %>%
    dplyr::filter(!is.na(logfc), !is.na(padj)) %>%
    dplyr::filter(logfc > fc.limit, padj < padj.limit) %>%
    dplyr::arrange(padj, dplyr::desc(logfc))
  
  unique(c(
    head(down$gene, n),
    head(up$gene, n)
  ))
}


# Wrapper around improved_make_volcano() function.
# This keeps volcano parameter handling in one place.
make_de_volcano <- function(results.df,
                            comp.params,
                            custom.labels = NULL,
                            remove.all.gene.labels = FALSE,
                            label.size = NULL) {
  improved_make_volcano(
    lmm.results = results.df,
    title = default_if_null(
      comp.params$volcano.title,
      default_if_null(comp.params$heading, comp.params$contrast.name)
    ),
    title.size = default_if_null(comp.params$title.size, 22),
    legend.title = default_if_null(comp.params$legend.title, "Expression"),
    fc.limit = default_if_null(comp.params$fc.limit, 1.5),
    custom.gene.labels = custom.labels,
    remove.controls = default_if_null(comp.params$remove.controls, FALSE),
    remove.genes = default_if_null(comp.params$remove.genes, NULL),
    remove.all.gene.labels = remove.all.gene.labels,
    legend.coordinates = default_if_null(comp.params$legend.coordinates, c(.99, .8)),
    x.lab = default_if_null(comp.params$x.lab, "Log2 Fold Change"),
    y.lab = default_if_null(comp.params$y.lab, "-Log10 adjusted p-value"),
    dotted.line.color = default_if_null(comp.params$dotted.line.color, "black"),
    alpha = default_if_null(comp.params$alpha, 1),
    nonDE.color = default_if_null(comp.params$nonDE.color, "gray60"),
    upDE.color = default_if_null(comp.params$upDE.color, "red"),
    downDE.color = default_if_null(comp.params$downDE.color, "blue"),
    label.size = default_if_null(
      label.size,
      default_if_null(comp.params$label.size, 8)
    ),
    label.color = default_if_null(comp.params$label.color, "black"),
    axis.tick.label.size = default_if_null(comp.params$axis.tick.label.size, 20),
    legend.text.size = default_if_null(comp.params$legend.text.size, 20),
    axis.title.size = default_if_null(comp.params$axis.title.size, 20)
  )
}



# StandR subsetting, design, and contrast helpers


# Given an annotation data frame and a list of fields/values, return TRUE for
# samples that match all requested fields.
#
# Example group.params:
#   list(class = "DKD", region = "tubule", segment = "PanCK_pos")
annotation_matches_group <- function(annotation.df, group.params) {
  if (is.null(group.params) || length(group.params) == 0) {
    return(rep(TRUE, nrow(annotation.df)))
  }
  
  keep <- rep(TRUE, nrow(annotation.df))
  
  for (field in names(group.params)) {
    if (!field %in% colnames(annotation.df)) {
      stop(
        "Field '", field, "' was not found in StandR colData.\n\n",
        "Available colData fields are:\n",
        paste(colnames(annotation.df), collapse = "\n")
      )
    }
    
    keep <- keep & annotation.df[[field]] %in% group.params[[field]]
  }
  
  keep
}


# Subset the StandR object to the AOIs/samples used in the current comparison.
# This mirrors the subset.group.1/subset.group.2 logic used for the default
# DSPWorkflow analysis.
subset_standr_object_for_comparison <- function(standr.object, comp.params) {
  annotation.df <- as.data.frame(colData(standr.object))
  
  keep.group.1 <- annotation_matches_group(
    annotation.df = annotation.df,
    group.params = comp.params$subset.group.1
  )
  
  keep.group.2 <- annotation_matches_group(
    annotation.df = annotation.df,
    group.params = comp.params$subset.group.2
  )
  
  # Keep samples belonging to either contrast group.
  keep <- keep.group.1 | keep.group.2
  
  if (!any(keep)) {
    stop(
      "No StandR samples matched subset.group.1 or subset.group.2 for contrast: ",
      comp.params$contrast.name
    )
  }
  
  standr.subset.object <- standr.object[, keep]
  
  cat(
    "**StandR samples retained:** ",
    sum(keep),
    " / ",
    length(keep),
    "\n\n",
    sep = ""
  )
  
  standr.subset.object
}


# Build the StandR/limma model matrix.
#
# Example:
#   design.annotation = "segment"
#   covariates = c("ruv_W1", "ruv_W2")
#
# becomes:
#   model.matrix(~0 + segment + ruv_W1 + ruv_W2, data = colData(...))
#
# The design annotation prefix is removed so columns like segmentPanCK_pos become
# PanCK_pos, making contrast definitions cleaner.
make_standr_design <- function(standr.object, comp.params) {
  standr.params <- comp.params$standr
  
  object.metadata <- as.data.frame(colData(standr.object))
  object.metadata <- droplevels(object.metadata)
  
  design.annotation <- standr.params$design.annotation
  covariates <- default_if_null(standr.params$covariates, NULL)
  
  if (is.null(design.annotation)) {
    stop(
      "standr$design.annotation is NULL for contrast: ",
      comp.params$contrast.name
    )
  }
  
  # Make sure all requested design fields exist in colData.
  required.fields <- c(design.annotation, covariates)
  missing.fields <- setdiff(required.fields, colnames(object.metadata))
  
  if (length(missing.fields) > 0) {
    stop(
      "The following design fields were not found in StandR colData: ",
      paste(missing.fields, collapse = ", "),
      "\n\nAvailable colData fields are:\n",
      paste(colnames(object.metadata), collapse = "\n")
    )
  }
  
  design.terms <- c(design.annotation, covariates)
  
  # Equivalent to ~0 + design.annotation + covariate1 + covariate2.
  # intercept = FALSE creates one coefficient per annotation level.
  design.formula <- stats::reformulate(
    termlabels = design.terms,
    intercept = FALSE
  )
  
  design <- stats::model.matrix(
    design.formula,
    data = object.metadata
  )
  
  # Clean names like segmentPanCK_pos -> PanCK_pos.
  if (isTRUE(default_if_null(standr.params$clean.design.names, TRUE))) {
    clean.pattern <- default_if_null(
      standr.params$clean.design.pattern,
      paste0("^", design.annotation)
    )
    
    clean.replacement <- default_if_null(
      standr.params$clean.design.replacement,
      ""
    )
    
    colnames(design) <- gsub(
      pattern = clean.pattern,
      replacement = clean.replacement,
      x = colnames(design)
    )
  }
  
  # Duplicate names make contrasts ambiguous, so stop if cleaning caused a clash.
  if (anyDuplicated(colnames(design))) {
    stop(
      "Cleaning the design column names produced duplicate column names:\n",
      paste(colnames(design), collapse = "\n")
    )
  }
  
  if (isTRUE(default_if_null(standr.params$print.design.columns, TRUE))) {
    cat("**StandR design columns**\n\n")
    cat(paste(colnames(design), collapse = ", "))
    cat("\n\n")
  }
  
  design
}


# Build the limma contrast matrix.
#
# Example:
#   contrast.numerator = "PanCK_pos"
#   contrast.denominator = "PanCK_neg"
#
# becomes:
#   PanCK_pos - PanCK_neg
make_standr_contrast <- function(design, comp.params) {
  standr.params <- comp.params$standr
  
  numerator <- standr.params$contrast.numerator
  denominator <- standr.params$contrast.denominator
  
  if (is.null(numerator) || is.null(denominator)) {
    stop(
      "standr$contrast.numerator and standr$contrast.denominator must both be supplied for contrast: ",
      comp.params$contrast.name
    )
  }
  
  # Make sure the requested contrast columns actually exist in the design matrix.
  missing.cols <- setdiff(
    c(numerator, denominator),
    colnames(design)
  )
  
  if (length(missing.cols) > 0) {
    stop(
      "The following contrast columns were not found in the StandR design matrix: ",
      paste(missing.cols, collapse = ", "),
      "\n\nAvailable design columns are:\n",
      paste(colnames(design), collapse = "\n")
    )
  }
  
  contrast.expression <- paste(numerator, "-", denominator)
  
  contrast.matrix <- limma::makeContrasts(
    contrasts = contrast.expression,
    levels = colnames(design)
  )
  
  colnames(contrast.matrix) <- default_if_null(
    standr.params$contrast.name,
    comp.params$contrast.name
  )
  
  contrast.matrix
}



# StandR GSEA log-count helper


# Create the expression matrix that will be passed into gsea_preranked_list()
# for the StandR method.
#
# Why this is needed:
#   gsea_preranked_list() calculates signal-to-noise from log-normalized
#   expression values.
#
# For the default method, you already have:
#   lmm.input$subset.object@assayData$log_q
#
# For StandR, the input object starts from counts, so we use the limma-voom
# expression matrix:
#   v$E
#
# v$E is normalized logCPM. Then, optionally, we remove RUV covariates and/or
# batch effects while preserving the biological contrast of interest.
make_standr_gsea_log_counts <- function(v,
                                        standr.object,
                                        comp.params) {
  standr.params <- comp.params$standr
  
  # Pull sample metadata for the same samples present in the voom matrix.
  object.metadata <- as.data.frame(colData(standr.object))
  
  # Reorder metadata to match the columns of v$E.
  # This matters because removeBatchEffect() assumes rows of metadata correspond
  # to columns of the expression matrix.
  object.metadata <- object.metadata[colnames(v$E), , drop = FALSE]
  object.metadata <- droplevels(object.metadata)
  
  # This is the biological field used by gsea_preranked_list(), for example:
  #   segment, class, or region.
  contrast.field <- default_if_null(
    comp.params$gsea.contrast.field,
    comp.params$region.col
  )
  
  if (!contrast.field %in% colnames(object.metadata)) {
    stop(
      "The GSEA contrast field was not found in StandR colData: ",
      contrast.field,
      "\n\nAvailable colData fields are:\n",
      paste(colnames(object.metadata), collapse = "\n")
    )
  }
  
  # This design tells removeBatchEffect() what biological variation to preserve.
  # Without this, batch/covariate removal could accidentally remove the contrast
  # signal that GSEA is supposed to rank.
  preserve.design <- stats::model.matrix(
    stats::reformulate(contrast.field, intercept = FALSE),
    data = object.metadata
  )
  
  # RUV covariates are usually technical/unwanted variation.
  # For GSEA, we can remove their effect from v$E while preserving the contrast.
  covariate.matrix <- NULL
  
  if (isTRUE(default_if_null(standr.params$remove.ruv.for.gsea, TRUE))) {
    covariates <- default_if_null(standr.params$covariates, NULL)
    
    if (!is.null(covariates)) {
      missing.covariates <- setdiff(
        covariates,
        colnames(object.metadata)
      )
      
      if (length(missing.covariates) > 0) {
        stop(
          "The following StandR covariates were not found in colData: ",
          paste(missing.covariates, collapse = ", ")
        )
      }
      
      covariate.df <- object.metadata[, covariates, drop = FALSE]
      
      # RUV W factors should be numeric. If they are not, stop clearly.
      non.numeric.covariates <- names(covariate.df)[
        !vapply(covariate.df, is.numeric, logical(1))
      ]
      
      if (length(non.numeric.covariates) > 0) {
        stop(
          "The following GSEA covariates are not numeric: ",
          paste(non.numeric.covariates, collapse = ", "),
          "\nremoveBatchEffect(covariates = ...) expects numeric covariates."
        )
      }
      
      covariate.matrix <- as.matrix(covariate.df)
    }
  }
  
  # Optionally remove slide/batch effects from the GSEA expression matrix.
  #
  # Default recommendation:
  #   remove.batch.for.gsea = FALSE
  #
  # Reason:
  #   slide_name can be confounded with biology. If a slide contains only DKD
  #   or only normal samples, removing slide as a fixed batch effect can remove
  #   real biology. The duplicateCorrelation() step in the DE model handles
  #   slide/block correlation separately.
  batch <- NULL
  
  if (isTRUE(default_if_null(standr.params$remove.batch.for.gsea, FALSE))) {
    batch.field <- standr.params$block.field
    
    if (is.null(batch.field)) {
      stop(
        "standr$remove.batch.for.gsea is TRUE, but standr$block.field is NULL."
      )
    }
    
    if (!batch.field %in% colnames(object.metadata)) {
      stop(
        "The requested batch/block field was not found in StandR colData: ",
        batch.field,
        "\n\nAvailable colData fields are:\n",
        paste(colnames(object.metadata), collapse = "\n")
      )
    }
    
    batch <- object.metadata[[batch.field]]
  }
  
  # If there is nothing to remove, return voom-normalized logCPM directly.
  if (is.null(batch) && is.null(covariate.matrix)) {
    return(v$E)
  }
  
  # Remove unwanted technical variation while preserving contrast.field.
  adjusted.log.counts <- limma::removeBatchEffect(
    x = v$E,
    batch = batch,
    covariates = covariate.matrix,
    design = preserve.design
  )
  
  return(adjusted.log.counts)
  
  cat("Created StandR log counts")
}



# DE method runners


# Run your original DSPWorkflow::diffExpr() method.
run_default_de <- function(lmm.input, comp.params) {
  results.list <- diffExpr(
    object = lmm.input$subset.object,
    analysis.type = default_if_null(comp.params$analysis.type, "Within Groups"),
    region.col = comp.params$region.col,
    regions = comp.params$regions,
    group.col = comp.params$group.col,
    groups = comp.params$groups,
    n.cores = default_if_null(comp.params$n.cores, parallel::detectCores())
  )
  
  results.df <- standardize_default_results(results.list$results)
  
  list(
    results = results.df,
    results.raw = results.list
  )
}


# Run the StandR/edgeR/limma-voom method.
run_standr_de <- function(standr.object, comp.params) {
  standr.params <- comp.params$standr
  
  # The StandR object should be loaded outside the function and passed in.
  # Example:
  #   standr.ruv.object <- readRDS(...)
  #   run_all_de_comparisons(..., standr.object = standr.ruv.object)
  if (is.null(standr.object)) {
    stop(
      "StandR was requested, but no standr.object was supplied. ",
      "Load it with readRDS(), then pass standr.object = standr.ruv.object."
    )
  }
  
  # Restrict StandR object to samples in this comparison.
  standr.subset.object <- subset_standr_object_for_comparison(
    standr.object = standr.object,
    comp.params = comp.params
  )
  
  # Convert SummarizedExperiment/SpatialExperiment-style object to edgeR DGEList.
  ruv.dge.object <- SE2DGEList(standr.subset.object)
  
  # Build model matrix, for example:
  #   ~0 + segment + ruv_W1 + ruv_W2
  design <- make_standr_design(
    standr.object = standr.subset.object,
    comp.params = comp.params
  )
  
  # Build contrast matrix, for example:
  #   PanCK_pos - PanCK_neg
  contr.matrix <- make_standr_contrast(
    design = design,
    comp.params = comp.params
  )
  
  # Remove genes that are too lowly expressed for reliable DE testing.
  if (isTRUE(default_if_null(standr.params$filter.by.expr, TRUE))) {
    keep <- edgeR::filterByExpr(
      y = ruv.dge.object,
      design = design
    )
  } else {
    keep <- rep(TRUE, nrow(ruv.dge.object))
    names(keep) <- rownames(ruv.dge.object)
  }
  
  cat(
    "**StandR genes retained after filterByExpr:** ",
    sum(keep),
    " / ",
    length(keep),
    "\n\n",
    sep = ""
  )
  
  if (!any(keep)) {
    stop(
      "No genes were retained by filterByExpr for contrast: ",
      comp.params$contrast.name
    )
  }
  
  removed.genes <- rownames(ruv.dge.object)[!keep]
  
  # Apply gene filter.
  ruv.dge.object.gene.filter <- ruv.dge.object[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  # Estimate dispersion using edgeR before voom.
  ruv.dge.object.gene.filter <- edgeR::estimateDisp(
    y = ruv.dge.object.gene.filter,
    design = design,
    robust = default_if_null(standr.params$estimate.disp.robust, TRUE)
  )
  
  object.metadata <- as.data.frame(colData(standr.subset.object))
  
  # Optional blocking variable for repeated or correlated samples.
  # In your case this is often slide_name.
  block <- NULL
  
  if (!is.null(standr.params$block.field)) {
    block.field <- standr.params$block.field
    
    if (!block.field %in% colnames(object.metadata)) {
      stop(
        "The requested StandR block.field was not found in colData: ",
        block.field,
        "\n\nAvailable colData fields are:\n",
        paste(colnames(object.metadata), collapse = "\n")
      )
    }
    
    block <- object.metadata[[block.field]]
    
    if (any(is.na(block))) {
      stop(
        "The StandR block field contains NA values: ",
        block.field
      )
    }
    
    # duplicateCorrelation only makes sense if at least one block has repeated
    # samples. If every block occurs once, run without duplicateCorrelation.
    if (!any(duplicated(block))) {
      warning(
        "No repeated values were found in block.field = '",
        block.field,
        "'. Running the StandR limma model without duplicateCorrelation."
      )
      
      block <- NULL
    }
  }
  
  plot.voom <- default_if_null(standr.params$plot.voom, FALSE)
  
  if (!is.null(block)) {
    # First voom pass without correlation.
    # This estimates the mean-variance trend.
    v0 <- limma::voom(
      ruv.dge.object.gene.filter,
      design = design,
      plot = plot.voom
    )
    
    # Estimate within-block correlation from the first voom object.
    corfit1 <- limma::duplicateCorrelation(
      object = v0,
      design = design,
      block = block
    )
    
    cat(
      "**StandR duplicateCorrelation, first pass:** ",
      round(corfit1$consensus.correlation, 4),
      "\n\n",
      sep = ""
    )
    
    # Second voom pass using the first-pass correlation estimate.
    v <- limma::voom(
      ruv.dge.object.gene.filter,
      design = design,
      block = block,
      correlation = corfit1$consensus.correlation,
      plot = plot.voom
    )
    
    # Re-estimate within-block correlation using the second voom object.
    corfit2 <- limma::duplicateCorrelation(
      object = v,
      design = design,
      block = block
    )
    
    cat(
      "**StandR duplicateCorrelation, second pass:** ",
      round(corfit2$consensus.correlation, 4),
      "\n\n",
      sep = ""
    )
    
    # Fit linear model using final voom object and second-pass correlation.
    fit <- limma::lmFit(
      v,
      design = design,
      block = block,
      correlation = corfit2$consensus.correlation
    )
  } else {
    # If there is no valid block variable, run standard voom + lmFit.
    v <- limma::voom(
      ruv.dge.object.gene.filter,
      design = design,
      plot = plot.voom
    )
    
    corfit1 <- NULL
    corfit2 <- NULL
    
    fit <- limma::lmFit(
      v,
      design = design
    )
  }
  
  # Apply the specified contrast to the fitted model.
  fit.contrast <- limma::contrasts.fit(
    fit = fit,
    contrasts = contr.matrix
  )
  
  # Empirical Bayes moderation.
  efit <- limma::eBayes(
    fit = fit.contrast,
    robust = default_if_null(standr.params$ebayes.robust, TRUE)
  )
  
  # Extract all genes for the contrast.
  results.full <- limma::topTable(
    fit = efit,
    coef = colnames(contr.matrix)[1],
    number = Inf,
    sort.by = default_if_null(standr.params$sort.by, "P")
  )
  
  # Convert limma column names to the common gene/logfc/pval/padj format.
  results.df <- standardize_standr_results(results.full)
  
  # Build StandR-specific expression matrix for GSEA.
  # This is the key update:
  #   gsea_preranked_list() will receive voom-normalized logCPM,
  #   optionally adjusted for RUV and/or batch.
  standr.gsea.log.counts <- make_standr_gsea_log_counts(
    v = v,
    standr.object = standr.subset.object,
    comp.params = comp.params
  )
  
  # Metadata to use with standr.gsea.log.counts.
  # Its rownames must match the sample IDs/column names of gsea.log.counts.
  standr.gsea.annotation <- as.data.frame(colData(standr.subset.object))
  standr.gsea.annotation <- standr.gsea.annotation[
    colnames(standr.gsea.log.counts),
    ,
    drop = FALSE
  ]
  
  out <- list(
    results = results.df,
    results.full = results.full,
    design = design,
    contrast.matrix = contr.matrix,
    keep = keep,
    removed.genes = removed.genes,
    corfit1 = corfit1,
    corfit2 = corfit2,
    gsea.log.counts = standr.gsea.log.counts,
    gsea.annotation = standr.gsea.annotation
  )
  
  # Optionally retain model objects for debugging.
  # These can be large, so default is FALSE.
  if (isTRUE(default_if_null(standr.params$return.fit, FALSE))) {
    out$v <- v
    out$fit <- fit
    out$fit.contrast <- fit.contrast
    out$efit <- efit
  }
  
  if(is.null(out$gsea.log.counts)){
    
    cat("gsea.log.counts not recorded in run_standr_de", "\n\n")
    
  }
  
  out
}



# Output directory helper


# Determine where GSEA input files should be saved.
# This lets your Quarto params use whichever folder field you prefer.
get_gsea_input_dir <- function(params) {
  if (!is.null(params$gsea.input.folder)) {
    return(params$gsea.input.folder)
  }
  
  if (!is.null(params$gsea.folder)) {
    return(file.path(params$gsea.folder, "input_lists"))
  }
  
  if (!is.null(params$results.folder)) {
    return(file.path(params$results.folder, "gsea", "input_lists"))
  }
  
  if (!is.null(params$data.folder)) {
    return(file.path(params$data.folder, "gsea", "input_lists"))
  }
  
  stop(
    "Could not determine where to write GSEA input files. ",
    "Please provide params$gsea.input.folder, params$gsea.folder, ",
    "params$results.folder, or params$data.folder."
  )
}



# Main comparison runner


run_de_comparison <- function(comp.params,
                              normalized.object,
                              params,
                              standr.object = NULL) {
  heading <- default_if_null(
    comp.params$heading,
    comp.params$contrast.name
  )
  
  # Determine which DE method(s) to run for this comparison.
  # Example:
  #   comp.params$de.method = c("default", "standr")
  de.methods <- comp.params$de.method
  
  cat("\n\n#### ", heading, "\n\n", sep = "")
  cat("**DE methods:** ", paste(de.methods, collapse = ", "), "\n\n", sep = "")
  
  # Subset the normalized GeoMx object once.
  # This is used for:
  #   1. the default diffExpr method
  #   2. the sample-count summary table
  #   3. the default-method GSEA signal-to-noise matrix
  lmm.input <- subset_object_for_lmm(
    object = normalized.object,
    subset.group.1 = comp.params$subset.group.1,
    subset.group.2 = comp.params$subset.group.2
  )
  
  # Build a table showing how many samples/AOIs are in the comparison.
  summary.cols <- union(
    names(comp.params$subset.group.1),
    names(comp.params$subset.group.2)
  )
  
  summary.table.df <- pData(lmm.input$subset.object) %>%
    dplyr::select(dplyr::all_of(summary.cols))
  
  summary.table <- as.data.frame(
    table(summary.table.df),
    stringsAsFactors = FALSE
  )
  
  cat("**Sample numbers per annotation group**\n\n")
  cat(knitr::kable(summary.table), sep = "\n")
  cat("\n\n")
  
  # Store method-specific outputs here.
  method.outputs <- stats::setNames(
    vector("list", length(de.methods)),
    de.methods
  )
  

  # Run each requested DE method

  
  for (de.method in de.methods) {
    cat("\n\n##### ", de.method, "\n\n", sep = "")
    
    result.file <- file.path(
      params$de.folder,
      paste0(comp.params$contrast.name, 
             "_", de.method, 
             "_de.results.csv"))
    
    result.RDS <- file.path(params$de.folder, 
                            paste0(comp.params$contrast.name, 
                                   "_", de.method, "_de.results.RDS"))
    
    # Reuse existing DE results unless overwrite.results = TRUE.
    if (
      file.exists(result.RDS) &&
      !isTRUE(default_if_null(comp.params$overwrite.results, FALSE))
    ) {
      
      #results.df <- read.csv(
      #  result.file,
      #  stringsAsFactors = FALSE,
      #  check.names = FALSE)
      
      method.output <- readRDS(result.RDS)
      
      results.df <- method.output$results
      
      cat("*Loaded existing DE results from: ", result.file, "*\n\n", sep = "")
    } else {
      # Run default DSPWorkflow diffExpr.
      if (de.method == "default") {
        method.output <- run_default_de(
          lmm.input = lmm.input,
          comp.params = comp.params
        )
      }
      
      # Run StandR/limma-voom.
      if (de.method == "standr") {
        method.output <- run_standr_de(
          standr.object = standr.object,
          comp.params = comp.params
        )
      }
      
      results.df <- method.output$results
      
      # Save standardized DE results.
      if (isTRUE(default_if_null(comp.params$write.results, TRUE))) {
        dir.create(
          params$de.folder,
          showWarnings = FALSE,
          recursive = TRUE
        )
        
        write.csv(
          results.df,
          file = result.file,
          row.names = FALSE
        )
        
        saveRDS(
          method.output, 
          file = result.RDS
        )
        
        cat("*Wrote DE results to: ", result.file, "*\n\n", sep = "")
      }
    }
    

    # DEG summary

    
    fc.limit <- default_if_null(comp.params$fc.limit, 1.5)
    padj.limit <- default_if_null(comp.params$padj.limit, 0.05)
    
    total.down.reg.degs <- sum(
      results.df$logfc < -fc.limit & results.df$padj < padj.limit,
      na.rm = TRUE
    )
    
    total.up.reg.degs <- sum(
      results.df$logfc > fc.limit & results.df$padj < padj.limit,
      na.rm = TRUE
    )
    
    cat(
      "**DEGs at |log2FC| > ",
      fc.limit,
      " and padj < ",
      padj.limit,
      ":** ",
      total.up.reg.degs,
      " up, ",
      total.down.reg.degs,
      " down\n\n",
      sep = ""
    )
    

    # Volcano plot

    
    custom.labels <- default_if_null(
      comp.params$custom.gene.labels,
      NULL
    )
    
    if (isTRUE(default_if_null(comp.params$label.top.genes, FALSE))) {
      custom.labels <- get_top_gene_labels(
        results.df = results.df,
        fc.limit = fc.limit,
        padj.limit = padj.limit,
        n = default_if_null(comp.params$n.labels, 10)
      )
    }
    
    volcano.plot <- make_de_volcano(
      results.df = results.df,
      comp.params = comp.params,
      custom.labels = custom.labels,
      remove.all.gene.labels = default_if_null(
        comp.params$remove.all.gene.labels,
        FALSE
      ),
      label.size = default_if_null(comp.params$label.size, 8)
    )
    
    print(volcano.plot)
    
    volcano.plot.nolabel <- NULL
    
    if (
      isTRUE(default_if_null(comp.params$make.nolabel.volcano, FALSE)) ||
      isTRUE(default_if_null(comp.params$export.volcano.nolabel, FALSE))
    ) {
      volcano.plot.nolabel <- make_de_volcano(
        results.df = results.df,
        comp.params = comp.params,
        custom.labels = NULL,
        remove.all.gene.labels = TRUE,
        label.size = default_if_null(comp.params$nolabel.label.size, 4)
      )
    }
    
    # Create volcano directory if needed.
    if (
      isTRUE(default_if_null(comp.params$export.volcano, FALSE)) ||
      isTRUE(default_if_null(comp.params$export.volcano.nolabel, FALSE))
    ) {
      volcano.dir <- file.path(params$de.folder, "volcano")
      
      dir.create(
        volcano.dir,
        showWarnings = FALSE,
        recursive = TRUE
      )
    }
    
    if (isTRUE(default_if_null(comp.params$export.volcano, FALSE))) {
      volcano.file <- file.path(
        params$de.folder,
        "volcano",
        paste0(comp.params$contrast.name, "_", 
               de.method, "_volcano_plot.png")
      )
      
      ggsave(
        filename = volcano.file,
        plot = volcano.plot,
        width = default_if_null(comp.params$volcano.width, 14),
        height = default_if_null(comp.params$volcano.height, 10)
      )
      
      cat("*Volcano plot saved to: ", volcano.file, "*\n\n", sep = "")
    }
    
    if (isTRUE(default_if_null(comp.params$export.volcano.nolabel, FALSE))) {
      volcano.nolabel.file <- file.path(
        params$de.folder,
        "volcano",
        paste0(comp.params$contrast.name, "_", de.method, 
               "_volcano_plot_nolabel.png")
      )
      
      ggsave(
        filename = volcano.nolabel.file,
        plot = volcano.plot.nolabel,
        width = default_if_null(comp.params$volcano.width, 14),
        height = default_if_null(comp.params$volcano.height, 10)
      )
      
      cat(
        "*No-label volcano plot saved to: ",
        volcano.nolabel.file,
        "*\n\n",
        sep = ""
      )
    }
    
    # Save method-specific outputs for later GSEA export and return object.
    method.output$method <- de.method
    method.output$result.file <- result.file
    method.output$results <- results.df
    method.output$volcano <- volcano.plot
    method.output$volcano.nolabel <- volcano.plot.nolabel
    
    method.outputs[[de.method]] <- method.output
  }
  

  # GSEA input 
  #
  # This uses your existing gsea_preranked_list() function for both methods.
  #
  # Key difference:
  #   default uses normalized.object log_q
  #   standr uses voom-normalized, optionally RUV/batch-adjusted logCPM

  
  gsea.outputs <- list()
  
  if (isTRUE(default_if_null(comp.params$export.gsea.input, FALSE))) {
    gsea.dir <- get_gsea_input_dir(params)
    
    dir.create(
      gsea.dir,
      showWarnings = FALSE,
      recursive = TRUE
    )
    
    for (de.method in names(method.outputs)) {
      
      if (de.method == "default") {
        # Default GSEA input:
        # use the same normalized log expression matrix you were already using.
        gsea.annotation <- lmm.input$subset.annotation
        gsea.log.counts <- lmm.input$subset.object@assayData$log_q
      }
      
      if (de.method == "standr") {
        
        ### Issue when DE results already exist for StandR. Need to export as a RDS file so that method.outputs can be loaded, which will contain gsea.annotation and gsea.log.counts
        
        # StandR GSEA input:
        # use voom-normalized logCPM generated inside run_standr_de().
        # This may also have RUV covariates removed depending on params.
        gsea.annotation <- method.outputs[[de.method]]$gsea.annotation
        gsea.log.counts <- method.outputs[[de.method]]$gsea.log.counts
        
        if(is.null(gsea.annotation)){
          
          
          cat(paste0(comp.params$heading, " gsea.annotation is NULL for StandR"), 
              "\n\n")
          
        }
        
        if(is.null(gsea.log.counts)){
          
          cat(paste0(comp.params$heading, " gsea.log.counts is NULL for StandR"), 
              "\n\n")
          
        }
        
        if (is.null(gsea.annotation) || is.null(gsea.log.counts)) {
          stop(
            "StandR GSEA annotation/log counts were not found for contrast: ",
            comp.params$contrast.name,
            ". This usually means StandR results were loaded from cache, ",
            "so the voom object was not regenerated. Set overwrite.results = TRUE ",
            "or rerun StandR so gsea.log.counts can be created."
          )
        }
      }
      
      # This is your existing signal-to-noise ranking function.
      gsea.preranked.df <- gsea_preranked_list(
        contrast.field = default_if_null(
          comp.params$gsea.contrast.field,
          comp.params$region.col
        ),
        contrast.levels = default_if_null(
          comp.params$gsea.contrast.levels,
          comp.params$regions
        ),
        annotation = gsea.annotation,
        log.counts = gsea.log.counts
      )
      
      # Save both the preranked data frame and a companion RDS containing
      # annotation and method metadata.
      gsea.input.list <- list(
        preranked_df = gsea.preranked.df,
        annotation_df = gsea.annotation,
        de_method = de.method,
        contrast_name = comp.params$contrast.name,
        contrast_field = default_if_null(
          comp.params$gsea.contrast.field,
          comp.params$region.col
        ),
        contrast_levels = default_if_null(
          comp.params$gsea.contrast.levels,
          comp.params$regions
        )
      )
      
      gsea.csv.file <- file.path(
        gsea.dir,
        paste0(
          comp.params$contrast.name,
          de.method,
          "_gsea_preranked_input.csv"
        )
      )
      
      write.csv(
        gsea.preranked.df,
        file = gsea.csv.file,
        row.names = FALSE
      )
      
      gsea.rds.file <- file.path(
        gsea.dir,
        paste0(
          comp.params$contrast.name,
          de.method,
          "_gsea_input_list.RDS"
        )
      )
      
      saveRDS(
        gsea.input.list,
        file = gsea.rds.file
      )
      
      cat(
        "*GSEA preranked input for ",
        de.method,
        " saved to: ",
        gsea.csv.file,
        "*\n\n",
        sep = ""
      )
      
      cat(
        "*GSEA input list for ",
        de.method,
        " saved to: ",
        gsea.rds.file,
        "*\n\n",
        sep = ""
      )
      
      gsea.outputs[[de.method]] <- list(
        preranked_df = gsea.preranked.df,
        input_list = gsea.input.list,
        csv_file = gsea.csv.file,
        rds_file = gsea.rds.file
      )
    }
  }
  

  # Return all useful objects invisibly

  
  first.method <- de.methods[1]
  
  out <- list(
    contrast.name = comp.params$contrast.name,
    de.methods = de.methods,
    summary.table = summary.table,
    method.results = method.outputs,
    results = method.outputs[[first.method]]$results,
    volcano = method.outputs[[first.method]]$volcano,
    volcano.nolabel = method.outputs[[first.method]]$volcano.nolabel,
    gsea.outputs = gsea.outputs
  )
  
  # Also allow convenient direct access:
  #   de.outputs$contrast_name$default
  #   de.outputs$contrast_name$standr
  for (de.method in names(method.outputs)) {
    out[[de.method]] <- method.outputs[[de.method]]
  }
  
  invisible(out)
}



# Run all comparisons


run_all_de_comparisons <- function(de.comparisons,
                                   normalized.object,
                                   params,
                                   standr.object = NULL) {
  de.outputs <- stats::setNames(
    vector("list", length(de.comparisons)),
    vapply(
      de.comparisons,
      function(x) x$contrast.name,
      character(1)
    )
  )
  
  for (i in seq_along(de.comparisons)) {
    de.outputs[[i]] <- run_de_comparison(
      comp.params = de.comparisons[[i]],
      normalized.object = normalized.object,
      params = params,
      standr.object = standr.object
    )
  }
  
  de.outputs
}
