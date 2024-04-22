subset_for_lmm <- function(object, 
                              subset.list){ 
  
  # Set up the object to subset
  subset.object <- object
  
  # Subset the object based on the given annotations
  for(column in names(subset.list)){ 
    
    subset.indices <- pData(subset.object)[[column]] %in% subset.list[[column]]
    subset.object <- subset.object[, subset.indices]
    
    # Factor the columns with relevant annotations
    pData(subset.object)[[column]] <- factor(pData(subset.object)[[column]])
    
  }
  
  # Factor the slide column
  pData(subset.object)[["slide_name"]] <- 
    factor(pData(subset.object)[["slide_name"]])
  
  # Create log2 counts
  assayDataElement(object = subset.object, elt = "log_q") <-
    assayDataApply(subset.object, 2, FUN = log, base = 2, elt = "q_norm")
  
  # Gather the log counts and annotation to return
  log.counts <- subset.object@assayData$log_q
  annotation.df <- pData(subset.object)
  
  # Replace all bad characters in column names
  annotation.df <- annotation.df %>%
    rename_all(~str_replace_all(., " ", "_"))
  
  return(list("subset.object" = subset.object, 
              "log.counts" = log.counts, 
              "annotation" = annotation.df))
  
}


run_limma <- function(counts, 
                       annotation, 
                       include.slide, 
                       within.slide, 
                       contrast, 
                       contrast.levels){
  
  # Create the DGE object
  DGE.list <- DGEList(counts = counts, 
                      samples = annotation)
  
  if(include.slide == FALSE){ 
    # Create the LM model design
    design <- model.matrix(formula(paste0("~ 0 + ", contrast)), 
                           data = DGE.list$samples)
      
  } else {
    
    if(within.slide == TRUE){ 
      # For within slide we use a random slope in the mixed effect
      
      # Create the LM model design with slide as a mixed effect
      design <- model.matrix(formula(paste0("~ 1 + ", 
                                            contrast, 
                                            " + (1 + " , 
                                            contrast, 
                                            " | slide_name)")), 
                             data = DGE.list$samples)
      
    } else{
      # For between slide we use slide in the mixed effect, no random slope
      
      # Create the LM model design with slide as a mixed effect
      design <- model.matrix(formula(paste0("~ 1 + ", 
                                            contrast, 
                                            " + (1 | slide_name)")), 
                             data = DGE.list$samples)
    }
    
  }
  
  
  # Create the fit for the model
  fit <- lmFit(DGE.list$counts, design)
  
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


make_volcano <- function(lmm.results, 
                         title, 
                         legend.title, 
                         x.axis.title){ 
  
  ## Make a volcano plot for the comparison
  
  # Define the columns for the volcano plot data
  #logfc.column.name <- paste0("logFC_", comparison)
  #padj.column.name <- paste0("adj.pval", comparison)
  
  #results$logfc <- results[[logfc.column.name]]
  #results$padj <- results[[padj.column.name]]
  
  # Create a column for direction of DEGs
  lmm.results$de_direction <- "NONE"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                         lmm.results$logfc > 1.0] <- "UP"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                         lmm.results$logfc < -1.0] <- "DOWN"
  
  # Create a label for DEGs
  lmm.results$deglabel <- ifelse(lmm.results$de_direction == "NONE", 
                                 NA, 
                                 lmm.results$gene)
  
  # Compute the scale for the volcano x-axis
  log2.scale <- max(abs(lmm.results$logfc))
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("steelblue4", "grey", "violetred4")
  names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
  
  # Make the volcano plot
  volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                             y = -log10(padj), 
                                             col = de_direction, 
                                             label = deglabel)) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    xlim(-7.5, 7.5) + 
    labs(x = x.axis.title,
         y = "-log10 adjusted p-value", 
         title = title) + 
    geom_point(size = 2) +
    scale_color_manual(legend.title, 
                       values = contrast.level.colors) + 
    geom_text_repel(max.overlaps = Inf) + 
    xlim(-log2.scale-1, log2.scale+1) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list("volcano.plot" = volcano.plot))
  
  }