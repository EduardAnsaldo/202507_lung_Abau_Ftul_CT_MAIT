   #Visualization Functions
scatterplot <- function (results, group1, group2, local_figures_path, FC_threshold, p_value_threshold, cluster = 'all_clusters', my_colors = c('green4', 'darkorchid4', 'gray'), max_overlaps = 15, label_size = 5, label_threshold = 10000, distance_from_diagonal_threshold = 0.5, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, pt_size = 1.3, ...) {    

    # Set colors for the plot
    names(my_colors) <- c("DOWN", "UP", "NO")

    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }
    # If genes_to_plot is provided, plot only those genes as described in the first paragraph
    
    ## prepare for visualization
    results_scatter <- results |>  
        drop_na(pvalue) |>
        mutate(
            log10_pval = log10(padj+10^-300)*-1,
            distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2))) |>            
        mutate(
            genes_to_label_first = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                (distance_from_diagonal > distance_from_diagonal_threshold) &
                ((!!sym(paste0('Avg_', group2)) > 100) | (!!sym(paste0('Avg_', group1)) > 100)),
                genes, NA
            ),
            genes_to_label_second = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                ((!!sym(paste0('Avg_', group2)) > label_threshold) | (!!sym(paste0('Avg_', group1)) > label_threshold)) &
                is.na(genes_to_label_first),
                genes, NA
            ),
            genes_to_label = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold) &
                is.na(genes_to_label_first) &
                is.na(genes_to_label_second),
                genes, NA
            ),
            diffexpressed = case_when(
                log2FoldChange >= FC_threshold & padj < p_value_threshold ~ 'UP',
                log2FoldChange <= -1 * FC_threshold & padj < p_value_threshold ~ "DOWN",
                TRUE ~ 'NO'
            )
        ) |>
        mutate(
            diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))
        ) |>
        arrange(diffexpressed)

      # Replace genes to label with provided gene list if applicable
      if (!is.null(genes_to_plot)) {
        results_scatter <- results_scatter %>%
            mutate(
                genes_to_label_first = ifelse(
                    genes %in% genes_to_plot,
                    genes, NA
                ),
                genes_to_label_second = NA,
                genes_to_label = NA)
              }
    
    # Scatterplot
    limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
    limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
    mylims <- max(limx, limy)*6
       
    scatter_plot <- results_scatter |> 
        ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
            geom_point(size=pt_size, stroke = 0) +
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = max_overlaps,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -0.4, 1.25),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 1.25, -0.4),
                fontface = 'italic',
                aes(label = genes_to_label_first,segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -0.4, 1.25),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 1.25, -0.4),
                fontface = 'italic',
                aes(label = genes_to_label_second, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -0.4, 1.25),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 1.25, -0.4),
                fontface = 'italic',
                aes(label = genes_to_label, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste0(test_type, ' DEGs in ', str_replace(cluster,pattern = '_',replace = ' ') ),
                    x=paste0(axis_test, ' in ', group1),
                    y=paste0(axis_test, ' in ',  group2))+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13))+
       scale_x_log10(limits =  c(0.5, mylims), expand = expansion(mult = c(0.01, 0.1)))+
       scale_y_log10(limits =  c(0.5, mylims), expand = expansion(mult = c(0.01, 0.1)))

    ggsave(plot = scatter_plot, filename = paste0(test_type,'_scatter_DEG_in_', cluster, '.pdf'), path = local_figures_path)
    print(scatter_plot)
    return(scatter_plot)    
}

volcano_plot <- function (results, group1, group2, cluster, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size = 5, my_colors = c('green4', 'darkorchid4', 'gray'), test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, pt_size = 1.5, ...) {
    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }

    # Set colors for the plot
    names(my_colors) <- c("DOWN", "UP", "NO")

    nudge_x <- 0
    nudge_y <- 0

    ## prepare for visualization
    results_volcano <- results |>  
        drop_na(pvalue) |>
        mutate(
            log10_pval = log10(padj+10^-300)*-1,
            distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2))) |>            
        mutate(
            genes_to_label_UP = ifelse(
                (log2FoldChange >= FC_threshold) &
                (padj < p_value_threshold),
                genes, NA
            ),
            genes_to_label_DOWN = ifelse(
                (log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold),
                genes, NA
            ),
            diffexpressed = case_when(
                log2FoldChange >= FC_threshold & padj < p_value_threshold ~ 'UP',
                log2FoldChange <= -1*FC_threshold & padj < p_value_threshold ~ "DOWN",
                TRUE ~ 'NO'
            )
        ) |>
        mutate(
            diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))
            
        ) |>
        arrange(diffexpressed)

      # Replace genes to label with provided gene list if applicable
      if (!is.null(genes_to_plot)) {
        results_volcano <- results_volcano %>%
            mutate(
                genes_to_label_UP = ifelse((
                    genes %in% genes_to_plot) & (log2FoldChange >= FC_threshold),
                    genes, NA
                ),
                genes_to_label_DOWN = ifelse((
                    genes %in% genes_to_plot) & (log2FoldChange <= -1 * FC_threshold),
                    genes, NA
                )
            )  
            nudge_x <- 3      
            nudge_y <- 3
      }

    # Remove non-significant genes that would bias the plot visualization
    initial_number_of_genes <- nrow(results_volcano)
    max_FC_up_significant <- results_volcano %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
    min_FC_up_significant <- results_volcano %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }  
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3    
    }
    results_volcano <- results_volcano %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant))) |>
        arrange(diffexpressed)
    final_number_of_genes <- nrow(results_volcano)
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

    volcano_plot <- results_volcano |> 
        arrange(desc(padj)) |>
        ggplot(aes(x=log2FoldChange, y=log10_pval,  col=diffexpressed)) +
        geom_point(size=pt_size, stroke = 0) +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            nudge_x = nudge_x,
            nudge_y = nudge_y,
            aes(label = genes_to_label_UP,segment.size=0.5, segment.alpha=0.8, segment.curvature=0),
            fontface = 'italic') +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            nudge_x = -1*nudge_x,
            nudge_y = nudge_y,
            fontface = 'italic',
            aes(label = genes_to_label_DOWN, segment.size=0.5, segment.alpha=0.8, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        geom_vline(xintercept= FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste0(test_type, ' DEGs in ', str_replace(cluster,pattern = '_',replace = ' ') ),
                    x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                    y= '-Log10 Adj. p-value')+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13),
            )         +
            scale_y_continuous(n.breaks = 8, expand = expansion(mult = c(0.01, 0.1))) +
            scale_x_continuous(n.breaks = 8)
    ggsave(plot = volcano_plot, filename = paste0(test_type, '_volcano_DEG_in_', cluster, '.pdf'), path = local_figures_path)
    print(volcano_plot)
    return(volcano_plot)
}


    
# Core pseudobulk differential expression analysis
pseudobulk_de <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters', 
                          path = './', FC_threshold = 0.3, p_value_threshold = 0.05, 
                          expression_threshold_for_gene_list = 20, minimum_cell_number = 10, 
                          genes_to_exclude = c(), ...) {

    
    requireNamespace('DESeq2', quietly = TRUE) || stop('DESeq2 package needed for this function to work. Please install it.', call. = FALSE)
    
    # Subset seurat object
    scRNAseq <- subset(scRNAseq, subset = (str_detect(!!as.name(comparison), group1) | str_detect(!!as.name(comparison), group2)))
    
    # Set Paths
    gene_lists_path <- here(path, 'gene_lists')
    dir.create(gene_lists_path, showWarnings = FALSE, recursive = TRUE)
    
    print(paste('Cluster', cluster))
    
    group1 <- fixed(group1)
    group2 <- fixed(group2)
    
    Idents(scRNAseq) <- comparison
    
    print('number of cells in group 1')
    print(scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow())
    print('number of cells in group 2')
    print(scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow())
    
    # Check there are enough cells
    n_group1 <- scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group1)) |> nrow()
    n_group2 <- scRNAseq@meta.data |> filter(str_detect(!!as.name(comparison), group2)) |> nrow()
    
    if (n_group1 < minimum_cell_number | n_group2 < minimum_cell_number) {
        return(list(
            all_count = 'Not enough cells',
            UP_count = 'Not enough cells',
            DOWN_count = 'Not enough cells',
            results = NULL
        ))
    }
    
    # Aggregate counts
    counts <- AggregateExpression(scRNAseq, group.by = c(comparison),
                                 assays = 'RNA',
                                 slot = 'counts',
                                 return.seurat = FALSE)
    
    counts <- counts$RNA |> 
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        dplyr::select(-any_of(genes_to_exclude)) |>
        column_to_rownames('genes')
    
    # Generate sample level metadata
    colData <- data.frame(samples = colnames(counts)) |>
        mutate(condition = ifelse(grepl(group1, samples), group1, group2))
    
    # Filter
    counts <- counts |> 
        mutate(row_sums = rowSums(counts)) |> 
        filter(row_sums >= 10) |> 
        dplyr::select(-row_sums)
    
    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))
    
    # Check for sufficient replicates
    if ((length(unique(colData$condition)) != 2) | 
        (nrow(colData |> filter(condition == group1)) < 2) | 
        (nrow(colData |> filter(condition == group2)) < 2)) {
        return(list(
            all_count = 'Not enough biological replicates per group',
            UP_count = 'Not enough biological replicates per group',
            DOWN_count = 'Not enough biological replicates per group',
            results = NULL
        ))
    }
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = colData,
                                  design = ~condition)
    dds$condition <- factor(dds$condition, levels = c(group1, group2))
    
    # DESeq2 QC
    rld <- rlog(dds, blind = TRUE)
    local_figures_path <- here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)
    
    DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition')
    ggsave(filename = paste0('Pseudobulk_PCA_', cluster, '.pdf'), path = local_figures_path)
    
    PCA_table <- DESeq2::plotPCA(rld, ntop = 500, intgroup = 'condition', returnData = TRUE)
    write.csv(PCA_table, file = here(path, paste('PCA_pseudobulk', cluster, group2, 'vs', group1, '.csv', sep = '_')))
    
    # Run DESeq2
    dds <- DESeq(dds)
    resultsNames(dds)
    
    # Generate results object
    results <- results(dds) |> as.data.frame()
    
    # Get Normalized Counts
    normalized_counts <- counts(dds, normalized = TRUE)
    normalized_counts <- normalized_counts |>
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        rowwise() |>
        mutate(
            !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
            !!paste0('Avg_', group1) := mean(c_across(contains(group1)))
        ) |>
        ungroup()
    
    # Add gene annotations
    annotations <- read.csv(here('scripts', 'annotations.csv'))
    results <- results |>
        rownames_to_column('genes') |>
        left_join(y = unique(annotations[, c('gene_name', 'description')]),
                 by = c('genes' = 'gene_name')) |>
        left_join(y = normalized_counts, by = c('genes' = 'genes'))
    
    # Filter results
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) |
             (!!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list)) &
            (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold)
    ) |>
        arrange(padj)
    
    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold)
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1 * FC_threshold)
    
    # Write results to CSV files
    write.csv(results |> arrange(padj), 
             file = here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), 
             file = here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), 
             file = here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep = '_')))
    
    # Return counts and results
    list(
        all_count = nrow(results_filtered),
        UP_count = nrow(results_filtered_UP),
        DOWN_count = nrow(results_filtered_DOWN),
        results = results
    )
}



# Wrapper function for complete pseudobulk analysis with plots
pseudobulk <- function(scRNAseq, comparison, group1, group2, cluster = 'all_clusters', 
                       path = './', FC_threshold = 0.3, p_value_threshold = 0.05, 
                       expression_threshold_for_gene_list = 20, minimum_cell_number = 10, 
                       run_pathway_enrichment = NULL, genes_to_exclude = c(), ...) {
    
    local_figures_path <- here(path, 'figures')
    dir.create(local_figures_path, showWarnings = FALSE, recursive = TRUE)

  group1 <- fixed(group1)
  group2 <- fixed(group2)
    
    # Run core DE analysis
    de_results <- pseudobulk_de(
        scRNAseq = scRNAseq,
        comparison = comparison,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        expression_threshold_for_gene_list = expression_threshold_for_gene_list,
        minimum_cell_number = minimum_cell_number,
        genes_to_exclude = genes_to_exclude,
        ...
    )
    
    # Check if analysis was successful
    if (is.null(de_results$results)) {
        return(list(
            all_count = de_results$all_count,
            UP_count = de_results$UP_count,
            DOWN_count = de_results$DOWN_count
        ))
    }
    
    # Generate scatter and volcano plots
    scatterplot_output <- scatterplot(
        results = de_results$results,
        group1 = group1,
        group2 = group2,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        cluster = cluster,
        test_type = 'Pseudobulk',
        ...
    )
    
    volcanoplot_output <- volcano_plot(
        results = de_results$results,
        group1 = group1,
        group2 = group2,
        cluster = cluster,
        local_figures_path = local_figures_path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        test_type = 'Pseudobulk',
        ...
    )
    
    # Overrepresentation analysis
    run_DEG_functional_analysis(
        results = de_results$results,
        method = run_pathway_enrichment,
        grouping_var = cluster,
        path = path,
        FC_threshold = FC_threshold,
        p_value_threshold = p_value_threshold,
        group1 = group1,
        group2 = group2,
        ...
    )
    
    return(list(
        all_count = de_results$all_count,
        UP_count = de_results$UP_count,
        DOWN_count = de_results$DOWN_count,
        results = de_results$results,
        scatterplot = scatterplot_output,
        volcanoplot = volcanoplot_output
    ))
}

# Wilcox DE analysis

DEG_FindMarkers_RNA_assay <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 30, run_pathway_enrichment = TRUE) {

    # Set colors for the plot
    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")
    
    # Set Paths
    gene_lists_path <- here(path, 'gene_lists/')
    local_figures_path <- here(path, 'figures/')
    dir.create(gene_lists_path)
    dir.create(local_figures_path)
    print(paste('Cluster',cluster))

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    Idents(scRNAseq) <- comparison

    print('number of cells in group 1')
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group1 )) |> nrow())
    print('number of cells in group 2')
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow())
    
    #Check there are enough cells
    if ((scRNAseq@meta.data |> 
            filter(str_detect( !!as.name(comparison) , group1 )) |> 
            nrow()  < minimum_cell_number) | 
            (scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow()  < minimum_cell_number)) {
        DEG_count <- 'Not enough cells'
        DEG_UP_count <- 'Not enough cells'
        DEG_DOWN_count <- 'Not enough cells'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }

    results <- FindMarkers(object = scRNAseq, ident.1 = group1, ident.2 = group2, assay = 'RNA', slot = 'data', test.use = 'wilcox')

    scRNAseq_CPM <- scRNAseq |> AggregateExpression(group.by=c(comparison),
                                        assays = 'RNA',
                                        return.seurat=TRUE,
                                        normalization.method='RC',
                                        scale.factor = 1e6)

    counts_CPM <- scRNAseq_CPM |> GetAssayData(assay = 'RNA', layer = 'data') |>
                            as.data.frame() |>
                            rownames_to_column(var = 'gene') |>
                            mutate(
                                !!paste0('Avg_', group2) := !!as.name(group2),
                                !!paste0('Avg_', group1) := !!as.name(group1),
                            ) 
 

    #Add gene annotations:
    annotations <- read.csv(here('scripts', 'annotations.csv'))    
    results <- results |>
                    rownames_to_column('genes') |>
                    rename(
                                log2FoldChange = avg_log2FC,
                                padj = p_val_adj,
                                pvalue = p_val
                            ) |>
                            mutate(log2FoldChange = log2FoldChange*-1) |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = counts_CPM, by = c('genes' = 'gene'))
    results_filtered <- filter(results, padj < p_value_threshold & (log2FoldChange >= FC_threshold | log2FoldChange <= -1*FC_threshold)) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >=  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <=  -1*FC_threshold)

    # Write results to CSV files
    write.csv(results_filtered |> arrange(padj), file=here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    results_scatter <- scatterplot(results = results, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = max_overlaps, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Wilcox')
    volcano_plot(results_scatter, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = max_overlaps, label_size = label_size, label_threshold = label_threshold, test_type ='Wilcox')

    ########## Overrepresentation analysis ##########
    if (run_pathway_enrichment) {
        Metascape_functional_analysis(results,  grouping_var = cluster, , group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
        if (!is.null(pathways_of_interest)) {
            pathways_of_interest_analysis(results = results, pathways_of_interest = pathways_of_interest,  cluster = cluster, path = path, group1 = group1, group2 = group2, comparison = comparison)
        }
    }

    ########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]                    
            print(genes_to_plot)
            # Generate scatter and volcano plots
            results_scatter <- scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results = results, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Wilcox')
            volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results_scatter = results_scatter, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, test_type = 'Wilcox')

        }
    }    

    return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
}

DEG_FindMarkers_SCT_assay <- function (scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 30, run_pathway_enrichment = TRUE) {

    # Set colors for the plot
    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")
    
    # Set Paths
    gene_lists_path <- here(path, 'gene_lists/')
    local_figures_path <- here(path, 'figures/')
    dir.create(gene_lists_path)
    dir.create(local_figures_path)
    print(paste('Cluster',cluster))

    group1 <- fixed(group1)
    group2 <- fixed(group2)

    Idents(scRNAseq) <- comparison

    print('number of cells in group 1')
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group1 )) |> nrow())
    print('number of cells in group 2')
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow())
    
    #Check there are enough cells
    if ((scRNAseq@meta.data |> 
            filter(str_detect( !!as.name(comparison) , group1 )) |> 
            nrow()  < minimum_cell_number) | 
            (scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow()  < minimum_cell_number)) {
        DEG_count <- 'Not enough cells'
        DEG_UP_count <- 'Not enough cells'
        DEG_DOWN_count <- 'Not enough cells'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }

    results <- FindMarkers(object = scRNAseq, ident.1 = group1, ident.2 = group2, assay = 'SCT', slot = 'data', test.use = 'wilcox', recorrect_umi = !is_integrated_subset)

    scRNAseq_CPM <- scRNAseq |> AggregateExpression(group.by=c(comparison),
                                        assays = 'SCT',
                                        slot = 'counts')
 
    # Calculating CPMs    
    counts_CPM <- scRNAseq_CPM$SCT |> 
                            as.data.frame() |>
                            mutate(across(where(is.numeric), ~ .x / sum(.x) * 1e6)) |> 
                            rownames_to_column(var = 'gene') |>
                            mutate(  
                                !!paste0('Avg_', group2) := !!as.name(group2),# it
                                !!paste0('Avg_', group1) := !!as.name(group1)
                            )


    #Add gene annotations:
    annotations <- read.csv(here('scripts', 'annotations.csv'))    
    results <- results |>
                    rownames_to_column('genes') |>
                    rename(
                                log2FoldChange = avg_log2FC,
                                padj = p_val_adj,
                                pvalue = p_val
                            ) |>
                            mutate(log2FoldChange = log2FoldChange*-1) |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = counts_CPM, by = c('genes' = 'gene'))
    results_filtered <- filter(results, padj < p_value_threshold & (log2FoldChange >= FC_threshold | log2FoldChange <= -1*FC_threshold)) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >= FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1*FC_threshold)

    # Write results to CSV files
    write.csv(results_filtered |> arrange(padj), file=here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'Wilcox', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    results_scatter <- scatterplot(results = results, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = max_overlaps, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Wilcox')
    volcano_plot(results_scatter, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = max_overlaps, label_size = label_size, label_threshold = label_threshold, test_type ='Wilcox')

    ########## Overrepresentation analysis ##########
    if (run_pathway_enrichment) {
        Metascape_functional_analysis(results,  grouping_var = cluster, group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
        if (!is.null(pathways_of_interest)) {
            pathways_of_interest_analysis(results = results, pathways_of_interest = pathways_of_interest,  cluster = cluster, path = path, group1 = group1, group2 = group2, comparison = comparison)
        }
    }

    ########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]                    
            print(genes_to_plot)
            # Generate scatter and volcano plots
            results_scatter <- scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results = results, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Wilcox')
            volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results_scatter = results_scatter, group1 = group1, group2 = group2, cluster = cluster, my_colors = my_colors, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, test_type = 'Wilcox')

        }
    }    

    return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
}


# Bulk functions

bulk_analysis <- function (counts_table, comparison = 'Groups', group1, group2, cluster='', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.4, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, minimum_cell_number = 10, run_pathway_enrichment = FALSE, ...) {

    # Set colors for the plot
    # names(my_colors) <- c("DOWN", "UP", "NO")
    
    # Set Paths
    gene_lists_path <- here(path, 'gene_lists/')
    local_figures_path <- here(path, 'figures/')
    dir.create(gene_lists_path)
    dir.create(local_figures_path)
    print(paste('Cluster',cluster))

    group1 <- fixed(group1)
    group2 <- fixed(group2)
  
    counts <- tibble(counts_table) |> column_to_rownames('genes')

    # Run DE Analysis
    #Generate sample level metadata
    colData <- data.frame(samples=colnames(counts)) |>
                mutate(condition = ifelse(grepl(group1, samples), group1, group2))
    
    ## Filter
    counts <- counts |> mutate(row_sums=rowSums(counts)) |> filter(row_sums >= 10) |> dplyr::select(-row_sums)
    
    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))

    #Perform DESeq2
    if ((length(unique(colData$condition)) != 2 ) | (nrow(colData |> filter(condition == group1)) < 2) | (nrow(colData |> filter(condition == group2)) < 2)) {
        DEG_count <- 'Not enough biological replicates per group'
        DEG_UP_count <- 'Not enough biological replicates per group'
        DEG_DOWN_count <- 'Not enough biological replicates per group'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }
    
    #Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = counts,
        colData = colData,
        design = ~condition)
    dds$condition <- factor(dds$condition, levels = c(group1, group2))    

    ## DESeq2 QC
    rld <- rlog(dds, blind=TRUE) #rlog normalization
    DESeq2::plotPCA(rld, ntop=500, intgroup='condition') #PCA
    ggsave(filename=paste0('Bulk_PCA_', cluster, '.pdf'), path=local_figures_path) 
    PCA_table <- DESeq2::plotPCA(rld, ntop=500, intgroup='condition', returnData = T) #PCA table
    write.csv(PCA_table, file=here(path, paste( 'PCA_bulk', cluster, group2, 'vs', group1, '.csv', sep='_')))

    #################### Run DESeq2
    dds <- DESeq(dds)

    #Check the coefficients for the scRNAseq
    resultsNames(dds)

    #Generate results object
    results <- results(dds) |> as.data.frame()
    
    #Get Normalized Counts
    normalized_counts <- counts(dds, normalized = T)
    normalized_counts <- normalized_counts |>
        as.data.frame() |>
        rownames_to_column('genes') |>
        as_tibble() |>
        rowwise() |>
        mutate(
            !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
            !!paste0('Avg_', group1) := mean(c_across(contains(group1))),
        ) |>
        ungroup()

    #Add gene annotations:
    annotations <- read.csv(here('scripts', 'annotations.csv'))
    results <- results |>
                    rownames_to_column('genes') |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = normalized_counts, by = c('genes' = 'genes'))

    # Filter results 
    results_filtered <- filter(
        results,
        padj < p_value_threshold &
            ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) |
                (!!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list)
            ) &(log2FoldChange >= FC_threshold |
                log2FoldChange <= -1 * FC_threshold)) |> 
        arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >=   FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <= -1*FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj), file= here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    
    
    results_scatter <- scatterplot(results = results, group1 = group1, group2 = group2, cluster = cluster, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Bulk', ...)
    
    volcano_plot(results = results, group1 = group1, group2 = group2, cluster = cluster, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, test_type ='Bulk', ...)

    ########## Overrepresentation analysis ##########
    if (run_pathway_enrichment) {
        Metascape_functional_analysis(results,  grouping_var = cluster, , group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
        if (!is.null(pathways_of_interest)) {
            pathways_of_interest_analysis(results = results, pathways_of_interest = pathways_of_interest,  cluster = cluster, path = path, group1 = group1, group2 = group2, comparison = comparison)
        }
    }

    ########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]                    
            print(genes_to_plot)
            # Generate scatter and volcano plots
            results_scatter <- scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results = results, group1 = group1, group2 = group2, cluster = cluster, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, distance_from_diagonal_threshold = distance_from_diagonal_threshold, test_type = 'Bulk', ...)
            volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results = results, group1 = group1, group2 = group2, cluster = cluster, local_figures_path = local_figures_path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, max_overlaps = 15, label_size = label_size, label_threshold = label_threshold, test_type ='Bulk', ...)

        }
    }   
    
    return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))    
}


# Cluster Annotation Functions
annotate_seurat_with_SingleR_Eduard <- function(
    seurat,
    local_path,
    database = c("ImmGen"),
    annotation_basis = c("cluster_fine", "cell_coarse", "cell_fine"),
    split_by_groups = TRUE
) {
    # Load required libraries
    require(SingleR)
    require(Seurat)
    require(dplyr)
    require(scCustomize) # for DimPlot_scCustom, if used

    # Select database
    if (database == "ImmGen") {
        ref <- ImmGenData(ensembl = FALSE)
    } else {
        stop("Currently only 'ImmGen' database is supported.")
    }

    DefaultAssay(seurat) <- 'RNA'

    # Annotation logic
    if (annotation_basis == "cluster_fine") {
        predictions <- SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.fine,
            cluster = seurat$seurat_clusters
        )
        row.names <- rownames(predictions)
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels)
        predictions_tbl$cluster <- row.names
        annotations <- seurat@meta.data |>
            left_join(predictions_tbl, by = join_by('seurat_clusters' == 'cluster')) |>
            pull(labels)
        seurat$labels_per_cluster_fine <- annotations
        Idents(seurat) <- 'labels_per_cluster_fine'
        p <- DimPlot(seurat, label = FALSE, label.size = 2.5)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cluster_SingleR_annotations_fine','.pdf'), path = local_path, width = 8, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot(seurat, label = TRUE, label.size = 2.5, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cluster_SingleR_annotations_fine_by_group','.pdf'), path = local_path, width = 10, height = 5)
        }
    } else if (annotation_basis == "cluster_coarse") {
        predictions <- SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.main,
            cluster = seurat$seurat_clusters
        )
        row.names <- rownames(predictions)
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels)
        predictions_tbl$cluster <- row.names
        annotations <- seurat@meta.data |>
            left_join(predictions_tbl, by = join_by('seurat_clusters' == 'cluster')) |>
            pull(labels)
        seurat$labels_per_cluster_coarse <- annotations
        Idents(seurat) <- 'labels_per_cluster_coarse'
        p <- DimPlot(seurat, label = FALSE, label.size = 2.5)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cluster_SingleR_annotations_coarse','.pdf'), path = local_path, width = 8, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot(seurat, label = TRUE, label.size = 2.5, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cluster_SingleR_annotations_coarse_by_group','.pdf'), path = local_path, width = 10, height = 5)
        }
    } else if (annotation_basis == "cell_coarse") {
        predictions <- SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.main
        )
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels) |>
            rename(labels_per_cell_coarse = labels)
        seurat$labels_per_cell_coarse <- predictions_tbl |> pull(labels_per_cell_coarse)
        Idents(seurat) <- 'labels_per_cell_coarse'
        p <- DimPlot_scCustom(seurat, label = FALSE)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cell_SingleR_annotations_coarse','.pdf'), path = local_path, width = 5, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot_scCustom(seurat, label = FALSE, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cell_SingleR_annotations_coarse_by_group','.pdf'), path = local_path, width = 6, height = 5)
        }
    } else if (annotation_basis == "cell_fine") {
        predictions <- SingleR(
            test = as.SingleCellExperiment(seurat),
            assay.type.test = 1,
            ref = ref,
            labels = ref$label.fine
        )
        predictions_tbl <- predictions |>
            as_tibble() |>
            dplyr::select(labels) |>
            rename(labels_per_cell_fine = labels)
        seurat$labels_per_cell_fine <- predictions_tbl |> pull(labels_per_cell_fine)
        Idents(seurat) <- 'labels_per_cell_fine'
        p <- DimPlot_scCustom(seurat, label = FALSE)
        print(p)
        ggsave(plot = p, filename = paste0('UMAP_cell_SingleR_annotations_fine','.pdf'), path = local_path, width = 26, height = 5)
        if (split_by_groups) {
            p1 <- DimPlot_scCustom(seurat, label = FALSE, split.by = 'Groups')
            ggsave(plot = p1, filename = paste0('UMAP_cell_SingleR_annotations_fine_by_group','.pdf'), path = local_path, width = 30, height = 5)
        }
    } else {
        stop("annotation_basis must be one of 'cluster_fine', 'cluster_coarse', 'cell_coarse', or 'cell_fine'.")
    }

    DefaultAssay(seurat) <- 'SCT'
    return(seurat)
}

# Find and save top marker genes per cluster in a Seurat object

top_genes_per_cluster <- function (seurat, n_genes_to_plot = 3, object_annotations = '', tables_path = 'results/tables/', figures_path = 'results/figures/', results_path = 'results/', run_pathway_enrichment = NULL, ...) {
    # Function to find top genes per cluster in a Seurat object and save results
    # Args:
    #   seurat: Seurat object
    #   object_annotations: String to append to output file names
    #   tables_path: Path to save the tables
    #   figures_path: Path to save the figures

    sequential_palette_dotplot <- hcl.colors(n = 20,'YlGn',rev = T)
    
    # Set the identity class for clustering
    Idents(seurat) <- 'seurat_clusters'

    seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    # saveRDS(seurat.markers, file = 'seurat.markers.rds')



    #Add gene annotations:
    annotations <- read.csv(here("scripts/annotations.csv"))
    seurat.markers <- seurat.markers |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('gene' = 'gene_name')) |>
                            mutate(cluster = fct_inseq(cluster))

    #Top10 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 10) -> top10

    #Top25 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 50) -> top50

    #Top100 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = 100) -> top100

    #Topn markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) |>
        slice_head(n = n_genes_to_plot) -> topn
    
    # Save the top markers to files
    write.table(top100,file=here(tables_path, paste0('top100', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)
    # write.table(top25,file=here(path,'top25',object_annotations, ".tsv"), sep="\t",row.names = FALSE)
    write.table(top10,file=here(tables_path, paste0('top10', '_',object_annotations, ".tsv")), sep="\t",row.names = FALSE)

    top100_genes_per_cluster <- top100 %>%
        group_by(cluster) %>%
        summarise(genes = str_flatten_comma(gene))

    write.table(top100_genes_per_cluster,
                file = here(tables_path, paste0('top100_gene_names_per_cluster_', object_annotations, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    gene_list_plot <- topn |> pull(gene)

    gene_list_plot <- gene_list_plot |> unique() |> rev()
    plot1 <- DotPlot_scCustom(seurat,
                    features = gene_list_plot,
                    colors_use=sequential_palette_dotplot,
                    flip_axes = T,
                    dot.scale = 8,
                    dot.min = 0,
                    scale.min = 0,
                    scale.max = 80,
                    x_lab_rotate = T,
                    y_lab_rotate = F) +
        theme(axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            legend.title = element_text(size = 18))
    
    metascape_results <- NULL
    ClusterProfiler_results <- NULL
        # Run pathway enrichment analysis
    if (!is.null(run_pathway_enrichment)) {
        if ('Metascape' %in% run_pathway_enrichment) {
            metascape_results <- Metascape_functional_analysis_cluster_identification(seurat, top100, identities = 'seurat_clusters', path=results_path, object_annotations = object_annotations)                        
        }
        if ('ClusterProfiler' %in% run_pathway_enrichment) {
            ClusterProfiler_results <- GO_functional_analysis_cluster_identification(seurat, seurat.markers, path=results_path, object_annotations = object_annotations, top_gene_number = 100, ...)             
        }

    }

    return(list(plot = plot1, ClusterProfiler_results = ClusterProfiler_results, metascape_results = metascape_results, topn = topn, top100 = top100) )

}

extract_cell_counts <- function(seurat, grouping_var, figures_path, tables_path, object_annotations='') {    
    cell_counts <- FetchData(seurat, vars = c(englue("{{grouping_var}}"), "Samples", "Groups"))
    cell_counts <- arrange(cell_counts, Samples)

    counts <- cell_counts %>% add_count(Samples, name='total_cell_count_by_sample') 

    counts <- counts %>% 
        dplyr::count(  {{grouping_var}},  , Samples, Groups,  total_cell_count_by_sample,name='cluster_count')  |> 
            mutate(frequency_within_sample=cluster_count*100/total_cell_count_by_sample)  |> 
            mutate(Samples = as.character(Samples)) |> 
            arrange(Samples, desc(Samples)) 

    frequency_table <- counts |> 
        arrange(Samples) |>     
        pivot_wider(id_cols = {{grouping_var}},  names_from = 'Samples', values_from = frequency_within_sample)
    write.csv(frequency_table,file=here(tables_path, paste0(englue("frequency per {{grouping_var}} per condition "), object_annotations, ".csv")), row.names=F)

    count_table <- counts |> 
        arrange(Samples) |>     
        pivot_wider(id_cols = {{grouping_var}}, names_from = 'Samples', values_from = cluster_count)
    write.csv(count_table,file = here(tables_path, paste0(englue("counts per {{grouping_var}} per condition "), object_annotations, ".csv")),row.names=F)

    # # Barplot of proportion of cells in each cluster by sample
    # plot1 <- ggplot(seurat@meta.data) +
    #     geom_bar(aes(x=Groups, fill={{grouping_var}}), position=position_fill()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
    # ggsave(plot = plot1, filename = paste0(figures_path, englue("frequency per {{grouping_var}} per group"), object_annotations, ".pdf"))
    # print(plot1)

    # counts <- cell_counts %>% add_count(Groups, name='total_cell_count_by_sample')
    # counts <- counts %>% 
    #     dplyr::count({{grouping_var}}, Groups, total_cell_count_by_sample,name='frequency_within_sample')  |> 
    #         mutate(frequency_within_sample=frequency_within_sample*100/total_cell_count_by_sample)  

    # Barplot of proportion of cells in each cluster by sample
    plot2 <- ggplot(counts, aes(x={{grouping_var}} |> fct_reorder(frequency_within_sample) |> fct_rev(), y = frequency_within_sample, fill=Groups)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), show.legend = TRUE) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), aes(fill = Groups, alpha = 0.5)) +
        stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9) ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_classic() +
        labs(x = 'Cell type', y = 'Frequency (%)', title = englue('Frequency per {{grouping_var}}'))+
        guides(alpha = 'none')+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggsave(plot = plot2, filename =  paste0(englue('frequency_per_{{grouping_var}}_per_sample '), object_annotations, '.pdf'), width = 12, height = 6, path = figures_path)
    print(plot2)

#     # Tidyplot of proportion of cells in each cluster by sample
# #     counts <- mutate(
# #         {{grouping_var}} := fct_reorder({{grouping_var}}frequency_within_sample),
# #         {{grouping_var}} := fct_rev({{grouping_var}}frequency_within_sample)
# #   )
#     plot2 <- tidyplot(counts, x = {{grouping_var}}, y = frequency_within_sample, color = Groups)  |>
#       add_data_points_beeswarm() |>
#       add_mean_bar(alpha = 0.4) |>
#       add_sem_errorbar() |>
#       add_test_asterisks(test = 't_test', p.adjust.method = 'BH')  |>
#       save_plot(filename = here(figures_path, paste0(englue('frequency_per_{{grouping_var}}_per_sample_tidyplot '), object_annotations, '.pdf')))
     
}
calculate_D50 <- function (seurat, cell_grouping_var, replicate_var, replicate_group_var = NULL, results_path, figures_path, tables_path, colors  = NULL) {

    #Extracting TCR data for clusters of interest
    Idents(seurat) <- englue("{{cell_grouping_var}}")
    combined2 <- scRepertoire:::.expression2List(seurat, split.by ='ident')
    combined3 <- scRepertoire:::.expression2List(seurat, split.by ='orig.ident')
    grouping_var_levels <- levels(seurat@meta.data |> pull({{cell_grouping_var}}))
    replicate_var_levels <- levels(seurat@meta.data |> pull({{replicate_var}}))

    #Initiating results data frame
    results <- as.data.frame(matrix(nrow = 0,ncol = length(grouping_var_levels)))
    # colnames(results)
    rnames <- c()

    #Calculate D50
    for (HTO in replicate_var_levels) {
        
        result <- c()

        cell_type_HTO_data <- combined3[[1]]|>
                                        filter({{replicate_var}} == HTO) |>
                                        add_count(CTaa, sort=TRUE)
            #Calculating D50
            if (nrow(cell_type_HTO_data) < 20) {
                D50 <- NA 
            } else {
                L50 <- floor(nrow(cell_type_HTO_data)/2)
                number_unique_50 <- cell_type_HTO_data[1:L50,] %>% summarise(n_distinct(CTaa)) %>% as.numeric()
                number_unique_total <- cell_type_HTO_data[] %>% summarise(n_distinct(CTaa)) %>% as.numeric()
                D50 <- number_unique_50/number_unique_total
            }
            result <- c(result, D50)


        for (cell_type in grouping_var_levels) {
            #Extracting data for cell type and HTO
            cell_type_HTO_data <- combined2[[cell_type]]|>
                                        filter({{replicate_var}} == HTO) |>
                                        add_count(CTaa, sort=TRUE)

            #Calculating D50
            if (nrow(cell_type_HTO_data) < 20) {
                D50 <- NA 
            } else {
                L50 <- floor(nrow(cell_type_HTO_data)/2)
                number_unique_50 <- cell_type_HTO_data[1:L50,] %>% summarise(n_distinct(CTaa)) %>% as.numeric()
                number_unique_total <- cell_type_HTO_data[] %>% summarise(n_distinct(CTaa)) %>% as.numeric()
                D50 <- number_unique_50/number_unique_total
            }
            result <- c(result, D50)
        }
        results <- rbind(results, result)
        rnames <- c(rnames, HTO)
        # print(HTO)

    }
    colnames(results) <- paste0(c('All', grouping_var_levels), '_D50')
    results <- results %>% mutate({{replicate_var}} := rnames) |> arrange({{replicate_var}}) |> relocate({{replicate_var}})

    write.csv(results, file = englue("{tables_path}/D50_per_{{cell_grouping_var}}.csv"), row.names=FALSE)
    print(head(results))

    # Convert results to long format for plotting
    results_long <- results %>%
        pivot_longer(-{{replicate_var}}, names_to = englue("{{cell_grouping_var}}"), values_to = "D50") |>
        mutate({{cell_grouping_var}} := fct_inorder({{cell_grouping_var}}))


    
# Add replicate grouping variable information
    if (!(englue('{{replicate_var}}') == englue('{{replicate_group_var}}'))) {
        replicate_grouping_var_info <- seurat@meta.data |> 
            dplyr::select({{replicate_var}}, {{replicate_group_var}})       
        replicate_grouping_var_info <- replicate_grouping_var_info[!duplicated(replicate_grouping_var_info), ]
        results_long <- results_long |> 
            left_join(replicate_grouping_var_info, by = join_by({{replicate_var}} == {{replicate_var}}))
    }


    # Remove NA values for plotting
    results_long <- results_long %>% filter(!is.na(D50))

    # Plot: x axis is {{replicate_var}}, show only dots (no columns)
    plot1 <- ggplot(results_long, aes(x = {{replicate_group_var}}, y = D50, color = {{replicate_group_var}})) +
        geom_quasirandom(size = 1.5) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.7), width = 0.5, aes(fill = {{replicate_group_var}}, alpha = 0.5)) +
        stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2, position = position_dodge(width = 0.7)) +
        scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = c(0, 0.05))) +
        theme_classic() +
        labs(x = englue("{{replicate_group_var}}"), y = "D50 Diversity", title = englue("D50 per {{replicate_group_var}} by {{cell_grouping_var}}")) +
        guides(alpha = 'none')+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
            text = element_text(size = 14),
            axis.line = element_line(colour = "black")) 
    if (!is.null(colors)) {
        plot1 <- plot1 + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    }
    ggsave(plot = plot1, filename = paste0(englue('D50_per_{{cell_grouping_var}}.pdf')), width = length(levels(results_long |> pull({{cell_grouping_var}})))+2, height = 6, path = figures_path)
    print(plot1)
    return(plot1)
}

run_cnmf_results <- function (
    seurat,                    # Seurat object
    data_dir,                  # path to directory containing cNMF run folder (string)
    runname,                   # cNMF run name (string)
    k_used = 6,
    local_density_threshold = 0.1,
    do_cnmf = TRUE,            # whether to run system cnmf command
    sequential_palette = NULL, # palette for FeaturePlot_scCustom
    sequential_palette_dotplot = NULL, # palette for DotPlot_scCustom
    color_palette = NULL, # heatmap color palette
    output_path = here::here('results', 'cNMF', runname, 'outputs'),
    object_annotations = '',
    top_n = 100,                # top genes per program to extract
    topn_plot = 6,             # top genes used in dotplot
    run_pathway_enrichment = pathway_enrichment,       # run Metascape_overrepresentation_analysis
    ...
) {

    requireNamespace("dplyr")
    requireNamespace("tidyr")
    requireNamespace("tibble")
    requireNamespace("here")
    requireNamespace("ggplot2")
    requireNamespace("Seurat")  

    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

    # 1) optionally run cNMF consensus
    if (isTRUE(do_cnmf)) {
        cmd <- paste("cnmf consensus --output-dir", data_dir,
                                 "--name", runname,
                                 "--components", k_used,
                                 "--local-density-threshold", local_density_threshold,
                                 "--show-clustering", sep = " ")
        system(cmd)
    }

    dt_str <- gsub("\\.", "_", as.character(local_density_threshold))
    # 2) build file paths (match your original filenames)
    usage_file <- here(data_dir, runname, paste0(runname, ".usages", ".k_", k_used, ".dt_",  dt_str , ".consensus", ".txt"))
    spectra_score_file <- here(data_dir, runname, paste0(runname, ".gene_spectra_score", ".k_", k_used, ".dt_",  dt_str , ".txt"))
    spectra_tpm_file <- here(data_dir, runname, paste0(runname, ".gene_spectra_tpm", ".k_", k_used, ".dt_",  dt_str , ".txt"))

    # 3) read tables
    usage <- utils::read.table(usage_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
    spectra_score <- utils::read.table(spectra_score_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
    spectra_tpm <- utils::read.table(spectra_tpm_file, sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)                                                

    # 4) normalize usages (per-cell sums to 1)
    usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))

    # 5) attach normalized usages to seurat meta.data (preserve existing meta columns except X*)
    barcodes <- Seurat::Cells(seurat)
    seurat@meta.data <- seurat@meta.data %>%
        dplyr::select(-dplyr::starts_with('X')) %>%
        dplyr::mutate(barcode = barcodes) %>%
        dplyr::left_join(tibble::rownames_to_column(usage_norm, 'barcode'), by = 'barcode') %>%
        tibble::column_to_rownames('barcode')

    # 4) normalize usages (per-cell sums to 1)
    usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
    colnames(usage_norm) <- paste0('cNMF_program_', seq_len(ncol(usage_norm)))


    barcodes <- Cells(seurat)

    seurat@meta.data <- seurat@meta.data |>
        select(-starts_with('cNMF_program_')) |>
        # select(-barcode) |>
        mutate(barcode = barcodes) |>
        left_join(usage_norm |> rownames_to_column('barcode'), by  = 'barcode') |>
        column_to_rownames('barcode')
     

    # 5) Feature plot for programs (relies on FeaturePlot_scCustom existing in environment)
    if (!is.null(sequential_palette)) {
        p1 <- FeaturePlot_scCustom(seurat, features = colnames(usage_norm), colors_use = sequential_palette, num_columns = 3)
    } else {
        p1 <- FeaturePlot_scCustom(seurat, features = colnames(usage_norm), num_columns = 3)
    }
    print(p1)
    ggsave(plot = p1, filename = paste0('FeaturePlot_cNMF_k_', k_used, '_', object_annotations, '.pdf'), path = output_path, width = 15, height = 10)


    signature_violin_plot <- function (signature) {
        plot1 <- VlnPlot(seurat, features = paste0(signature, ''), group.by = 'Groups', pt.size = 0) + labs(title = signature) + NoLegend()
        # print(plot1)
        ggsave(plot = plot1, filename = paste0('VlnPlot_', k_used, signature, '_by_group_', object_annotations, '.pdf'), path = output_path, width = 5, height = 4)
        return(plot1)
    }
    plot_list <- map(colnames(usage_norm), signature_violin_plot)
    plots <- wrap_plots(plot_list)
    print(plots)

    #### Interpreting Top genes per gene expression program
    
    # 8) Extracting top_n genes per program
    top_colnames <- spectra_score |>
        tibble::rownames_to_column("program") |>
        tidyr::pivot_longer(-program, names_to = "gene", values_to = "score") |>
        mutate(program = fct_inseq(program)) |>
        arrange(program) |>
        dplyr::group_by(program) |>
        filter(score > 0) |>
        dplyr::slice_max(order_by = score, n = top_n, with_ties = FALSE) |>
        dplyr::arrange(desc(score), .by_group = TRUE) |> 
        dplyr::select(program, gene) |>
        dplyr::mutate(rank = dplyr::row_number()) |>
        ungroup() |>
        tidyr::pivot_wider(names_from = program, values_from = gene) |>
        dplyr::select(-rank)

    # spectra_score |>
    #     tibble::rownames_to_column("program") |>
    #     tidyr::pivot_longer(-program, names_to = "gene", values_to = "score") |>
    #     mutate(program = fct_inseq(program)) |>
    #     arrange(program) |>
    #     dplyr::group_by(program) |>
    #     filter(score > 0) |>
    #     select(program, gene) |>
    #     summarize(n = n()) |>
    #     print()


    #Add gene annotations:
    annotations <- read.csv(here("scripts/annotations.csv"))
    top_colnames_long <- top_colnames |> 
        pivot_longer(everything(), names_to = 'gene_expression_program', values_to = 'gene')  |>
        mutate(gene_expression_program = fct_inseq(gene_expression_program)) |>
        arrange(gene_expression_program) |>
        left_join(y= unique(annotations[,c('gene_name', 'description')]), by = c('gene' = 'gene_name'))
    write.table(top_colnames_long,file=here(output_path, paste0('top50_per_gene_program_k_', k_used, ".tsv")), sep="\t",row.names = FALSE)

    top_colnames <- top_colnames |>
        select(levels(top_colnames_long$gene_expression_program))

        # 8.2) Extracting the top genes per program for heatmap plotting
    genes_score <- spectra_score |>
        as.matrix() |> 
        t() |> 
        as.data.frame() |>
        rownames_to_column('gene')  |> 
        pivot_longer(-gene, names_to = 'GEP', values_to = 'loading')  |>
        mutate(GEP = fct_inseq(GEP)) |>
        group_by(GEP) |>
        # Z-score normalize the loadings within each GEP
        mutate(loading_zscore = loading) |>
        arrange(desc(loading_zscore), .by_group = TRUE) #|>
        # slice_max(order_by = loading_zscore, n = topn_plot, with_ties = FALSE)  

    genes_to_plot <- genes_score |>
        group_by(GEP) |> 
        slice_max(order_by = loading_zscore, n = topn_plot, with_ties = FALSE) |> 
        ungroup() |> 
        distinct(gene) |> 
        pull(gene)

    genes_score <- genes_score  |> 
        filter(gene %in% genes_to_plot)  |>
        mutate(gene = factor(gene, levels = genes_to_plot))

    # Calculate symmetric limits for the color scale based on z-scores
    max_abs <- max(abs(genes_score$loading_zscore), na.rm = TRUE)
    heatmap_plot <- ggplot(genes_score, aes(x = gene, y = GEP, fill = loading_zscore)) +
        geom_tile(linewidth = 0) +
        scale_fill_gradientn(
            colors = color_palette, 
            name = "z-score<br>loading",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            axis.text.y = element_text(),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = 'Top genes per\nGene Expression Program')
    print(heatmap_plot)

    # 10)  overrepresentation analysis per program
        local_path_pathway_enrichment <- here(output_path, paste0('pathway_enrichment_k_', k_used))
        dir.create(local_path_pathway_enrichment, recursive = TRUE, showWarnings = FALSE)
        over_representation_results <- NULL
        # Run pathway enrichment analysis
    if (!is.null(run_pathway_enrichment)) {
        if ('Metascape' %in% run_pathway_enrichment) {
            plot_list <- list()
            for (gene_expression_program in colnames(top_colnames)) {
                local_path_2 <- here(local_path_pathway_enrichment, gene_expression_program)
                dir.create(local_path_2, recursive = TRUE, showWarnings = FALSE)
                plot_list[[gene_expression_program]] <- Metascape_overrepresentation_analysis(top_colnames %>% dplyr::pull(gene_expression_program),
                            local_path = local_path_2,
                            group = gene_expression_program,
                            filename = paste0(gene_expression_program, '_'),
                            grouping_var = 'GEP',
                            nterms_to_plot_metascape = 20, ...)
                }
            }
        plots <- wrap_plots(plot_list) &
            theme(plot.title = element_text(size = 9, hjust = 0.5),
                axis.text.y = element_text(size = 6),
                axis.text.x = element_text(size = 6),
                axis.title.x  = element_text(size = 7),
            axis.title.y  = element_text(size = 7))
        print(plots)        
            
        if ('ClusterProfiler' %in% run_pathway_enrichment) {
            all_genes <- Features(seurat[['RNA']]) |>unique()
            gene_list <- top_colnames |> as.list() |>  map(~ .x[!is.na(.x)])
            over_representation_results <- GO_overrepresentation_analysis_multiple_lists(gene_list, all_genes, local_path_pathway_enrichment, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '',  drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = F, run_network = F, network_n_terms = 100, nterms_to_plot = 5, font_size = 8, ...)  
        }
    }

    return(list(seurat = seurat, top_colnames_long = top_colnames_long, top_colnames=top_colnames,ORA_results = over_representation_results, genes_score_table = genes_score, spectra_score = spectra_score) )
    }


# Dispatcher wrappers for pathway analyses (supports 1..3 methods)
# - DEG functional dispatcher run_DEG_functional_analysis
# - "ClusterProfiler" used as GO
# - if method is NULL or length 0 returns invisible(NULL)

run_DEG_functional_analysis <- function(results,
                                        method = NULL,
                                        ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "Metascape", "gProfiler2"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("GO_functional_analysis", mode = "function")) {
               GO_functional_analysis(results = results, ...)
             } else stop("GO_functional_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_functional_analysis", mode = "function")) {
               Metascape_functional_analysis(results = results, ...)
             } else stop("Metascape_functional_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_functional_analysis", mode = "function")) {
               gProfiler2_functional_analysis(results = results, ...)
             } else stop("gProfiler2_functional_analysis not found")
           },
           stop("Unknown DEG functional method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}

run_overrepresentation_analysis <- function(genes,
                                           method = NULL,
                                           ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "gProfiler2", "Metascape"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("GO_overrepresentation_analysis", mode = "function")) {
               GO_overrepresentation_analysis(genes, ...)
             } else stop("GO_overrepresentation_analysis_overrepresentation_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_overrepresentation_analysis", mode = "function")) {
               gProfiler2_overrepresentation_analysis(genes, ...)
             } else stop("gProfiler2_overrepresentation_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_overrepresentation_analysis", mode = "function")) {
               Metascape_overrepresentation_analysis(genes, ...)
             } else stop("Metascape_overrepresentation_analysis not found")
           },
           stop("Unknown overrepresentation method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}


# Heatmap plotting function for pathway genes
plot_pathways_heatmap <- function(genes_to_plot, seurat, pathway_name, color_palette = diverging_palette_2, grouping_var = 'Samples') {
        
    Aggregated_expression <- AggregateExpression(seurat, group.by = grouping_var, return.seurat = T)
    data_to_plot <- Aggregated_expression[['RNA']]$data |> 
        as.data.frame() |> 
        rownames_to_column('gene') |> 
        as_tibble() |> 
        filter(gene %in% genes_to_plot) |>
        mutate(gene = factor(gene, levels = genes_to_plot)) |>
        pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'expression') |>
        mutate(Sample = factor(Sample, levels = colnames(Aggregated_expression[['RNA']]))) |>
        group_by(gene) |>
        mutate(scaled_expression = scale(expression)[,1]) |>
        ungroup() |>
        mutate(Sample = fct_rev(factor(Sample)))

    # Calculate symmetric limits for the color scale
    max_abs <- max(abs(data_to_plot$scaled_expression), na.rm = TRUE)

    plot <- ggplot(data_to_plot, aes(x = gene, y = Sample, fill = scaled_expression)) +
        geom_tile(color = "grey60", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = color_palette, 
            name = "z-score<br>expression",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
            legend.title = element_text(),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = pathway_name)
    
    return(plot)
}

plot_pathways_heatmap2 <- function(genes_to_plot, seurat, pathway_name, color_palette = diverging_palette_2, grouping_var = 'Samples') {
        
    Aggregated_expression <- AggregateExpression(seurat, group.by = grouping_var, return.seurat = T)
    data_to_plot <- Aggregated_expression[['RNA']]$data |> 
        as.data.frame() |> 
        rownames_to_column('gene') |> 
        as_tibble() |> 
        filter(gene %in% genes_to_plot) |>
        mutate(gene = factor(gene, levels = genes_to_plot)) |>
        pivot_longer(cols = -gene, names_to = 'Sample', values_to = 'expression') |>
        mutate(Sample = factor(Sample, levels = colnames(Aggregated_expression[['RNA']]))) |>
        group_by(gene) |>
        mutate(scaled_expression = scale(expression)[,1]) |>
        ungroup() |>
        mutate(Sample = fct_rev(factor(Sample)))
    

    # Calculate symmetric limits for the color scale
    max_abs <- max(abs(data_to_plot$scaled_expression), na.rm = TRUE)

    plot <- ggplot(data_to_plot, aes(x = Sample, y = gene, fill = scaled_expression)) +
        geom_tile(color = "grey60", linewidth = 0.3) +
        scale_fill_gradientn(
            colors = color_palette, 
            name = "z-scored\nexpression",
            limits = c(-max_abs, max_abs)
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = 14)
        ) +
        labs(x = NULL, y = NULL, title = pathway_name)
    
    return(plot)
}

plot_clonal_repertoire <- function(seurat, group_name, figures_path, viridis_option = 'D') {
  # Prepare repertoire data
  repertoire <- seurat@meta.data |>
    filter(Groups == group_name & !is.na(CTaa)) |>
    arrange(desc(clonalFrequency))
  
  # Determine unique clone sizes and generate appropriate number of colors
  repertoire2 <- repertoire |> arrange(cloneSize)
  clone_sizes <- unique(repertoire2$cloneSize[!is.na(repertoire$cloneSize)]) |> as.character()
  n_colors <- length(clone_sizes)
  viridis_colors <- viridis(n = n_colors, option = viridis_option, direction = -1)
  print(clone_sizes)
  
  # Create color mapping
  repertoire <- repertoire |>
    mutate(color = case_when(
      str_equal(cloneSize, clone_sizes[1]) ~ viridis_colors[1],
      str_equal(cloneSize, clone_sizes[min(2, n_colors)]) ~ viridis_colors[min(2, n_colors)],
      str_equal(cloneSize, clone_sizes[min(3, n_colors)]) ~ viridis_colors[min(3, n_colors)],
      str_equal(cloneSize, clone_sizes[min(4, n_colors)]) ~ viridis_colors[min(4, n_colors)],
      str_equal(cloneSize, clone_sizes[min(5, n_colors)]) ~ viridis_colors[n_colors],
      is.na(cloneSize) ~ 'grey90'
    )) |>
    select(clonalFrequency, CTaa, color) |>
    distinct()
  
  # Calculate circle packing layout
  packing <- circleProgressiveLayout(repertoire$clonalFrequency, sizetype = 'area')
  repertoire <- repertoire |>
    mutate(x = packing$x,
           y = packing$y)
  ggplot_data_circles <- circleLayoutVertices(packing, npoints = 50)
  
  # Create plot
  plot <- ggplot(data = ggplot_data_circles) +
    geom_polygon(aes(x, y, group = id, fill = factor(id)), 
                 colour = "black", alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = repertoire$color, labels = repertoire$cloneSize) +
    geom_text(data = repertoire,
              aes(x, y, size = clonalFrequency, label = clonalFrequency)) +
    guides(fill = guide_legend(), size = "none") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = paste0("Clonal repertoire for group: ", group_name)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    coord_equal()
  ggsave(plot = plot, filename = paste0('clonal_repertoire_', group_name, '.pdf'), width = 8, height = 8, path = figures_path)

#print(plot)
return(plot)
}