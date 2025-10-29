scatterplot <- function (results, group1, group2, local_figures_path, FC_threshold, p_value_threshold, cluster = 'all_clusters', my_colors = c('green4', 'darkorchid4'), max_overlaps = 15, label_size = 5, label_threshold = 10000, distance_from_diagonal_threshold = 0.5, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, ...) {    

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
            log10_pval = log10(padj+10^-90)*-1,
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
            geom_point(size=1.3
            ) +
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = max_overlaps,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_first,segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_second, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
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
       scale_x_log10(limits =  c(0.5, mylims))+
       scale_y_log10(limits =  c(0.5, mylims))

    ggsave(plot = scatter_plot, filename = paste0(test_type,'_scatter_DEG_in_', cluster, '.pdf'), path = local_figures_path)
    print(scatter_plot)
    return(scatter_plot)    
}

volcano_plot <- function (results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size = 5, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, ...) {
    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }

## prepare for visualization
    results_volcano <- results |>  
        drop_na(pvalue) |>
        mutate(
            log10_pval = log10(padj+10^-90)*-1,
            distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2))) |>            
        mutate(
            genes_to_label_volcano = ifelse(
                (log2FoldChange >= FC_threshold | log2FoldChange <= -1 * FC_threshold) &
                (padj < p_value_threshold),
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
        results_volcano <- results_volcano %>%
            mutate(
                genes_to_label_volcano = ifelse(
                    genes %in% genes_to_plot,
                    genes, NA
                )
            )        
      }

    max_FC_up_significant <- results_volcano %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
    min_FC_up_significant <- results_volcano %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }  
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3    
    }
    results_volcano <- results_volcano %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
    final_number_of_genes <- nrow(results_volcano)
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

    volcano_plot <- results_volcano |> 
        arrange(desc(padj)) |>
        ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_volcano, col=diffexpressed)) +
        geom_point(size=1.5) +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            aes(segment.size=0.5, segment.alpha=0.8, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        geom_vline(xintercept=FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
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
            scale_y_continuous(n.breaks = 8) +
            scale_x_continuous(n.breaks = 8)
    ggsave(plot = volcano_plot, filename = paste0(test_type, '_volcano_DEG_in_', cluster, '.pdf'), path = local_figures_path)
    print(volcano_plot)
    return(volcano_plot)
}
pseudobulk <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, pathways_of_interest = NULL, genes_to_plot = NULL, expression_threshold_for_gene_list = 20, minimum_cell_number = 10, run_pathway_enrichment = NULL, genes_to_exclude = c(), ...) {

    requireNamespace('DESeq2', quietly = TRUE) || stop('DESeq2 package needed for this function to work. Please install it.', call. = FALSE)
 
    # Subset seurat object
    scRNAseq <- subset(scRNAseq, subset = (str_detect(!!as.name(comparison), group1) | str_detect(!!as.name(comparison), group2)))

    # Set colors for the plot
    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")
    # Set Paths
    gene_lists_path <- here(path, 'gene_lists')
    local_figures_path <- here(path, 'figures')
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
    
    # Aggregate counts
    counts <- AggregateExpression(scRNAseq, group.by=c(comparison),
                            assays='RNA',
                            slot='counts',
                            return.seurat=FALSE)

    counts <- counts$RNA |> as.data.frame() |>
                rownames_to_column('genes') |>
                as_tibble() |>
                dplyr::select(-any_of(genes_to_exclude)) |>
                column_to_rownames('genes')

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
    ggsave(filename=paste0('Pseudobulk_PCA_', cluster, '.pdf'), path=local_figures_path) 
    PCA_table <- DESeq2::plotPCA(rld, ntop=500, intgroup='condition', returnData = T) #PCA table
    write.csv(PCA_table, file = here(path, paste('PCA_pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_')))

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
    results_filtered_UP <- filter(results_filtered, log2FoldChange >=  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <=  -1*FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj), file= here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    scatterplot_output <- scatterplot(results, group1, group2, cluster, local_figures_path, FC_threshold, p_value_threshold, test_type = 'Pseudobulk', ...)
    volcanoplot_output <- volcano_plot(results, group1, group2, cluster, local_figures_path, FC_threshold, p_value_threshold, test_type ='Pseudobulk', ...)

    ########## Overrepresentation analysis ##########
    
    run_DEG_functional_analysis(results, method = run_pathway_enrichment, grouping_var = cluster, path= path, FC_threshold = FC_threshold, p_value_threshold = p_value_threshold, group1 = group1, group2 = group2, ...)

    return(list(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count,  results  = results, scatterplot = scatterplot_output, volcanoplot = volcanoplot_output))    
}