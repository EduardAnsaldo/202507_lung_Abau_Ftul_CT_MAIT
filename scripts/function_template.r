#Visualization Functions
## Scatter Plot Function
scatterplot <- function (results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 30, label_size, label_threshold, distance_from_diagonal_threshold, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, gene_list_name = NULL) {
    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }
    # If genes_to_plot is provided, plot only those genes as described in the first paragraph
    if (!is.null(genes_to_plot)) {
        results_scatter <- results |>  
            drop_na(pvalue) |>
            mutate(
                log10_pval = log10(padj+10^-90)*-1,
                distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2)),
                genes_to_label_first = ifelse(genes %in% genes_to_plot, genes, NA),
                diffexpressed = case_when(
                    log2FoldChange>=FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot) ~ 'UP', 
                    log2FoldChange<=-1*FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot)  ~ "DOWN",
                    genes %in% genes_to_plot & log2FoldChange >= 0 & padj < p_value_threshold ~ 'INTEREST UP',
                    genes %in% genes_to_plot & log2FoldChange < 0 & padj < p_value_threshold ~ 'INTEREST DOWN',
                    TRUE ~ 'NO')
            ) |>
            mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP', 'INTEREST DOWN', 'INTEREST UP'))) |>
            arrange(diffexpressed)

        my_colors <- c("green4", "darkorchid4", "gray", 'green3', 'darkorchid3')
        names(my_colors) <- c("DOWN", "UP", "NO", 'INTEREST DOWN', 'INTEREST UP')

        limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
        limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*6

        scatter_plot <- results_scatter |> 
            filter(!is.na(genes_to_label_first)) |>
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
                geom_point(size=1.3) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.35,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    max.time = 10,
                    max.iter = 10000000,
                    nudge_x = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', -1, 1),
                    nudge_y = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', 0.75, -0.75),
                    aes(label = genes_to_label_first,segment.size=0.5, segment.alpha=1, segment.curvature=0)) +
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

        plot_filename <- if (!is.null(gene_list_name)) {
            paste0(local_figures_path, test_type, '_scatter_', gene_list_name, '_DEG_in_', cluster, '.pdf')
        } else {
            paste0(local_figures_path, test_type, '_scatter_genes_of_interest_DEG_in_', cluster, '.pdf')
        }
        ggsave(plot = scatter_plot, filename = plot_filename)
        print(scatter_plot)
        return(results_scatter)
    }

    # Otherwise, do the original function as described
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

    ggsave(plot = scatter_plot, filename = paste0(local_figures_path, test_type,'_scatter_DEG_in_', cluster, '.pdf'))
    print(scatter_plot)
    return(results_scatter)    
}
## Volcano Plot Function
volcano_plot <- function (results_scatter, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, test_type = c('Wilcox', 'Pseudobulk', 'Bulk'), genes_to_plot = NULL, gene_list_name = NULL) {
    #Determine test type
    if (test_type == 'Pseudobulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Bulk') {
        axis_test <- 'Average Normalized Counts'
    } else if (test_type == 'Wilcox') {
        axis_test <- 'Average CPMs'
    }

    # If genes_to_plot is provided, plot only those genes as described in the other code blocks
    if (!is.null(genes_to_plot)) {
        results_volcano <- results_scatter |>
            mutate(
                genes_to_label_first = ifelse(genes %in% genes_to_plot, genes, NA),
                diffexpressed = case_when(
                    log2FoldChange>=FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot) ~ 'UP', 
                    log2FoldChange<=-1*FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot)  ~ "DOWN",
                    genes %in% genes_to_plot & log2FoldChange >= 0 & padj < p_value_threshold ~ 'INTEREST UP',
                    genes %in% genes_to_plot & log2FoldChange < 0 & padj < p_value_threshold ~ 'INTEREST DOWN',
                    TRUE ~ 'NO')
            ) |>
            mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP', 'INTEREST DOWN', 'INTEREST UP'))) |>
            arrange(diffexpressed)

        my_colors <- c("green4", "darkorchid4", "gray", 'green3', 'darkorchid3')
        names(my_colors) <- c("DOWN", "UP", "NO", 'INTEREST DOWN', 'INTEREST UP')

        # Filter values which are not significant but with high FC that would bias the plot visualization
        initial_number_of_genes <- nrow(results_volcano)
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

        volcano_plot_2 <- results_volcano |> 
            ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_first, col=diffexpressed)) +
            geom_point(size=1.5) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = Inf,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_volcano$diffexpressed == 'UP' | results_volcano$diffexpressed == 'INTEREST UP', 3 ,-3),
                nudge_y = ifelse(results_volcano$diffexpressed == 'UP', 3, 3),
                aes(segment.size=0.3, segment.alpha=1, segment.curvature=0)) +
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
                scale_y_continuous(n.breaks = 8) +
                scale_x_continuous(n.breaks = 8)
        plot_filename <- if (!is.null(gene_list_name)) {
            paste0(local_figures_path, test_type, '_volcano_', gene_list_name, '_DEG_in_', cluster, '.pdf')
        } else {
            paste0(local_figures_path, test_type, '_volcano_genes_of_interest_DEG_in_', cluster, '.pdf')
        }
        ggsave(plot=volcano_plot_2, filename = plot_filename)
        print(volcano_plot_2)
        return(invisible(NULL))
    }

    # Otherwise, do the original function as described
    initial_number_of_genes <- nrow(results_scatter)
    max_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
    min_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }  
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3    
    }
    results_volcano <- results_scatter %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
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
    ggsave(plot = volcano_plot, filename = paste0(local_figures_path, test_type, '_volcano_DEG_in_', cluster, '.pdf'))
    print(volcano_plot)
}

# Pseudobulk function

pseudobulk <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.4, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 10, run_gProfiler2 = FALSE) {

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

    counts <- counts$RNA |> as.data.frame()

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
    write.csv(PCA_table, file=paste(path, 'PCA_pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_'))

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
    results_filtered_UP <- filter(results_filtered, log2FoldChange >  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <  FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj), file= here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'pseudobulk', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    results_scatter <- scatterplot(results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, distance_from_diagonal_threshold, test_type = 'Pseudobulk')
    volcano_plot(results_scatter, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, test_type ='Pseudobulk')

    ########## Overrepresentation analysis ##########
    if (run_gProfiler2) {
        gProfiler2_functional_analysis(results,  cluster = cluster, , group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
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
            results_scatter <- scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, distance_from_diagonal_threshold, test_type = 'Pseudobulk')
            volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results_scatter, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, test_type = 'Pseudobulk')

        }
    }   
    
    return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))    
}                            

# Wilcox DE analysis

DEG_FindMarkers_RNA_assay <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 30, run_gProfiler2 = FALSE) {

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
    results_filtered <- filter(results, padj < p_value_threshold) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <  FC_threshold)

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
    if (run_gProfiler2) {
        gProfiler2_functional_analysis(results,  cluster = cluster, , group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
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

DEG_FindMarkers_SCT_assay <- function (scRNAseq, comparison, group1, group2, is_integrated_subset = FALSE, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.7, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 30, run_gProfiler2 = FALSE) {

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
    results_filtered <- filter(results, padj < p_value_threshold) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange > FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange < FC_threshold)

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
    if (run_gProfiler2) {
        gProfiler2_functional_analysis(results,  cluster = cluster, group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
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

bulk_analysis <- function (counts_table, comparison = 'Groups', group1, group2, cluster='', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.4, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4'), minimum_cell_number = 10, run_gProfiler2 = FALSE) {

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
    write.csv(PCA_table, file=paste(path, 'PCA_bulk', cluster, group2, 'vs', group1, '.csv', sep='_'))

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
    results_filtered_UP <- filter(results_filtered, log2FoldChange >  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <  FC_threshold)

    # Write results to CSV files
    write.csv(results |> arrange(padj), file= here(gene_lists_path, paste('ALL_GENES_DEG_Analysis', cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_UP |> arrange(desc(log2FoldChange)), file=here(gene_lists_path, paste('DEG_UP_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))
    write.csv(results_filtered_DOWN |> arrange(log2FoldChange), file=here(gene_lists_path, paste('DEG_DOWN_in', group2, cluster, 'bulk', group2, 'vs', group1, '.csv', sep='_')))

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)

    # Generate scatter and volcano plots
    results_scatter <- scatterplot(results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, distance_from_diagonal_threshold, test_type = 'Bulk')
    volcano_plot(results_scatter, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, test_type ='Bulk')

    ########## Overrepresentation analysis ##########
    if (run_gProfiler2) {
        gProfiler2_functional_analysis(results,  cluster = cluster, , group1 = group1, group2 = group2, path= path , FC_threshold = FC_threshold)
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
            results_scatter <- scatterplot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, distance_from_diagonal_threshold, test_type = 'Bulk')
            volcano_plot(genes_to_plot = genes_to_plot, gene_list_name = gene_list, results_scatter, group1, group2, cluster, my_colors, local_figures_path, FC_threshold, p_value_threshold, max_overlaps = 15, label_size, label_threshold, test_type = 'Bulk')

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

top_genes_per_cluster <- function (seurat, object_annotations = '', tables_path = 'results/tables/', figures_path = 'results/figures/', results_path = 'results/', run_pathway_enrichment = FALSE) {
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
                        by = c('gene' = 'gene_name'))

    #Top10 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) |>
        slice_head(n = 10) -> top10

    #Top25 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) |>
        slice_head(n = 25) |>
        ungroup() -> top25

    #Top100 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) |>
        slice_head(n = 100) -> top100

    #Top3 markers
    seurat.markers %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) |>
        slice_head(n = 3) -> top3
    
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

    gene_list_plot <- top3 |> pull(gene)

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
    
    if (run_pathway_enrichment) {
        # Run pathway enrichment analysis
        gProfiler2_functional_analysis_cluster_identification(seurat, top25, identities = 'seurat_clusters', path=results_path, object_annotations = object_annotations) 
    }

    return(plot1)

}

extract_cell_counts <- function(seurat, grouping_var, figures_path, tables_path, object_annotations='') {    
    cell_counts <- FetchData(seurat, vars = c(englue("{{grouping_var}}"), "Samples", "Groups"))
    cell_counts <- arrange(cell_counts, Samples)

    counts <- cell_counts %>% add_count(Samples, name='total_cell_count_by_sample') 

    counts <- counts %>% 
        dplyr::count(  {{grouping_var}},  , Samples, total_cell_count_by_sample,name='cluster_count')  |> 
            mutate(frequency_within_sample=cluster_count*100/total_cell_count_by_sample)  |> 
            mutate(Samples = as.character(Samples)) |> 
            arrange(Samples, desc(Samples)) 

    frequency_table <- counts |> 
        arrange(Samples) |>     
        pivot_wider(id_cols = {{grouping_var}},  names_from = 'Samples', values_from = frequency_within_sample)
    write.csv(frequency_table,file=paste0(tables_path, englue("frequency per {{grouping_var}} per condition "), object_annotations, ".csv"), row.names=F)

    count_table <- counts |> 
        arrange(Samples) |>     
        pivot_wider(id_cols = {{grouping_var}}, names_from = 'Samples', values_from = cluster_count)
    write.csv(count_table,file=paste0(tables_path, englue("counts per {{grouping_var}} per condition "), object_annotations, ".csv"),row.names=F)

    # # Barplot of proportion of cells in each cluster by sample
    # plot1 <- ggplot(seurat@meta.data) +
    #     geom_bar(aes(x=Groups, fill={{grouping_var}}), position=position_fill()) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
    # ggsave(plot = plot1, filename = paste0(figures_path, englue("frequency per {{grouping_var}} per group"), object_annotations, ".pdf"))
    # print(plot1)

    counts <- cell_counts %>% add_count(Groups, name='total_cell_count_by_sample')
    counts <- counts %>% 
        dplyr::count({{grouping_var}}, Groups, total_cell_count_by_sample,name='frequency_within_sample')  |> 
            mutate(frequency_within_sample=frequency_within_sample*100/total_cell_count_by_sample)  

    # Barplot of proportion of cells in each cluster by sample
    plot2 <- ggplot(counts, aes(x={{grouping_var}} |> fct_reorder(frequency_within_sample) |> fct_rev(), y = frequency_within_sample,fill=Groups)) +
        geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), show.legend = TRUE) +
        stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.9), aes(fill = Groups, alpha = 0.5)) +
        stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
            geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9) ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme_classic() +
        labs(x = 'Cell type', y = 'Frequency (%)', title = englue('Frequency per {{grouping_var}}'))+
        guides(alpha = 'none')+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggsave(plot = plot2, filename =  paste0(figures_path, englue('frequency_per_{{grouping_var}}_per_sample '), object_annotations, '.pdf'), width = 12, height = 6)
    print(plot2)
}




