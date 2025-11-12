## Functions
### Overrepresentation analysis function GO and MSigDB -- ClusterProfiler

#' GO Overrepresentation Analysis
#'
#' Performs Gene Ontology overrepresentation analysis using clusterProfiler
#'
#' @param significant_genes Character vector of significant gene symbols
#' @param all_genes Character vector of all genes in the universe
#' @param local_path Path to save output files
#' @param ontology GO ontology to use: 'ALL', 'BP', 'MF', or 'CC' (default: 'ALL')
#' @param minGSSize Minimum gene set size (default: 5)
#' @param maxGSSize Maximum gene set size (default: 400)
#' @param filename Prefix for output files (default: '')
#' @param group Group name for plot titles (default: '')
#' @param drop_levels Whether to drop GO levels (default: FALSE)
#' @param levels_to_drop GO levels to drop if drop_levels is TRUE (default: c())
#' @param simplify_function Function to use for simplification (default: min)
#' @param simplify_by Column to use for simplification (default: 'p.adjust')
#' @param simplify_terms Whether to simplify GO terms (default: TRUE)
#' @param run_network Whether to generate network plots (default: FALSE)
#' @param network_n_terms Number of terms to show in network (default: 100)
#' @param nterms_to_plot Number of terms to plot in dotplot (default: 50)
#' @param color_palette Color palette for plots (default: viridis)
#' @param organism OrgDb organism database (default: org.Mm.eg.db for mouse)
#' @param ... Additional arguments passed to other functions
#'
#' @return Invisibly returns NULL (results are saved to files)
#' @export
GO_overrepresentation_analysis <- function (significant_genes, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '', group = '', drop_levels = FALSE, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = TRUE, run_network = FALSE, network_n_terms = 100, nterms_to_plot = 50, color_palette = NULL, organism = org.Mm.eg.db, ...)  {

     # Input validation
     stopifnot("significant_genes must be a character vector" = is.character(significant_genes))
     stopifnot("all_genes must be a character vector" = is.character(all_genes))
     stopifnot("local_path must be a character string" = is.character(local_path) && length(local_path) == 1)
     stopifnot("ontology must be one of 'ALL', 'BP', 'MF', or 'CC'" = ontology %in% c('ALL', 'BP', 'MF', 'CC'))
     stopifnot("significant_genes cannot be empty" = length(significant_genes) > 0)

     if (drop_levels && length(levels_to_drop) == 0) {
          warning("drop_levels is TRUE but levels_to_drop is empty. No levels will be dropped.")
     }

     # Set color palette
     if (is.null(color_palette)) {
          color_scale <- viridis(n = 4, direction = -1)
     } else {
          color_scale <- color_palette
     }
     options(enrichplot.colours = color_scale)
     
     enrichment_results <- tryCatch({
          enrichGO(gene = significant_genes,
                    universe = all_genes,
                    keyType = "SYMBOL",
                    OrgDb = organism,
                    ont = ontology,
                    pAdjustMethod = "BH",
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    qvalueCutoff = 0.25)
     }, error = function(e) {
          stop("enrichGO failed: ", e$message, "\nCheck that gene symbols are valid for the specified organism.")
     })

     if (drop_levels == TRUE) {
            enrichment_results <- dropGO(enrichment_results, level = levels_to_drop)
     }
     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here(local_path, paste0(filename,'GO_OverRepresentation_analysis_results_', ontology, '.csv')))

     if (nrow(enrichment_results_table) > 1) {
       if (simplify_terms == TRUE) {
            ## Add similarity matrix to the termsim slot of enrichment result
          enrichment_results <- enrichplot::pairwise_termsim(enrichment_results, showCategory = dim(enrichment_results)[1])
          enrichment_results_unfiltered <-  enrichment_results
          enrichment_results <- clusterProfiler::simplify(enrichment_results, cutoff=0.7, by=simplify_by, select_fun=simplify_function)        
       }

          write.csv(as_tibble(enrichment_results), here(local_path, paste0(filename,'GO_OverRepresentation_analysis_results_filtered_', ontology, '.csv')))
          p1 <- dotplot(enrichment_results,
               showCategory=nterms_to_plot,
               title = paste0(filename,'GO ORA UP in ', group, ' - ', ontology),
            #    x = 'p.adjust',
               label_format = 60)
          print(p1)
          ggsave(plot = p1, filename = paste0(filename, 'GO overrepresentation_analysis_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## =  Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
          if (run_network) {
               p2 <- try(emapplot(enrichment_results_unfiltered, showCategory = network_n_terms) + ggtitle(paste0(filename, 'GO Overrepresentation analysis ', ontology)))
               print(p2)
               p2 <- try(emapplot(enrichment_results, showCategory = network_n_terms) + ggtitle(paste0(filename, 'GO Overrepresentation analysis filtered terms', ontology)))
               print(p2)
            #    ggsave(plot = p2, filename = paste0(filename, 'GO_overrepresentation_analysis_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)
            
          }
# #           p2 <- try(emapplot(enrichment_results, showCategory = (nrow(enrichment_results_table)-1)) + ggtitle(paste0(filename, 'Overrepresentation analysis ', ontology)))
#         #   p2 <- try(emapplot(enrichment_results, showCategory = (nrow(enrichment_results_table)-1)) + ggtitle(paste0(filename, 'Overrepresentation analysis ', ontology)))
#           p2 <- treeplot(enrichment_results, showCategory = 50) + ggtitle(paste0(filename, 'Overrepresentation analysis ', ontology))
#                print(p2)

            
          # ggsave(plot = p2, filename = paste0(filename, 'GO_overrepresentation_analysis_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)
     }
}

#' GO Gene Set Enrichment Analysis
#'
#' Performs Gene Ontology GSEA using clusterProfiler
#'
#' @param results Data frame with columns 'genes' and 'log2FoldChange'
#' @param local_path Path to save output files
#' @param ontology GO ontology to use: 'ALL', 'BP', 'MF', or 'CC' (default: 'ALL')
#' @param group Group name for plot titles
#' @param color_palette Color palette for plots (default: viridis)
#' @param organism OrgDb organism database (default: org.Mm.eg.db for mouse)
#'
#' @return Invisibly returns NULL (results are saved to files)
#' @export
GO_GSEA_analysis <- function (results, local_path, ontology = 'ALL', group, color_palette = NULL, organism = org.Mm.eg.db) {

     # Input validation
     stopifnot("results must be a data frame" = is.data.frame(results))
     stopifnot("results must contain 'genes' column" = "genes" %in% colnames(results))
     stopifnot("results must contain 'log2FoldChange' column" = "log2FoldChange" %in% colnames(results))
     stopifnot("ontology must be one of 'ALL', 'BP', 'MF', or 'CC'" = ontology %in% c('ALL', 'BP', 'MF', 'CC'))

     # Set color palette
     if (is.null(color_palette)) {
          color_scale <- viridis(n = 4, direction = -1)
     } else {
          color_scale <- color_palette
     }
     options(enrichplot.colours = color_scale)

          #### GSEA ####

     fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
     names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(genes)

     gsea_results <- tryCatch({
          gseGO(geneList     = fold_changes,
               OrgDb        = organism,
               ont          = ontology,
               keyType = "SYMBOL",
               minGSSize    = 5,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE,
            eps = 0)
     }, error = function(e) {
          stop("gseGO failed: ", e$message, "\nCheck that gene symbols are valid and fold changes are properly named.")
     })               

     gsea_results_table <- as_tibble(gsea_results) 
     write.csv(gsea_results_table, here(local_path, paste0('Gene_Set_Enrichment_Analysis_results_', ontology,'.csv')))

     if (nrow(gsea_results_table) > 1) {
          p3 <- dotplot(gsea_results,
               showCategory=50,
               title = paste0('GSEA analysis ', group, ' - ', ontology),
               label_format = 60)
          print(p3)
          ggsave(plot = p3, filename = paste0('GSEA_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## Add similarity matrix to thenes,  termsim slot of enrichment result
          gsea_results <- enrichplot::pairwise_termsim(gsea_results, showCategory = dim(gsea_results)[1])

        #   ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
        # try(emapplot(gsea_results) + ggtitle(paste0('GSEA analysis ', ontology)))
        #   ggsave(paste0('GSEA_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)

# ## Make individual GSEA plots
#           local_path2 <- here(local_path, paste0('Individual_GSEA_plot_', ontology))
#           unlink(local_path2, recursive = T)
#           dir.create(local_path2)

#           for (pathway in head(gsea_results$ID, 20)) {
#                pathway_name <- gsea_results |> as_tibble() |> filter(ID == pathway) |> pull(Description)
#                anno <- gsea_results[pathway, c("NES", "pvalue", "p.adjust")]
#                lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

#                try(p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 1, base_size = 13,title = paste0(pathway, ' ', pathway_name)))

#                try(x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75)
#                try(y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17)

#                try(p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5))
#                try(p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 2, base_size = 13))
#                try(p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 3, base_size = 13))

#                try(cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v'))    
#                try(ggsave(paste0('GSEA ',pathway_name , '.pdf'), path = local_path2, height = 10, width = 8))
#                }
           }
}

GO_functional_analysis <- function (results,  grouping_var, path='./', FC_threshold = 0.3, p_value_threshold = 0.05, group1 = '', group2 = '', run_GSEA = FALSE, ...) {
 
    #  results <- results[which(duplicated(results$genes) == F),]
#    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    #results <- results[which(duplicated(results$entrezid) == F),]
cluster <- grouping_var

####################################### UP ########################################

     local_path <- here(path, paste0('GO_functional_analysis_UP_in_', group2, '_', cluster))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> pull(genes)
     all_genes <- results |> arrange(padj) |> pull(genes)

     ######################## ORA ########################

#      GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'ALL', group = group2)
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group2, ...)
     # GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'MF')
     # GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'CC')

     #################### GSEA ####################
    if (run_GSEA) {
    #  GO_GSEA_analysis(results, local_path, ontology = 'ALL', group = group2)
     GO_GSEA_analysis(results, local_path, group = group2, ...)
     # GO_GSEA_analysis(results, local_path, ontology = 'MF')
     # GO_GSEA_analysis(results, local_path, ontology = 'CC')
    }


######################################## DOWN ########################################

     local_path <- here(path, paste0('GO_functional_analysis_UP_in_', group1, '_', cluster))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> pull(genes)
     all_genes <- results |> arrange(padj) |> pull(genes)

     ######################## ORA ########################

    #  GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'ALL', group = group1)
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, group = group1, ...)
     # GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'MF')
     # GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'CC')

    
     #################### GSEA ####################
     if (run_GSEA) {
        # GO_GSEA_analysis(results, local_path, ontology = 'ALL', group1)
     GO_GSEA_analysis(results, local_path, group1, ...)
     # GO_GSEA_analysis(results, local_path, ontology = 'MF')
     # GO_GSEA_analysis(results, local_path, ontology = 'CC')
     }

     return()
}



GO_functional_analysis_cluster_identification <- function (scRNAseq, results, identities = 'seurat_clusters', path='./', object_annotations, ...) {


    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
     local_path <- here(path, paste0('Cluster_identification_functional_analysis_GO_', object_annotations))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()
    #  all_genes_entrezid <- Features(scRNAseq[['RNA']]) |>unique() |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique()

    #  mouse_database <- msigdbr(species = 'Mus musculus',category = 'C8') |> dplyr::select(gs_name, entrez_gene)

for (cluster in unique(scRNAseq@meta.data |> pull(!!identities))) {

####################################### GO ########################################

     name <- paste0('Cluster ', cluster, ' - ')

     significant_results_cluster <- results |> 
          filter(cluster == cluster) |>
          arrange(desc(avg_log2FC)) |>
          pull(gene) |>
          unique()

     GO_overrepresentation_analysis(significant_results_cluster, all_genes, local_path, ontology = 'ALL',filename = name, group = name, nterms_to_plot = 20, ...)  
 
#################### msigdbr ####################

    #  significant_results_cluster  <-  significant_results_cluster |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique() 
    #  enrichment_results <- enricher(gene = significant_results_cluster, 
    #                 universe = all_genes_entrezid,
    #                 #keyType = "ENTREZID",
    #                 #OrgDb = org.Mm.eg.db, 
    #                 pAdjustMethod = "BH", 
    #                 minGSSize    = 5,
    #                 maxGSSize    = 2000,
    #                 qvalueCutoff = 0.05,
    #                 TERM2GENE = mouse_database)
                        
    #     enrichment_results_table <- as_tibble(enrichment_results)
    #     write.csv(enrichment_results_table, paste0(local_path, name,'MSigDbr_OverRepresentation_analysis_results_', 'C8','.csv'))

    #     if (nrow(enrichment_results_table) > 0) {
    #         ## Add similarity matrix to the termsim slot of enrichment result
    #         enrichment_results <- enrichplot::pairwise_termsim(enrichment_results)
            
    #         enrichment_results_table <- as_tibble(enrichment_results) 
    #         write.csv(enrichment_results_table, paste0(local_path, 'MSigDbr_OverRepresentation_analysis_results_', 'C8','.csv'))

    #         dotplot(enrichment_results,
    #             showCategory=50,
    #             title = paste0('MSigDbr Overrepresentation analysis ', 'C8'),
    #             label_format = 60)
    #         ggsave(paste0(name,'MSigDbr_overrepresentation_analysis_dotplot_', 'C8','.pdf'), width = 10, height = 18, path = local_path)

    #         ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
    #         emapplot(enrichment_results, showCategory = 60) + ggtitle(paste0('Overrepresentation analysis ', 'C8'))
    #         ggsave(paste0(name,'MSigDbr_overrepresentation_analysis_network_', 'C8','.pdf'), width = 14, height = 18, path = local_path)
    #     }
     
}
     return()
}
msigdbr_functional_analysis <- function (results,cluster,  path='./') { 
    
    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
        drop_na(entrezid)

        
    local_path <- paste0(path, 'MSigDb_functional_analysis_', cluster,'/')
    unlink(local_path, recursive = T)
    dir.create(local_path)

    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)

    for (database in c('H', 'C2', 'C3', 'C8')) {

        mouse_database <- msigdbr(species = 'Mus musculus',category = database) |> dplyr::select(gs_name, entrez_gene)

        #### ORA ####

        significant_genes <- filter(results, padj < 0.05) |> arrange(padj) |> pull(entrezid)
        all_genes <- results |> arrange(padj) |> pull(entrezid)

        enrichment_results <- enricher(gene = significant_genes, 
                        universe = all_genes,
                        #keyType = "ENTREZID",
                        #OrgDb = org.Mm.eg.db, 
                        pAdjustMethod = "BH", 
                        minGSSize    = 5,
                        maxGSSize    = 500,
                        qvalueCutoff = 0.05,
                        TERM2GENE = mouse_database)
                        
        enrichment_results_table <- as_tibble(enrichment_results)
        write.csv(enrichment_results_table, paste0(local_path, 'OverRepresentation_analysis_results_', database,'.csv'))

        if (nrow(enrichment_results_table) > 0) {
            ## Add similarity matrix to the termsim slot of enrichment result
            enrichment_results <- enrichplot::pairwise_termsim(enrichment_results)
            
            enrichment_results_table <- as_tibble(enrichment_results) 
            write.csv(enrichment_results_table, paste0(local_path, 'OverRepresentation_analysis_results_', database,'.csv'))

            dotplot(enrichment_results,
                showCategory=50,
                title = paste0('Overrepresentation analysis ', database),
                label_format = 60)
            ggsave(paste0('ORA_dotplot_', database,'.pdf'), width = 10, height = 18, path = local_path)

            ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
            emapplot(enrichment_results, showCategory = 60) + ggtitle(paste0('Overrepresentation analysis ', database))
            ggsave(paste0('ORA_network_', database,'.pdf'), width = 14, height = 18, path = local_path)
        }

        #### GSEA ####

        fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
        names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 5,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    TERM2GENE = mouse_database)
        
        gsea_results_table <- as_tibble(gsea_results)
        write.csv(gsea_results_table, paste0(local_path, 'GSEA_analysis_results_', database,'.csv'))

        if (nrow(gsea_results_table) > 1) {
            gsea_results_table <- as_tibble(gsea_results) 
            write.csv(gsea_results_table, paste0(local_path, 'GSEA_analysis_results_', database,'.csv'))

            dotplot(gsea_results,
                showCategory=50,
                title = paste0('GSEA analysis ', database),
                label_format = 60)
            ggsave(paste0('GSEA_dotplot_', database,'.pdf'), width = 10, height = 18, path = local_path)

            ## Add similarity matrix to thenes,  termsim slot of enrichment result
            gsea_results <- enrichplot::pairwise_termsim(gsea_results)

            ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
            emapplot(gsea_results, showCategory = 60) + ggtitle(paste0('GSEA analysis ', database))
            ggsave(paste0('GSEA_network_', database,'.pdf'), width = 14, height = 18, path = local_path)
        }
        
    }

    
    return()
}
pathways_of_interest_analysis <- function (results,  pathways_of_interest_list, cluster,path='./', FC_threshold = 0.3, p_value_threshold = 0.1, max_overlaps = 1000, label_size = 5, group1, group2, comparison) {

    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
    drop_na(entrezid)
        
    local_path <- paste0(path, 'pathways_of_interest_', cluster, '/')
    unlink(local_path, recursive = T)
    dir.create(local_path)


    # pathways_of_interest <- pathways_of_interest |>
    #     # pivot_longer(cols = everything(),values_to = 'genes', names_to = 'pathway', values_drop_na = TRUE) |>
    #     arrange(desc(pathway)) 

    pathways_of_interest_list <- map(pathways_of_interest_list, \(x) bitr(x, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID))

    #head(pathways_of_interest)

    fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
    names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

    for (term_of_interest in names(pathways_of_interest_list)) {

        term <- pathways_of_interest_list[[term_of_interest]]

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 4,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE,
                    TERM2GENE = term)
                    
        anno <- gsea_results[term_of_interest, c("NES", "pvalue", "p.adjust")]
        lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

        p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 1, base_size = 13,title = term_of_interest)

        x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75
        y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17

        p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5)
        p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 2, base_size = 13)
        p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 3, base_size = 13)     
                
        cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v')    
        ggsave(paste0('GSEA_', term_of_interest, '.pdf'), path = local_path, height = 10, width = 8)

        ########## scatterplot ##########

        ## prepare for visualization
        results <- results %>% 
                        
                        mutate(                        
                            genes_to_label2 = ifelse(entrezid %in% term$genes,  genes ,''),
                            add_label = ifelse(entrezid %in% term$genes,  'YES' ,'NO')
                            ) |>
                        mutate(add_label <- factor(add_label, levels = c('NO', 'YES'))) |>     
                        mutate(fold_change_direction = case_when(
                                                log2FoldChange>= 0 ~ 'UP',
                                                log2FoldChange<= 0 ~ "DOWN",
                                                TRUE ~ 'NO')) |>
                        arrange(add_label)
                         
                            
        my_colors <- c( "#32228b", "gray")
        names(my_colors) <- c("YES", "NO")
        
        # Scatterplot
        limx <- results |> pull(paste0('Avg_', group1)) |> max()
        limy <- results |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*5
        
        results |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label2, col = add_label))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.2,
                    show.legend = FALSE,
                    max.overlaps = max_overlaps,
                    max.time = 30,
                    max.iter = 10000000,
                    nudge_x = ifelse(results$fold_change_direction == 'UP', -0.75, 0.75),
                    nudge_y = ifelse(results$fold_change_direction == 'UP', 0.75, -0.75),
                    aes(segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            
            
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=term_of_interest,
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+       
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(local_path, 'Pseudobulk scatter ', term_of_interest, '.pdf'))
    }
    return()
}


pathways_of_interest_analysis2 <- function (results,  pathways_of_interest_table, cluster,path='./', FC_threshold = 0.3, p_value_threshold = 0.1, max_overlaps = 1000, label_size = 5, group1, group2, comparison) {

    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
    drop_na(entrezid)
        
    local_path <- paste0(path, 'pathways_of_interest_', cluster, '/')
    unlink(local_path, recursive = T)
    dir.create(local_path)


    # pathways_of_interest_table <- pathways_of_interest_table |>
    #     pivot_longer(cols = everything(),values_to = 'genes', names_to = 'pathway', values_drop_na = TRUE) |>
    #     arrange(desc(pathway)) 

    # pathways_of_interest_table$genes <- pathways_of_interest_table |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) 

#     pathways_of_interest_table <- pathways_of_interest_table |> filter(!is.na(genes) )

    #head(pathways_of_interest_table)

    fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
    names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

    for (term_of_interest in unique(pathways_of_interest_table$gs_name)) {

        term <- pathways_of_interest_table |> filter(gs_name == term_of_interest)

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 4,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE,
                    TERM2GENE = term)
                    
        anno <- gsea_results[term_of_interest, c("NES", "pvalue", "p.adjust")]
        lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

        p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 1, base_size = 13,title = term_of_interest)

        x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75
        y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17

        p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5)
        p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 2, base_size = 13)
        p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 3, base_size = 13)     
                
        cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v')    
        ggsave(paste0('GSEA_', term_of_interest, '.pdf'), path = local_path, height = 10, width = 8)

        ########## scatterplot ##########

        ## prepare for visualization
        results <- results %>% 
                        
                        mutate(                        
                            genes_to_label2 = ifelse(entrezid %in% term$genes,  genes ,''),
                            add_label = ifelse(entrezid %in% term$genes,  'YES' ,'NO')
                            ) |>
                        mutate(add_label <- factor(add_label, levels = c('NO', 'YES'))) |>     
                        mutate(fold_change_direction = case_when(
                                                log2FoldChange>= 0 ~ 'UP',
                                                log2FoldChange<= 0 ~ "DOWN",
                                                TRUE ~ 'NO')) |>
                        arrange(add_label)
                         
                            
        my_colors <- c( "#32228b", "gray")
        names(my_colors) <- c("YES", "NO")
        
        # Scatterplot
        limx <- results |> pull(paste0('Avg_', group1)) |> max()
        limy <- results |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*5
        
        results |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label2, col = add_label))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.2,
                    show.legend = FALSE,
                    max.overlaps = max_overlaps,
                    max.time = 30,
                    max.iter = 10000000,
                    nudge_x = ifelse(results$fold_change_direction == 'UP', -0.75, 0.75),
                    nudge_y = ifelse(results$fold_change_direction == 'UP', 0.75, -0.75),
                    aes(segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            
            
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=term_of_interest,
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+       
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(local_path, 'Pseudobulk scatter ', term_of_interest, '.pdf'))
    }
    return()
}