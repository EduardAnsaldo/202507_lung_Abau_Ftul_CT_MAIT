## Functions
### Overrepresentation analysis function GO and MSigDB -- ClusterProfiler
GO_overrepresentation_analysis <- function (significant_genes, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '', group = '', drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = T, run_network = F, network_n_terms = 100, nterms_to_plot = 50, font_size = 8,  ...)  {

     color_scale <- viridis(n = 4, direction = -1)
     options(enrichplot.colours = color_scale)
     
     enrichment_results <- enrichGO(gene = significant_genes, 
                    universe = all_genes,
                    keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db, 
                    ont = ontology, 
                    pAdjustMethod = "BH", 
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    qvalueCutoff = 0.25)

     if (drop_levels  == T ) {
            enrichment_results <- dropGO(enrichment_results, level = levels_to_drop)
        
     }

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here(local_path, paste0(filename,'GO_OverRepresentation_analysis_results_', ontology, '.csv')))

     if (nrow(enrichment_results_table) > 1) {
       if (simplify_terms == T) {
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
               label_format = 60, 
               font.size = font_size) 
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

GO_GSEA_analysis <- function (results, local_path, ontology = 'ALL', group) {
     
     color_scale <- viridis(n = 4, direction = -1)
     options(enrichplot.colours = color_scale)

          #### GSEA ####

     fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
     names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(genes)

     gsea_results <- gseGO(geneList     = fold_changes,
               OrgDb        = org.Mm.eg.db,
               ont          = ontology,
               keyType = "SYMBOL",
               minGSSize    = 5,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE,
            eps = 0)               

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



GO_functional_analysis_cluster_identification <- function (scRNAseq, results, path='./', object_annotations = '', top_gene_number = 50, ...) {

    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
     local_path <- here(path, paste0('Cluster_identification_functional_analysis_GO_', object_annotations))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()
    #  all_genes_entrezid <- Features(scRNAseq[['RNA']]) |>unique() |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique()

    gene_lists <- results %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) |>
        slice_head(n = top_gene_number)  |>
        ungroup() |>
        select(cluster, gene)  |>
        group_by(cluster) |>
        mutate(id = row_number()) |>
        ungroup() |>
        pivot_wider(names_from = cluster, values_from = gene)  |>
        select(-id) |>
        as.list() |>
        map(~ .x[!is.na(.x)])
        
    overrepresentation_results <- GO_overrepresentation_analysis_multiple_lists(gene_lists, all_genes = all_genes, local_path = local_path, grouping_var = 'Clusters', ...)  

    return(overrepresentation_results$plot)


####################################### GO ########################################

    #  name <- paste0('Cluster ', cluster, ' - ')

    #  significant_results_cluster <- results |> 
    #       filter(cluster == cluster) |>
    #       arrange(desc(avg_log2FC)) |>
    #       pull(gene) |>
    #       unique()

    #  GO_overrepresentation_analysis(significant_results_cluster, all_genes, local_path, ontology = 'ALL',filename = name, group = name, nterms_to_plot = 20, ...)  
 
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


GO_overrepresentation_analysis_multiple_lists <- function (gene_list, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 400, filename = '',  drop_levels = F, levels_to_drop = c(), simplify_function = min, simplify_by = 'p.adjust', simplify_terms = F, run_network = F, network_n_terms = 100, nterms_to_plot = 5, font_size = 8, grouping_var = '', ...)  {

     color_scale <- viridis(n = 4, direction = -1)
     options(enrichplot.colours = color_scale)
     
     enrichment_results <- compareCluster(
          geneCluster = gene_list,
          fun = "enrichGO",
          universe = all_genes,
          keyType = "SYMBOL",
          OrgDb = org.Mm.eg.db,
          ont = ontology,
          pAdjustMethod = "BH",
          minGSSize    = minGSSize,
          maxGSSize    = maxGSSize,
          qvalueCutoff = 0.2,
          pvalueCutoff = 0.05
     )

     if (drop_levels  == T ) {
            enrichment_results <- dropGO(enrichment_results, level = levels_to_drop)
        
     }
     enrichment_results <- enrichment_results |>
          mutate(Cluster = fct_inseq(Cluster)) 
     #   group_by(Cluster) |>
     #   arrange((p.adjust))# |>
       # ungroup()
     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, here(local_path, paste0(filename,'GO_OverRepresentation_analysis_gene_lists_results_', ontology, '.csv')))

     if (nrow(enrichment_results_table) > 1) {
       if (simplify_terms == T) {
            ## Add similarity matrix to the termsim slot of enrichment result
          enrichment_results <- enrichplot::pairwise_termsim(enrichment_results, showCategory = dim(enrichment_results)[1])
          enrichment_results_unfiltered <-  enrichment_results
          enrichment_results <- clusterProfiler::simplify(enrichment_results, cutoff=0.7, by=simplify_by, select_fun=simplify_function)        
       }

          write.csv(as_tibble(enrichment_results), here(local_path, paste0(filename,'GO_OverRepresentation_analysis_gene_lists_results__filtered_', ontology, '.csv')))
          p1 <- dotplot(enrichment_results,
               showCategory=nterms_to_plot,
               title = paste0(filename,'GO Overrepresentation analysis ', ' - ', ontology),
            #    x = 'p.adjust',
               label_format = 60, 
               font.size = font_size) +
               labs(x = grouping_var)
          print(p1)
          ggsave(plot = p1, filename = paste0(filename, 'GO overrepresentation_analysis_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)
     }
     return(list(plot = p1, results = enrichment_results))
}