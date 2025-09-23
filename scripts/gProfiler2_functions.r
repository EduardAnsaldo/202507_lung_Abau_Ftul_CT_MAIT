### Over representation analysis -- gProfiler2
gProfiler2_overrepresentation_analysis <- function (significant_genes_FC_ordered, local_path, group, cluster,  filename = '', p_value_threshold = 0.05) {

    significant_genes_FC_ordered <- list(significant_genes_FC_ordered)
    names(significant_genes_FC_ordered) <- paste0('UP in ', group, ' - ', cluster)
    enrichment_results <- gost(query = significant_genes_FC_ordered, 
                    organism = "mmusculus", ordered_query = FALSE, 
                    sources = c('GO:BP','GO:MF','GO:CC','KEGG','REAC','TF'),
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = p_value_threshold, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", as_short_link = FALSE, highlight = TRUE)

    if (!(is.null(enrichment_results[['result']]) )) {

     # Manhattan plot of results
     plot  <- gostplot(enrichment_results, capped = T, interactive=T)
     htmlwidgets::saveWidget((plot), paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.html'))
     enrichment_results2 <- enrichment_results
     enrichment_results2[['result']] <- enrichment_results2[['result']] |> filter((source %in% c('GO:BP', 'GO:CC', 'GO:MF') & highlighted) | !(source %in% c('GO:BP', 'GO:CC', 'GO:MF')))
     plot  <- gostplot(enrichment_results2, capped = T, interactive=T)
     htmlwidgets::saveWidget((plot), paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2_summary_pathways', '.html'))
     rm(enrichment_results2)
     # print(plot)
      
     term_table <- enrichment_results[['result']] |> 
       filter((source %in% c('GO:BP', 'GO:CC', 'GO:MF') & highlighted) | !(source %in% c('GO:BP', 'GO:CC', 'GO:MF'))) |>
       mutate(term_name = str_c(source, term_name, sep = '-')) |>
       arrange(p_value) |>
       head(20) 
      
     # p <- gostplot(enrichment_results, capped = T, interactive=F)
     # publish_gostplot(p , highlight_terms = term_table$term_id, 
     #                   width = NA, height = NA, filename = paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.pdf'))
      
     plot2 <- term_table |>       
     ggplot(aes(x = -log10(p_value), y = fct_reorder(term_name, p_value, .desc = T), fill = -log10(p_value))) +
       geom_col(width = 0.7) +
       scale_fill_viridis_c(option = 'mako', direction = -1)+
       labs(x = '-log10(p-value)', y = '', title = paste0('UP in ', group, ' - ', cluster)) +
       theme_minimal() +
       theme(axis.text.y = element_text(size = 12), title = element_text(size = 16), plot.title.position = 'plot', legend.position = 'none', axis.text.x = element_text(size = 12))
     print(plot2)
     ggsave(plot = plot2, filename = paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.pdf'), width = 12, height = 6)
      
      # Saving results table
     enrichment_results2 <- enrichment_results[['result']] |> dplyr::select(-c('parents'))
     write.csv(enrichment_results2, paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.csv'))

    }else {
     p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A'))+ xlab(NULL)
     print(p1)
    }
    
}

gProfiler2_functional_analysis <- function (results,  cluster, group2, group1, path='./', FC_threshold, p_value_threshold = 0.05) {

     results <- results[which(duplicated(results$genes) == F),]


####################################### UP ########################################

     local_path <- here(path, paste0('Functional_analysis_UP_gProfiler2_', cluster, '_', '/'))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> arrange(desc(log2FoldChange)) |> pull(genes)
     
     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path ,  group = group2, cluster = cluster, filename = '', p_value_threshold = p_value_threshold)

     }else {
          p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A'))+ xlab(NULL)
          print(p1)
     }

######################################## DOWN ########################################

     local_path <- here(path, paste0('Functional_analysis_DOWN_gProfiler2_', cluster, '_', '/'))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> arrange(log2FoldChange) |> pull(genes)

     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path , group = group1, cluster = cluster, filename = '', p_value_threshold = p_value_threshold)
     }
     else {
          p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A'))+ xlab(NULL)
          print(p1)
     }
}

gProfiler2_functional_analysis_cluster_identification <- function (scRNAseq, results, identities = 'seurat_clusters', path='./', object_annotations = '') {

#     color_scale <- viridis(n = 4, direction = -1)
#     options(enrichplot.colours = color_scale)
     local_path <- paste0(path, object_annotations, '_Cluster_identification_functional_analysis_gProfiler2/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()
     #all_genes_entrezid <- Features(scRNAseq[['RNA']]) |>unique() |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique()

     #mouse_database <- msigdbr(species = 'Mus musculus',category = 'C8') |> dplyr::select(gs_name, entrez_gene)

for (cluster in levels(scRNAseq@meta.data |> pull(!!as.name(identities)))) {

     name <- paste0('Cluster ', cluster, ' - ')

     significant_results_cluster <- results |> 
          filter(cluster == {{cluster}}) |>
          pull(gene) |>
          unique()
    
     gProfiler2_overrepresentation_analysis(significant_genes_FC_ordered = significant_results_cluster, local_path = local_path, group = '', filename =  name,  cluster = cluster, p_value_threshold = 0.05)
 

     }
}