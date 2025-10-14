### Over representation analysis -- Metascape
Metascape_overrepresentation_analysis <- function (significant_genes_FC_ordered, path, group, grouping_var,  filename = '') {

#     # Create local folder
#     local_path <- here(path, 'Metascape_overrepresentation_analysis')
#     unlink(path, recursive = T)
#     dir.create(local_path, recursive = T, showWarnings = F)
     local_path <- path
    
    # Create folder in msbio_v3-2
    project_name <- here() |> basename()
    msbio_path <- paste0('~/msbio_v3-2/', 'data/', project_name, '/', basename(path))
    msbio_path_short <- paste0('data/', project_name, '/', basename(path))
    dir.create(msbio_path, recursive = T)
    msbio_file <- paste0(msbio_path, '/', filename, '_', grouping_var, '_UP_in_', group, '_gene_list.txt')
    msbio_file_path_short <- paste0(msbio_path_short, '/', filename, '_', grouping_var, '_UP_in_', group, '_gene_list.txt')

    # Write gene list to file in msbio folder
    writeLines(significant_genes_FC_ordered, msbio_file)

    # Run bash command
    bash_command <- paste0("cd ~/msbio_v3-2/; ls; bin/ms.sh -up -o", msbio_path_short, " ", msbio_file_path_short)
    system(bash_command)
    
    # Copy output back to current working directory
    # Note: You may need to adjust this based on what output files are created
    output_files <- list.files(path = msbio_path, full.names = TRUE)
    file.copy(output_files, local_path, recursive = TRUE)

   if (!file.exists(here(local_path, 'metascape_result.xlsx'))) {
     p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A')) + theme(text = element_text(size = 24)) + xlab(NULL)
          print(p1)
     return()
   }       
     
  # Read results
  enrichment_results <- read_excel(here(local_path, 'metascape_result.xlsx'), sheet = 'Enrichment')

  # Filter and prepare for plotting
    enrichment_results <- enrichment_results |>
      filter(str_detect(GroupID, 'Summary') ) |>
      #select(Description, `Log(q-value)`) |>
      select(Description, LogP) |>
      rename(term = Description, log_q_value = LogP)  |>
      mutate(log_q_value = -1*log_q_value) |> 
      # rename(term = Description, log_q_value = LogP)  |>
      head(n = 20)

    plot2 <- enrichment_results |>       
        ggplot(aes(x = log_q_value, y = fct_reorder(term, log_q_value , .desc = F), fill = log_q_value)) +
          geom_col(width = 0.7) +
          scale_fill_binned(palette = hcl.colors(n = 6, 'YlOrRd', rev = T), breaks = c(2, 4, 6, 10, 20)) +
          labs(x = '-log10(p-value)', y = '', title = paste0('UP in ', group, ' - ', grouping_var)) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 9), title = element_text(size = 16), plot.title.position = 'plot', legend.position = 'none', axis.text.x = element_text(size = 12))
        print(plot2)
         ggsave(plot = plot2, filename = paste0(filename, 'Pathway_enrichment_analysis_metascape', '.pdf'), width = 12, height = 6, path = local_path)
  
}


Metascape_functional_analysis <- function (results, grouping_var, group2, group1, path='./', FC_threshold, p_value_threshold = 0.05) {

     results <- results[which(duplicated(results$genes) == F),]

####################################### UP ########################################

     local_path <- here(path, paste0('Functional_analysis_metascape_UP_', grouping_var))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> arrange(desc(log2FoldChange)) |> pull(genes)
     
     if (length(significant_genes) > 2) {
          Metascape_overrepresentation_analysis(significant_genes, path =  local_path , group = group2, grouping_var = grouping_var, filename = '')

     }else {
          p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A')) + theme(text = element_text(size = 24)) + xlab(NULL)
          print(p1)
     }

######################################## DOWN ########################################

     local_path <- here(path, paste0('Functional_analysis_metascape_DOWN_', grouping_var))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> arrange(log2FoldChange) |> pull(genes)

     if (length(significant_genes) > 2) {
          Metascape_overrepresentation_analysis(significant_genes, path = local_path , group = group1, grouping_var = grouping_var, filename = '')
     }
     else {
          p1 <- ggplot()+theme_void()+ geom_text(aes(0,0,label='N/A')) + theme(text = element_text(size = 24)) + xlab(NULL)
          print(p1)
     }
}

Metascape_functional_analysis_cluster_identification <- function (seurat, results, identities = 'seurat_clusters', path='./', object_annotations = '') {

     local_path <- here(path, paste0(object_annotations, '_Cluster_identification_functional_analysis'))
     unlink(local_path, recursive = T)
     dir.create(local_path)

for (cluster in levels(seurat@meta.data |> pull(!!as.name(identities)))) {

     name <- paste0('Cluster_', cluster)

     significant_results_cluster <- results |> 
          filter(cluster == {{cluster}}) |>
          pull(gene) |>
          unique()
     

  Metascape_overrepresentation_analysis(significant_genes_FC_ordered = significant_results_cluster, path = local_path, group = '', filename =  name,  grouping_var = cluster) 

  }
}