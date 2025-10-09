create_group_palette <- function(samples, groups, base_colors = c(
                "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid", 
                "seagreen", "gold", "lavenderblush", "tomato")) {
                # Helper to get shades for a base color
                get_shades <- function(color_name, n) {
                        shades <- c(                            
                                color_name,
                                paste0(color_name, "1"),
                                paste0(color_name, "2"),
                                paste0(color_name, "3"),
                                paste0(color_name, "4")
                        )
                        # If more than 5, use grDevices::colors() to find similar colors
                        if (n > 5) {
                                all_colors <- colors()
                                similar <- grep(color_name, all_colors, value = TRUE)
                                # Remove already used shades
                                similar <- setdiff(similar, shades)
                                shades <- c(shades, head(similar, n - 5))
                        }
                        head(shades, n)
                }
                palette <- setNames(rep(NA, length(samples)), samples)
                for (grp in groups) {
                        grp_samples <- samples[grepl(grp, samples)]
                        n <- length(grp_samples)
                        # Pick base color for group (from base_colors, in order of groups)
                        idx <- which(groups == grp)
                        base_color <- base_colors[idx]
                        shades <- get_shades(base_color, n)
                        palette[grp_samples] <- shades[seq_along(grp_samples)]
                }
                return(palette)
        } 


plot_circos_clonotypes <- function(
    clonotype_data_plot,
    clonotype_data_plot_distinct,
    figures_path,
    grouping_variable,
    variables_to_color_by = NULL,
    file_type = "pdf",
    circle_margin = 0.2,
    major_ticks = 500,
    cex = 0.6,
    samples,
    groups,
    color_palette = NULL,
    alpha_col = 0.8,
    base_colors = c(
        "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid", 
        "seagreen", "gold", "lavenderblush", "tomato")
) {
    cluster <- 'all'
    if (is.null(color_palette)) {
        grid_cols <- create_group_palette(samples, groups, base_colors)
    } else {
        grid_cols <- color_palette
    }

    file_name <- paste0(figures_path, 'Circos_clonotypes_per_', grouping_variable, '.', file_type)
    if (file_type == "pdf") {
        pdf(file = file_name)
    } else if (file_type == "png") {
        png(file = file_name, width = 800, height = 800)
    } else {
        stop("file_type must be either 'pdf' or 'png'")
    }

    circos.par(gap.degree = 2, track.height = 0.1, cell.padding = c(0, 0, 0, 0), circle.margin = circle_margin)
    circos.initialize(xlim = clonotype_data_plot_distinct)

    circos.track(ylim = c(0,1),
        panel.fun = function(x, y) {
            print(CELL_META$xrange[[1]])
            if (CELL_META$cell.width < 45) {
                circos.text(
                    CELL_META$xcenter, 
                    CELL_META$cell.ylim[2] + mm_y(9), 
                    adj = c(0 , 1),
                    CELL_META$sector.index,
                    facing = 'clockwise', 
                    niceFacing = TRUE, 
                    cex  = cex
                )
            } else {
                circos.text(
                    CELL_META$xcenter, 
                    CELL_META$cell.ylim[2] + mm_y(11),
                    CELL_META$sector.index,
                    facing = 'bending.inside', 
                    niceFacing = TRUE,
                    cex = cex + 0.2
                )
            }
            if (CELL_META$xrange[[1]] > major_ticks * 2) {
                circos.axis(
                    labels.cex = 0.5, 
                    minor.ticks =  0,
                    major.tick = 1,
                    labels.facing = 'clockwise',
                    major.at = seq(major_ticks, CELL_META$xrange[[1]], by = major_ticks)
                )
            }
            highlight.sector(CELL_META$sector.index, col = grid_cols[CELL_META$sector.index])
        }
    )

    done <- c()
    for (origin in levels(clonotype_data_plot$group)) {
        for (target in levels(clonotype_data_plot$group)) {
            if (origin == target | target %in% done) {
                next
            } else {
                table_one <- clonotype_data_plot |>
                    filter(group == origin & group_counts != 0)  |>
                    dplyr::select(all_of(c('clonotype', 'group', 'coordinates')))
                table_two <- clonotype_data_plot |>
                    filter(group == target & group_counts != 0) |>
                    dplyr::select(all_of(c('clonotype', 'group', 'coordinates')))
                link_table  <-  inner_join(table_one, table_two, by = 'clonotype') |> column_to_rownames(var = 'clonotype')
                for (clonotype1 in rownames(link_table)) {
                    group1 <- link_table[[clonotype1, 'group.x']]
                    group2 <- link_table[[clonotype1, 'group.y']]
                    coordinates1 <- as.vector(link_table[[clonotype1, 'coordinates.x']])
                    coordinates2 <- as.vector(link_table[[clonotype1, 'coordinates.y']])
                    if (!is.null(variables_to_color_by)) {
                        for (variable in variables_to_color_by) {
                            if (str_detect(origin, variable) | str_detect(target, variable)) {
                                color = alpha(grid_cols[variable], alpha_col)
                                break
                            } else {
                                color = alpha("gray", 0.25)                                                               
                            }                        
                        }
                    } else {
                        color = alpha(grid_cols[origin], alpha_col)                                                               
                    }
                    circos.link(
                        origin, 
                        coordinates1,
                        target,
                        coordinates2,
                        col = color
                    )
                }
            }
        }
        done <- c(done, origin)
    }
    title(paste0('Clonotype Overlap ', cluster, ' per ', grouping_variable))
    dev.off()
    circos.clear()
}

overlap_circos_and_tables <- function(
    seurat,
    grouping_variable = 'Samples',
    results_path = './',
    variables_to_color_by = NULL,
    cell_types_column = 'cell_types',
    write_table = TRUE,
    figures_path,
    circle_margin = 0.2,
    major_ticks = 500,
    cex = 0.6,
    samples,
    groups,
    color_palette = NULL,
    alpha_col = 0.8,
    base_colors = c(
        "skyblue", "darkorange", "wheat", "darkolivegreen", "orchid", 
        "seagreen", "gold", "lavenderblush", "tomato")
) {
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    
    combined2 <- scRepertoire:::.expression2List(seurat, split.by = 'orig.ident')
    
    clonotype_data <- combined2[[1]] |>
        as_tibble() |>
        dplyr::select(all_of(c('CTaa', 'CTgene', grouping_variable, cell_types_column))) |>
        add_count(CTaa, !!as.name(grouping_variable), sort = TRUE, name = 'counts_per_condition') |>
        group_by(CTaa, !!as.name(grouping_variable)) |>
        summarize(across(everything(), Mode), .groups = 'drop') |>
        arrange(desc(counts_per_condition)) |>
        pivot_wider(
            id_cols = c('CTaa'),
            names_from = !!as.name(grouping_variable),
            values_from = counts_per_condition,
            unused_fn = Mode,
            values_fill = 0
        ) |>
        ungroup() |>
        mutate(clonotype = paste0('clonotype ', as.character(row_number())))
    
    if (write_table) {
        write_csv(clonotype_data, file = paste0(results_path, 'TCR_data_per_group.csv'))
    }
    
    group_levels <- levels(pull(seurat@meta.data, grouping_variable))
    
    clonotype_data_plot <- clonotype_data |>
        pivot_longer(
            c(group_levels),
            names_to = 'group',
            values_to = 'group_counts'
        ) |>
        group_by(group) |>
        mutate(sequence_count_by_grouping_variable = sum(group_counts)) |>
        arrange(desc(group_counts)) |>
        mutate(clonotype_position_on_circos_by_grouping_variable = cumsum(group_counts)) |>
        mutate(clonotype_start_position_by_grouping_variable = c(
            0,
            clonotype_position_on_circos_by_grouping_variable[-length(clonotype_position_on_circos_by_grouping_variable)]
        )) |>
        ungroup() |>
        rowwise() |>
        mutate(coordinates = list(c(
            clonotype_start_position_by_grouping_variable,
            clonotype_position_on_circos_by_grouping_variable
        ))) |>
        ungroup() |>
        mutate(group = fct(group, levels = group_levels))
    
    clonotype_data_plot_distinct <- clonotype_data_plot |>
        dplyr::select(all_of(c('group', 'sequence_count_by_grouping_variable'))) |>
        distinct() |>
        arrange(group) |>
        mutate(origin = 0) |>
        column_to_rownames('group') |>
        relocate(origin) |>
        as.matrix()
    
    plot_circos_clonotypes(
        clonotype_data_plot = clonotype_data_plot,
        clonotype_data_plot_distinct = clonotype_data_plot_distinct,
        figures_path = figures_path,
        grouping_variable = grouping_variable,
        variables_to_color_by = variables_to_color_by,
        file_type = 'pdf',
        circle_margin = circle_margin,
        major_ticks = major_ticks,
        cex = cex,
        samples = samples,
        groups = groups,
        alpha_col = alpha_col,
        base_colors = base_colors,
        color_palette = color_palette
    )

    plot_circos_clonotypes(
        clonotype_data_plot = clonotype_data_plot,
        clonotype_data_plot_distinct = clonotype_data_plot_distinct,
        figures_path = figures_path,
        grouping_variable = grouping_variable,
        variables_to_color_by = variables_to_color_by,
        file_type = 'png',
        circle_margin = circle_margin,
        major_ticks = major_ticks,
        cex = cex,
        samples = samples,
        groups = groups,
        base_colors = base_colors,
        color_palette = color_palette
    )

    return(clonotype_data)
}