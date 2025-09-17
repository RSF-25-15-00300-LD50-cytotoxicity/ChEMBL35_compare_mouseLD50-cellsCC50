library(tidyverse)
library(patchwork)

# Input
cells <- read_tsv(".../output/cells_completed.tab") |> rename(source = progeny)
data  <- read_tsv(".../output/cell_data.tab")
cell_data <- data |> select(assay, molregno, tid) |> inner_join(cells)
# Read the normalizeed prices of the cells, values were mainly obtained manually from ATCC website (https://www.atcc.org/) on 15/09/2025
cell_norm_price_raw <- read_tsv(".../output/cells_price_normd.tab")

# Summarize the data by cells' source
source_summary <- cell_data |> select(source) |> group_by(source) |> summarize(n_in_source = n()) |> arrange(desc(n_in_source)) |> mutate(source = fct_reorder(source, n_in_source))
# Summarize the data by cells' type
type_summary <- cell_data |> select(type) |> group_by(type) |> summarize(n_in_type = n()) |> arrange(desc(n_in_type)) |> mutate(type = fct_reorder(type, n_in_type))
# Summarize the data by the developmental state of origins of the cells
dvlp_summary <- cell_data |> select(developmental) |> group_by(developmental) |> summarize(n_in_stage = n()) |> arrange(desc(n_in_stage))
# Summarize the data by cell, divide into two parts: top cells (each cell mentioned in > ~1% of total records) and other cells
cell_summary <- cell_data |> select(tid) |> group_by(tid) |> summarize(n_in_cell = n()) |> arrange(desc(n_in_cell))

# Combine cell summaries to express the graph structure for the further visualization
cell_summary_all <- cell_summary |> inner_join(cells |> select(tid, cell_name, source, type, developmental)) |>
                                                        inner_join(source_summary) |>
                                                        inner_join(type_summary) |>
                                                        inner_join(dvlp_summary) |>
                                                        select(tid, cell_name, developmental, type, source,  n_in_cell, n_in_stage, n_in_type, n_in_source) |>
                                                        mutate(developmental = if_else(is.na(developmental), "probably not", "developmental")) |> # human readable variant
                                                        mutate(count_group = if_else(n_in_cell / nrow(cell_data) > .01, "top", "others")) |> # 67 is too many cells, other cells could be agregated
                                                        group_by(count_group) |>
                                                        mutate(n_in_count_group = sum(n_in_cell)) |>
                                                        ungroup()

# Update some distinct summaries
dvlp_summary <- cell_summary_all |> select(developmental, n_in_stage) |> distinct()
cell_summary_short <- cell_summary_all |> select(cell_name, n_in_cell, count_group, n_in_count_group) |>
                                          mutate(cell_name = if_else(count_group == "others", count_group, cell_name),
                                                  n_in_cell = if_else(count_group == "others", n_in_count_group, n_in_cell)) |>
                                          select(cell_name, n_in_cell) |>
                                          distinct() |>
                                          arrange(desc(n_in_cell))

# Add the count group to the cells' prices
cell_norm_price <- cell_norm_price_raw |> right_join(cell_summary_all) |>
                                      select(tid, cell_name, normd_price, count_group) |>
                                      mutate(count_group = if_else(count_group == "top", cell_name, "others")) |>
                                      filter(!is.na(count_group))
# Calculate confidence interval for cells' prices, SEE: https://r-graph-gallery.com/4-barplot-with-error-bar.html
mult_prices <- cell_norm_price |> filter(count_group == "others" & !is.na(normd_price)) |> pull(normd_price)
others_mean <- mean(mult_prices)
sd <- sd(mult_prices)
se = sd(mult_prices) / sqrt(length(mult_prices))
alpha<-.05
t=qt((1-alpha)/2 + .5, length(mult_prices)-1)
ci=t*se
others_price <- tibble(tid = NA, cell_name = "others", normd_price = others_mean, ci = ci)
cell_price <- cell_norm_price |> mutate(ci = NA) |> filter(count_group != "others") |> select(-count_group)
cell_summary_short <- cell_summary_short |> left_join(bind_rows(cell_price, others_price)) |>
                                            mutate(cell_name = fct_reorder(cell_name, n_in_cell))

# Visualize summaries
# Color palette for cells:
cell_color <- c("#b71c1c", "#fb8c00", "#78909c", "#f44336", "#e57373", "#e57373", "#64b5f6", "#ff8a80", "#ff8a80", "#ff8a80", "#ffcc80", "#ff8a80", "#ffe0b2", "#ffcdd2", "#ffcdd2", "#ffcdd2", "#ffe0b2", "#ffcdd2") |> rev()
# Prices, SEE: https://r-graph-gallery.com/4-barplot-with-error-bar.html
price_plot <- ggplot(cell_summary_short) +
      geom_segment( aes(x=cell_name, xend=cell_name, y=0, yend=normd_price), color="black" ) +
      geom_point( aes(x=cell_name, y=normd_price, fill=cell_name), shape=21, size=4, color = "black" ) +
      scale_fill_manual(values = cell_color) +
      geom_errorbar( aes(x=cell_name, ymin=normd_price-ci, ymax=normd_price+ci), width=0.4, colour="black", alpha=1, linewidth=0.6 ) +
      geom_hline( yintercept = 0, linewidth=0.2 ) +
      theme_classic() +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none") +
      labs(y = "Normalized Price")
price_plot
records_plot <- ggplot(cell_summary_short) +
      geom_bar( aes(x=cell_name, y=n_in_cell, fill=cell_name), stat="identity" ) +
      scale_fill_manual(values = cell_color) +
      scale_y_continuous( transform = "sqrt", breaks = c(40, 100, 200, 400, 1500) ) +
      geom_hline( yintercept = 40, linetype = 5, linewidth=0.2 ) +
      geom_hline( yintercept = 100, linetype = 5, linewidth=0.2 ) +
      geom_hline( yintercept = 200, linetype = 5, linewidth=0.2 ) +
      geom_hline( yintercept = 400, linetype = 5, linewidth=0.2 ) +
      geom_hline( yintercept = 1500, linetype = 5, linewidth=0.2 ) +
      theme_classic() +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
      labs(x = "Cell", y = "Number of Measurements")
records_plot
price_plot / records_plot
# Developmental stage
dvlp_plot <- ggplot(dvlp_summary) +
              geom_bar( aes(x=developmental, y=n_in_stage, fill=developmental), stat="identity" ) +
              scale_fill_manual( values = c("#8FB9A8", "#765D69") ) +
              scale_y_continuous(breaks = c(233, 3389)) +
              theme_classic() +
              theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                    legend.position = "none") +
              labs(y = "Number of Measurements")
dvlp_plot
# Cell type
type_plot <- ggplot(type_summary) +
              geom_bar( aes(x=type, y=n_in_type, fill=type), stat="identity" ) +
              scale_fill_manual( values = c("#64b5f6", "#ffcc80", "#e57373") ) +
              scale_y_continuous(breaks = c(174, 903, 2545)) +
              theme_classic() +
              theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(),
                    legend.position = "none") +
              labs(y = "Number of Measurements")
type_plot
# Cell source
source_plot <- ggplot(source_summary) +
              geom_bar( aes(x=source, y=n_in_source, fill=n_in_source), stat="identity" ) +
              scale_fill_gradient(low = "#b2dfdb", high = "#004d40") +
              scale_y_continuous( transform = "sqrt", breaks = c(1, 10, 30, 100, 150, 300, 900, 2000) ) +
              geom_hline( yintercept = 10, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 30, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 100, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 150, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 300, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 900, linetype = 5, linewidth=0.2 ) +
              geom_hline( yintercept = 2000, linetype = 5, linewidth=0.2 ) +
              theme_classic() +
              theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                    legend.position = "none") +
              labs(y = "Number of Measurements")
source_plot
# Compose the plots using patchwork
cell_plot <- ( price_plot / records_plot ) | ( (dvlp_plot | type_plot) / source_plot )
cell_plot
# Save the composed plot
ggsave(".../output/cell_description__17-09-2025.png", plot = cell_plot, width = 10, height = 6, units = "in", dpi = 300)


# Prepare the data for Sankey diagram. Numerous realizations of Sankey diagram are available in R, but at the first glance they do not provide enough control.
# So, D3.js could be to build Sankey diagram, SEE: https://observablehq.com/@d3/sankey/2
# D3.js takes the data as JSON-object which includes array of Nodes and array of Links OR links in csv-format 
# Nodes: name, category
# Links: start, end, value +type (to distinguish between cancer and not cancer in this case), color of link and colors of nodes
# Colors to fill the nodes
color_code <- tibble( entity = c('Huh-7', 'MT4', 'others', 'HepG2', 'THP-1', 'Huh-7.5', 'MRC5', 'RD', 'A549', 'Caco-2', 'HEK-293T', 'HeLa', 'MT2', 'TZM', 'L02', 'CCRF-CEM', 'HEK293', 'HepG2 2.2.15',
                                  "developmental", "probably not", 
                                  "Liver", "Immune", "Lung", "Reproduction", "Kidney", "Gut", "Musculoskeletal", "Breast", "Skin", "Oral and Craniofacial", "Nervous system", "Eye", "Heart and Vascular", "Pancreas", "Spleen"),
                      node_color  = c("#b71c1c", "#fb8c00", "#78909c", "#f44336", "#e57373", "#e57373", "#64b5f6", "#ff8a80", "#ff8a80", "#ff8a80", "#ffcc80", "#ff8a80", "#ffe0b2", "#ffcdd2", "#ffcdd2", "#ffcdd2", "#ffe0b2", "#ffcdd2",
                                  "#8FB9A8", "#765D69",
                                  "#004d40", "#64988f", "#98c7c2", "#98c7c2", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#a6d4cf", "#b2dfdb"
                                  ))
# Links:
stage_links_cell <- cell_summary_all |> rename(start = developmental) |>
                                            mutate(end = if_else(count_group == "others", count_group, cell_name)) |>
                                            uncount(n_in_cell) |>
                                            group_by(start, end, type) |>
                                            mutate(value = n()) |>
                                            ungroup() |>
                                            select(start, end, value, type) |>
                                            distinct() |>
                                            mutate(start_type = "stage") |>
                                            mutate(end_type = "cell") |>
                                            inner_join(color_code, by = c("start" = "entity")) |>
                                            inner_join(color_code, by = c("end" = "entity")) |>
                                            rename(source_color = node_color.x, target_color = node_color.y)
cell_links_source <- cell_summary_all |> rename(start = cell_name, end = source, value = n_in_cell) |>
                                            mutate(start = if_else(count_group == "others", count_group, start)) |>
                                            select(start, end, value, type) |>
                                            uncount(value) |>
                                            group_by(start, end, type) |>
                                            mutate(value = n()) |>
                                            ungroup() |>
                                            distinct() |>
                                            select(start, end, value, type) |>
                                            mutate(start_type = "cell") |>
                                            mutate(end_type = "source") |>
                                            inner_join(color_code, by = c("start" = "entity")) |>
                                            inner_join(color_code, by = c("end" = "entity")) |>
                                            rename(source_color = node_color.x, target_color = node_color.y)
# Gather data
cell_links <- bind_rows(stage_links_cell, cell_links_source) |> mutate(link_color = case_match(
                                                                                type,
                                                                                "cancer" ~ "#e57373",
                                                                                "finite" ~ "#64b5f6",
                                                                                "transformed" ~ "#ffcc80"
                                                                                )) |>
                                                              rename(source = start, target = end)
# Write CSV
write_csv(cell_links, ".../output/cell-links_sankey.csv")