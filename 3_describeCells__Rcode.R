library(tidyverse)
library(ggsankey)
library(jsonlite)

# Input
cells <- read_tsv(".../output/cells_completed.tab") |> rename(source = progeny)
data  <- read_tsv(".../output/cell_data.tab")
cell_data <- data |> select(assay, molregno, tid) |> inner_join(cells)

# Summarize the data by cells' source
source_summary <- cell_data |> select(source) |> group_by(source) |> summarize(n_in_source = n()) |> arrange(desc(n_in_source))
# Summarize the data by cells' type
type_summary <- cell_data |> select(type) |> group_by(type) |> summarize(n_in_type = n()) |> arrange(desc(n_in_type))
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
                                                        mutate(developmental = if_else(is.na(developmental), "probably not developmental", "developmental")) |> # human readable variant
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

# Prepare the data for Sankey diagram. Numerous realizations of Sankey diagram are available in R, but at the first glance they do not provide enough control.
# So, D3.js could be to build Sankey diagram, SEE: https://observablehq.com/@d3/sankey/2
# D3.js takes the data as JSON-object which includes array of Nodes and array of Links OR links in csv-format 
# Nodes: name, category
# Links: start, end, value +type (to distinguish between cancer and not cancer in this case)
# Nodes:
cell_nodes_cell <- cell_summary_all |> rename(name = cell_name) |>
                                        mutate(name = if_else(count_group == "others", count_group, name)) |>  
                                        mutate(category = "cell") |>
                                        mutate(layer = 1) |>
                                        select(name, category, layer) |>
                                        distinct()
cell_nodes_stage <- cell_summary_all |> rename(name = developmental) |>        
                                        mutate(category = "stage") |>
                                        mutate(layer = 2) |>
                                        select(name, category, layer) |>
                                        distinct()
cell_nodes_type <- cell_summary_all |> rename(name = type) |>        
                                        mutate(category = "type") |>
                                        mutate(layer = 3) |>
                                        select(name, category, layer) |>
                                        distinct()
cell_nodes_source <- cell_summary_all |> rename(name = type) |>        
                                        mutate(category = "source") |>
                                        mutate(layer = 4) |>
                                        select(name, category, layer) |>
                                        distinct()
cell_nodes <- bind_rows(cell_nodes_cell, cell_nodes_stage, cell_nodes_type, cell_nodes_source)
# Links:
cell_links_cell <- cell_summary_all |> rename(start = cell_name, end = developmental, value = n_in_cell) |>
                                            mutate(start = if_else(count_group == "others", count_group, start)) |>
                                            select(start, end, value, type) |>
                                            uncount(value) |>
                                            group_by(start, end, type) |>
                                            mutate(value = n()) |>
                                            ungroup() |>
                                            distinct() |>
                                            select(start, end, value, type)
cell_links_stage <- cell_summary_all |> rename(start = developmental) |>
                                            mutate(end = type) |>
                                            uncount(n_in_cell) |>
                                            group_by(start, end, type) |>
                                            mutate(value = n()) |>
                                            ungroup() |>
                                            select(start, end, value, type) |>
                                            distinct()
cell_links_type <- cell_summary_all |> rename(end = source) |>
                                            mutate(start = type) |>
                                            uncount(n_in_cell) |>
                                            group_by(start, end, type) |>
                                            mutate(value = n()) |>
                                            ungroup() |>
                                            select(start, end, value, type) |>
                                            distinct()
cell_links <- bind_rows(cell_links_cell, cell_links_stage, cell_links_type)
# Gather data to JSON
cell_json <- list(cell_nodes, cell_links) |> toJSON()

# Write JSON
write_json(cell_json, ".../output/cell_sankey.json")
# Write CSV
write_tsv(cell_links |> mutate(color = if_else(type == "cancer", "#D8737F", "#475C7A")), ".../output/cell-links_sankey.tsv")