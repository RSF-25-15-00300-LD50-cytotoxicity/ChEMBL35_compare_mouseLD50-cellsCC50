library(tidyverse)
library(RMariaDB)
library(DBI)

# Connect
mysql_password = '***'
con <- dbConnect(
  drv = RMariaDB::MariaDB(),
  dbname = 'chembl_35',
  username = 'root',
  password = mysql_password,
  host = NULL, 
  port = 3306
)
#Check the connection
dbListTables(con)

# Input, cells
hcl <- read_tsv(".../output/cells_completed.tab") |> mutate(tid = as.integer(tid))
# Input, assays
assay_hcl <- read_tsv(".../output/hcl_tox_assays.tab")
assay_mouse <- read_tsv(".../output/mouse_tox_assays.tab")
assay_hcl_vec <- assay_hcl |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_hcl_vec <- str_c("'", assay_hcl_vec, "'")
assay_mouse_vec <- assay_mouse |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_mouse_vec <- str_c("'", assay_mouse_vec, "'")
# Input, distinct records from assays
hcl_records <- read_tsv(".../output/hcl_tox_assays.tab")
mouse_records <- read_tsv(".../output/mouse_tox_assays.tab")
# Formulate queries
assay_cell__query <- dbSendQuery(con, str_glue("SELECT a.chembl_id, a.tid, a.assay_type, a.assay_test_type, a.assay_category,
                                                  a.doc_id, s.src_id, s.src_url, a.description, d.pubmed_id, d.doi,
                                                  d.patent_id, d.journal, d.year, d.title, d.abstract, s.src_short_name, s.src_description, s.src_comment
                                                  FROM assays a
                                                  JOIN docs d
                                                  JOIN source s
                                                  WHERE a.chembl_id IN ({assay_hcl_vec}) AND 
                                                    a.doc_id = d.doc_id AND 
                                                    d.src_id = s.src_id"))
assay_cell__description <- dbFetch(assay_cell__query) |> distinct() |> rename(assay = chembl_id)
dbClearResult(assay_cell__query)
assay_mouse__query <- dbSendQuery(con, str_glue("SELECT a.chembl_id, a.tid, a.assay_type, a.assay_test_type, a.assay_category,
                                                  a.doc_id, s.src_id, s.src_url, a.description, d.pubmed_id, d.doi,
                                                  d.patent_id, d.journal, d.year, d.title, d.abstract, s.src_short_name, s.src_description, s.src_comment
                                                  FROM assays a
                                                  JOIN docs d
                                                  JOIN source s
                                                  WHERE a.chembl_id IN ({assay_mouse_vec}) AND 
                                                    a.doc_id = d.doc_id AND 
                                                    d.src_id = s.src_id"))
assay_mouse__description <- dbFetch(assay_mouse__query) |> distinct() |> rename(assay = chembl_id)
dbClearResult(assay_mouse__query)

# Join cells and assays
data_hcl <- hcl |> left_join(assay_cell__description) |> select(-clo_id, -efo_id, -cl_lincs_id,
                                                            -cell_ontology_id, -cell_description, -cell_source_tissue,
                                                            -progeny, -developmental, -type_ref)
sources_hcl <- data_hcl |> select(pubmed_id, doi, journal, year, title) |> distinct()

# Rename mouse
data_mouse <- assay_mouse__description

# Add mouse data to the mouse assays -> empty
mouse_hcl_intersection <- data_hcl |> select(doi) |> inner_join(data_mouse |> select(doi))

# Check years
hcl_years <- data_hcl |> select(assay, year) |> mutate(object = 'cell')
mouse_years <- data_mouse |> select(assay, year) |> mutate(object = 'mouse')
years <- bind_rows(hcl_years, mouse_years) |> group_by(year, object) |>
                    summarise(n=n(), .groups = "drop")
# Plot years
assay_year_plot <- ggplot(years) +
                    geom_segment( aes(x=year, xend=year, y=0, yend=n), color="black" ) +
                    geom_point( aes(x=year, y=n, fill=object), shape=21, size=4, color = "black" ) +
                    scale_x_continuous(breaks = c(1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022, 2025), labels = c("1975", "1980", "1985", "1990", "1995", "2000", "2005", "2010", "2015", "2020", "2\n0\n2\n2", "2025")) +
                    scale_y_sqrt(breaks = c(1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130) ) +
                    scale_fill_manual(values = c("#8FB9A8", "#765D69")) +
                    theme_classic() +
                    labs(y = "Number of assays")
assay_year_plot

# Check years in relation to the numbers of studied compounds
hcl_years_compounds <- hcl_years |> inner_join(hcl_records) |> select(year, object, molregno) |> distinct()
mouse_years_compounds <- mouse_years |> inner_join(mouse_records) |> select(year, object, molregno) |> distinct()
years_compounds <- bind_rows(hcl_years_compounds, mouse_years_compounds) |>
                    group_by(year, object) |>
                    summarise(n=n(), .groups = "drop")
# Plot years with distinct compounds counted
assay_year_compound_plot <- ggplot(years_compounds) +
                    geom_segment( aes(x=year, xend=year, y=0, yend=n), color="black" ) +
                    geom_point( aes(x=year, y=n, fill=object), shape=21, size=4, color = "black" ) +
                    scale_x_continuous(breaks = c(1975, 1980, 1985, 1990, 1995, 2000, 2005, 2008, 2010, 2015, 2020, 2022, 2025), labels = c("1975", "1980", "1985", "1990", "1995", "2000", "2005", "2\n0\n0\n8", "2010", "2015", "2020", "2\n0\n2\n2", "2025") ) +
                    scale_y_sqrt(breaks = c(1, 5, 10, 50, 100, 200, 400, 600, 800, 1000, 1200, 1400) ) +
                    scale_fill_manual(values = c("#8FB9A8", "#765D69")) +
                    theme_classic() +
                    labs(y = "Number of assays")
assay_year_compound_plot


# Export the results
write_tsv(data_hcl, ".../output/cells_&_sources.tab")
write_tsv(sources_hcl, ".../output/sources.tab")
ggsave(".../output/n_of_experiments__by_year&object.png", plot = assay_year_plot, width = 7, height = 5, units = "in", dpi = 300)
ggsave(".../output/n_of_compounds__by_year&object.png", plot = assay_year_compound_plot, width = 7, height = 5, units = "in", dpi = 300)

# Close the connection
dbDisconnect(con)