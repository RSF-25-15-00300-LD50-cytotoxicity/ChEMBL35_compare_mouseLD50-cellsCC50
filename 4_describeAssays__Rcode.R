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

# Input cells
hcl <- read_tsv(".../output/cells_completed.tab") |> mutate(tid = as.integer(tid))
# Input assays
assay_hcl <- read_tsv(".../output/hcl_tox_assays.tab")
assay_mouse <- read_tsv(".../output/mouse_tox_assays.tab")
assay_hcl_vec <- assay_hcl |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_hcl_vec <- str_c("'", assay_hcl_vec, "'")
assay_mouse_vec <- assay_mouse |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_mouse_vec <- str_c("'", assay_mouse_vec, "'")
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
data <- hcl |> left_join(assay_cell__description) |> select(-clo_id, -efo_id, -cl_lincs_id,
                                                            -cell_ontology_id, -cell_description, -cell_source_tissue,
                                                            -progeny, -developmental, -type_ref)
sources <- data |> select(pubmed_id, doi, journal, year, title) |> distinct()

# Export the data
write_tsv(data, ".../output/cells_&_sources.tab")
write_tsv(sources, ".../output/sources.tab")

# Close the connection
dbDisconnect(con)