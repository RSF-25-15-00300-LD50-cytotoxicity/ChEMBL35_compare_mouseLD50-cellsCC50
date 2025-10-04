library(tidyverse)
library(RMariaDB)
library(DBI)

# Connect
mysql_password = '*****'
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

# Input, assays
assay_hcl <- read_tsv(".../output/hcl_tox_assays.tab")
assay_mouse <- read_tsv(".../output/mouse_tox_assays.tab")
assay_hcl_vec <- assay_hcl |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_hcl_vec <- str_c("'", assay_hcl_vec, "'")
assay_mouse_vec <- assay_mouse |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_mouse_vec <- str_c("'", assay_mouse_vec, "'")
# Formulate queries
assayParameters_cell__query <- dbSendQuery(con, str_glue("SELECT a.chembl_id, a.description, p.standard_type, p.standard_relation, p.standard_value, p.standard_units, p.standard_text_value, p.comments
                                                  FROM assays a
                                                  LEFT JOIN assay_parameters p
                                                  ON a.assay_id = p.assay_id
                                                  WHERE a.chembl_id IN ({assay_hcl_vec})"))
assay_cell__parameters <- dbFetch(assayParameters_cell__query) |> distinct() |> rename(assay = chembl_id) |> mutate(object = "cell")
dbClearResult(assayParameters_cell__query)
assayParameters_mouse__query <- dbSendQuery(con, str_glue("SELECT a.chembl_id, a.description, p.standard_type, p.standard_relation, p.standard_value, p.standard_units, p.standard_text_value, p.comments
                                                  FROM assays a
                                                  LEFT JOIN assay_parameters p
                                                  ON a.assay_id = p.assay_id
                                                  WHERE a.chembl_id IN ({assay_mouse_vec})"))
assay_mouse__parameters <- dbFetch(assayParameters_mouse__query) |> distinct() |> rename(assay = chembl_id) |> mutate(object = "mouse")
dbClearResult(assayParameters_mouse__query)

# Check the results
assay_cell__parameters      # -> empty
assay_mouse__parameters     # -> 5 assays have defined route of administration, 16 assays have defined administered dose in addition to the LD50 provided in the distinct field
# NB this work was done using the data from ChEMBL v35, ChEMBL v36 is out now and the results could be different.
# Some observations:
# Route of administration has game changing significance in the case of MOUSE and should be considered in the actual modeling.
# CELL assays also may have important parameters hidden in their textual descriptions: for example, time of exposure, since different compounds may have different mechanism of action requring different time to lead to the effect.

# Check assay descriptions
assay_cell__description <- assay_cell__parameters |> select(assay, description, object) |> distinct()
assay_mouse__description <- assay_mouse__parameters |> select(assay, description, object) |> distinct()
# Mouse, very basic check-up, could be improved and deepened
mouse_descriptions_voc <- assay_mouse__description |> pull(description) |> str_c(collapse = " ") |> str_split(" ") |> unlist()
mouse_descriptions_voc <- tibble(word = mouse_descriptions_voc) |> group_by(word) |> summarise(n = n(), .groups = "drop") |> arrange(desc(n)) # -> 457 unique words
# NB radiation is mentioned in some descriptionss, thus, some of the tested compounds could be evaluated as protectors, actually.
# \\b - word's boundary
# Words, which are probably describe the route of administration
iv_pattern <- "(\\biv\\b)|(\\bintravenous\\b)|(i\\.v\\.)|(\\bintravenal\\b)|(\\bintravenously\\b)|(\\bintra venous\\b)"
gastro_pattern <- "(p\\.o\\.)|(\\bintragastrically\\b)|(\\bpo\\b)|(\\bperoral\\b)|(\\boral\\b)|(\\borally\\b)|(\\bp\\.o\\.\\b)|(\\bperorally\\b)|(\\bperorla\\b)|(\\b\\(po\\)\\b)|(\\bOral\\b)|(\\bmice\\(po\\)\\b)"
subcutaneous_pattern <- "(\\bsc\\b)|(\\bsubcutaneously\\b)|(\\bsubcutaneous\\b)"
intraperitoneal_pattern <- "(\\bintraperitoneal\\b)|(\\bintraperitoneally\\b)|(\\bip\\b)|(\\b\\(ip\\)\\b)|(\\bip\\)\\b)|(\\bmice\\(ip\\)\\b)|(i\\.p\\.)"
assay_mouse__checked <- assay_mouse__description |> mutate(iv = if_else(str_detect(description, iv_pattern), 1, 0),
                                                            gi = if_else(str_detect(description, gastro_pattern), 1, 0),
                                                            sc = if_else(str_detect(description, subcutaneous_pattern), 1, 0),
                                                            ip = if_else(str_detect(description, intraperitoneal_pattern), 1, 0)) |>
                                                    rowwise() |>
                                                    mutate(unmarked = if_else(sum(iv + gi + sc + ip)>0, 0, 1)) |>
                                                    ungroup()
# Cells, very basic check-up, could be improved and deepened
cell_descriptions_voc <- assay_cell__description |> pull(description) |> str_c(collapse = " ") |> str_split(" ") |> unlist()
cell_descriptions_voc <- tibble(word = cell_descriptions_voc) |> group_by(word) |> summarise(n = n(), .groups = "drop") |> arrange(desc(n)) # -> 324 unique words
hrs12_pattern <- "12 hrs"
hrs20_pattern <- "20 hrs"
hrs24_pattern <- "24 hrs"
hrs30_pattern <- "30 hrs"
hrs44_pattern <- "44 hrs"
hrs48_pattern <- "48 hrs"
hrs72_pattern <- "(72 hrs)|(3 days)"
hrs96_pattern <- "(96 hrs)|(4 days)"
hrs120_pattern <- "5 days"
hrs144_pattern <- "6 days"
hrs168_pattern <- "7 days"
assay_cell__checked <- assay_cell__description |> mutate(hrs12 = if_else(str_detect(description, hrs12_pattern), 1, 0),
                                                            hrs20 = if_else(str_detect(description, hrs20_pattern), 1, 0),
                                                            hrs24 = if_else(str_detect(description, hrs24_pattern), 1, 0),
                                                            hrs30 = if_else(str_detect(description, hrs30_pattern), 1, 0),
                                                            hrs44 = if_else(str_detect(description, hrs44_pattern), 1, 0),
                                                            hrs48 = if_else(str_detect(description, hrs48_pattern), 1, 0),
                                                            hrs72 = if_else(str_detect(description, hrs72_pattern), 1, 0),
                                                            hrs96 = if_else(str_detect(description, hrs96_pattern), 1, 0),
                                                            hrs120 = if_else(str_detect(description, hrs120_pattern), 1, 0),
                                                            hrs144 = if_else(str_detect(description, hrs144_pattern), 1, 0),
                                                            hrs168 = if_else(str_detect(description, hrs168_pattern), 1, 0)) |>
                                                    rowwise() |>
                                                    mutate(unmarked = if_else(sum(hrs12 + hrs20 + hrs24 + hrs30 + hrs44 + hrs48 + hrs72 + hrs96 + hrs120 + hrs144 + hrs168)>0, 0, 1)) |>
                                                    ungroup()
# Summarize the discovered assay parameters
# Mouse, by assay
assay_mouse_route <- assay_mouse__checked |> pivot_longer(cols = 4:8, names_to = "category", values_to = "mark") |>
                                              group_by(category) |>
                                              summarize(count_assays = sum(mark), .groups = "drop")
# Mouse, by compound
assay_mouse_route_compound <- assay_mouse__checked |> inner_join(assay_mouse) |>
                                                      pivot_longer(cols = 4:8, names_to = "category", values_to = "mark") |>
                                                      group_by(category) |>
                                                      summarize(count_compounds = sum(mark), .groups = "drop")
# Mouse
assay_mouse_route_summary <- assay_mouse_route |> inner_join(assay_mouse_route_compound) |>
                                                  rename(route = category) |>
                                                  mutate(pct_assays = (count_assays / assay_mouse_route |> pull(count_assays) |> sum() * 100) |> round(digits = 1)) |>
                                                  mutate(pct_compounds = (count_compounds / assay_mouse_route_compound |> pull(count_compounds) |> sum() * 100) |> round(digits = 1)) |>
                                                  select(route, count_assays, count_compounds, pct_assays, pct_compounds)

# Cell, by assay
assay_cell_exposure <- assay_cell__checked |> pivot_longer(cols = 4:15, names_to = "category", values_to = "mark") |>
                                              group_by(category) |>
                                              summarize(count_assays = sum(mark), .groups = "drop")
# Cell, by compound
assay_cell_exposure_compound <- assay_cell__checked |> inner_join(assay_hcl) |>
                                                      pivot_longer(cols = 4:15, names_to = "category", values_to = "mark") |>
                                                      group_by(category) |>
                                                      summarize(count_compounds = sum(mark), .groups = "drop")
# Cell
assay_cell_exposure_summary <- assay_cell_exposure |> inner_join(assay_cell_exposure_compound) |>
                                                      rename(exposure = category) |>
                                                      mutate(pct_assays = (count_assays / assay_cell_exposure |> pull(count_assays) |> sum() * 100) |> round(digits = 1)) |>
                                                      mutate(pct_compounds = (count_compounds / assay_cell_exposure_compound |> pull(count_compounds) |> sum() * 100) |> round(digits = 1)) |>
                                                      select(exposure, count_assays, count_compounds, pct_assays, pct_compounds)

# Export
write_tsv(assay_mouse_route_summary, ".../output/mouse_tox_routes.tab")
write_tsv(assay_cell_exposure_summary, ".../output/cells_tox_exposures.tab")

# Close the connection
dbDisconnect(con)