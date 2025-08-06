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

# 50594 - mouse
# Get the IDs for mouse LD50
mouse_LD50_cmpnd_id__query <- dbSendQuery(con, 'SELECT cs.molregno, cs.canonical_smiles FROM compound_structures cs
													JOIN molecule_hierarchy mh
													JOIN activities act
													JOIN assays a
													WHERE cs.molregno = mh.parent_molregno AND
														  act.molregno = mh.molregno AND
														  act.standard_type = "LD50" AND
														  act.standard_relation = "=" AND
														  a.assay_id = act.assay_id AND
														  a.assay_type = "T" AND
														  a.tid = "50594"')
mouse_LD50_cmpnd_id__result <- dbFetch(mouse_LD50_cmpnd_id__query) |> distinct() |>
									mutate(entity = 'mouse')
dbClearResult(mouse_LD50_cmpnd_id__query)
# Get the IDs for cell-line toxicity
hcl_tox_cmpnd_id__query <- dbSendQuery(con, 'SELECT cs.molregno, cs.canonical_smiles FROM compound_structures cs
													JOIN molecule_hierarchy mh
													JOIN activities act
													JOIN assays a
													JOIN target_dictionary td
													WHERE cs.molregno = mh.parent_molregno AND
														  act.molregno = mh.molregno AND
														  act.standard_type = "CC50" AND
														  act.standard_relation = "=" AND
														  a.assay_id = act.assay_id AND
														  a.assay_type = "T" AND
														  td.target_type = "CELL-LINE" AND
														  td.organism = "Homo sapiens" AND
														  a.tid = td.tid')
hcl_tox_cmpnd_id__result <- dbFetch(hcl_tox_cmpnd_id__query) |> distinct() |>
									mutate(entity = 'cell')
dbClearResult(hcl_tox_cmpnd_id__query)
# Get the IDs for compounds studied as drugs
drug_cmpnd_id__query <- dbSendQuery(con, 'SELECT molregno FROM molecule_dictionary
																					WHERE max_phase IS NOT NULL')
drug_cmpnd_id__result <- dbFetch(drug_cmpnd_id__query) |> mutate(type = 'clinical_candidate')
dbClearResult(drug_cmpnd_id__query)

# Result
structures <- bind_rows(hcl_tox_cmpnd_id__result, mouse_LD50_cmpnd_id__result) |>
														left_join(drug_cmpnd_id__result) |>
														replace_na(list(type = " ")) |>
														rowwise() |>
														mutate(description = str_c(entity, type, sep = " ") |> str_trim()) |>
														ungroup() |>
														select(molregno, canonical_smiles, description)


# Save the result
write_tsv(structures, ".../output/mouse-vs-cl_ldcc50.tab")

# Close the connection
dbDisconnect(con)