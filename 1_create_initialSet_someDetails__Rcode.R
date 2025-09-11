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
mouse_LD50_cmpnd_id__query <- dbSendQuery(con, 'SELECT a.chembl_id, cs.molregno, cs.canonical_smiles, act.standard_value, act.standard_units FROM compound_structures cs
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
mouse_LD50_cmpnd_id__result <- dbFetch(mouse_LD50_cmpnd_id__query) |>
									mutate(entity = 'mouse')
dbClearResult(mouse_LD50_cmpnd_id__query)
# Select the immediatelly interpretable and sufficiently large subset for the further investigation and comparison
mouse_LD50_mgkg <- mouse_LD50_cmpnd_id__result |> filter(standard_units == "mg.kg-1")
# Assays
mouse_tox_assays <- mouse_LD50_mgkg |> select(chembl_id, molregno) |> rename(assay = chembl_id)
write_tsv(mouse_tox_assays, ".../output/mouse_tox_assays.tab")
# Delete assays from the main data body for convinience
mouse_LD50_mgkg <- mouse_LD50_mgkg |> select(-chembl_id)
# Make summary of the selected set
mouse_LD50_mgkg_summary <- mouse_LD50_mgkg |> group_by(molregno) |> summarise(n = n(), p_med_value = log10(1/median(standard_value)) )
# Summarize the whole mouse set
mouse_LD50_summary <- mouse_LD50_cmpnd_id__result |> group_by(molregno) |> summarize(n = n())
mouse_n_max <- mouse_LD50_summary |> pull(n) |> max()
mouse_n_mean <- mouse_LD50_summary |> pull(n) |> mean()
mouse_n_min <- mouse_LD50_summary |> pull(n) |> min()
# Prepare the mouse set for the visualization
mouse_LD50_cmpnd_id__result <- mouse_LD50_cmpnd_id__result |> select(-standard_value, -standard_units) |> distinct()

# Get the IDs for cell-line toxicity
hcl_tox_cmpnd_id__query <- dbSendQuery(con, 'SELECT a.chembl_id, cs.molregno, cs.canonical_smiles, act.standard_value, act.standard_units, cp.full_mwt FROM compound_structures cs
													JOIN molecule_hierarchy mh
													JOIN compound_properties cp
													JOIN activities act
													JOIN assays a
													JOIN target_dictionary td
													WHERE cs.molregno = mh.parent_molregno AND
															cs.molregno = cp.molregno AND
														  act.molregno = mh.molregno AND
														  act.standard_type = "CC50" AND
														  act.standard_relation = "=" AND
														  a.assay_id = act.assay_id AND
														  a.assay_type = "T" AND
														  td.target_type = "CELL-LINE" AND
														  td.organism = "Homo sapiens" AND
														  a.tid = td.tid')
hcl_tox_cmpnd_id__result <- dbFetch(hcl_tox_cmpnd_id__query) |>
									mutate(entity = 'cell')
dbClearResult(hcl_tox_cmpnd_id__query)
# Select the immediatelly interpretable and sufficiently large subset for the further investigation and comparison
hcl_tox_nm <- hcl_tox_cmpnd_id__result |> filter(standard_units == "nM")
# Assays
hcl_tox_assays <- hcl_tox_nm |> select(chembl_id, molregno) |> rename(assay = chembl_id)
write_tsv(hcl_tox_assays, ".../output/hcl_tox_assays.tab")
# Delete assays from the main data body for convinience
hcl_tox_nm <- hcl_tox_nm |> select(-chembl_id)
# Make summary of the selected set
hcl_tox_nm_summary <- hcl_tox_nm |> group_by(molregno) |> summarise(n = n(), p_med_value = log10(1 / ( median(standard_value * (10^-9) * full_mwt) ) ))
# Summarize the whole cell set
hcl_tox_summary <- hcl_tox_cmpnd_id__result |> group_by(molregno) |> summarize(n = n())
hcl_n_max <- hcl_tox_summary |> pull(n) |> max()
hcl_n_mean <-hcl_tox_summary |> pull(n) |> mean()
hcl_n_min <- hcl_tox_summary |> pull(n) |> min()
# Prepare the cell set for the visualization
hcl_tox_cmpnd_id__result <- hcl_tox_cmpnd_id__result |> select(-standard_value, -standard_units, -full_mwt) |> distinct()

# Intersect selected subsets of mouse and cells
hcl_mouse_tox <- hcl_tox_nm_summary |> inner_join(mouse_LD50_mgkg_summary, by = "molregno") |>
																		select(p_med_value.x, p_med_value.y) |>
																		rename(cell_tox = p_med_value.x, mouse_tox = p_med_value.y)
# Plot
ggplot(data=hcl_mouse_tox, aes(x=cell_tox, y=mouse_tox)) +
						geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
						geom_point() +
						theme_classic()
# Save the result
write_tsv(hcl_mouse_tox, ".../output/mouse-vs-hcl_intersection.tab")

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