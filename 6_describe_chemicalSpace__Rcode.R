library(tidyverse)
library(rcdk)
library(fingerprint)
library(umap)
library(praznik)
library(mirt)

# NB, the same procedure could be applied to the other types of molecular descriptors.
# In this work we will also use MNA-descriptors, their generator is proprietary and will not be included.
# It's important that MNA descriptors are generated without dictionary, thus the results obtained using them could differ significantly from those obtained using fixed-length (166-bit) MACCS keys

# Input chemical structures (SMILES) and gather the related data
cs_hcl_mouse <- read_tsv(".../output/mouse-vs-cl_ldcc50.tab") |>
						group_by(canonical_smiles) |>
						mutate(description = if_else(description == "cell clinical_candidate", "cell__clinical_candidate", description)) |>
						mutate(description = if_else(description == "mouse clinical_candidate", "mouse__clinical_candidate", description)) |>
						mutate(molregnos = unique(molregno) |> sort() |> str_c(collapse = ", ") |> as.factor()) |>
						mutate(entities = str_c(description, collapse = ", ")) |>
						ungroup() |>
						distinct() |>
						rowwise() |>
						mutate(entity = str_split(entities, ", ") |> unlist() |> unique() |> sort() |> str_c(collapse = "|")) |>
						ungroup() |>
						select(molregnos, canonical_smiles, entity)
# Parse SMILES to mols using CDK for R, it is possible to use RDkit and other libs via some means.
# This is not tidy.
mols <- parse.smiles(cs_hcl_mouse[,2] |> pull())
# Get MACCS (interpretable) fingerprints (binary feature vectors) for mols
maccs <- lapply(mols, get.fingerprint, type='maccs') |> fp.to.matrix()
# Add the descriptors (features) to the data
data <- cs_hcl_mouse |> bind_cols(maccs) |> rename_with(~ paste0('maccs_', seq_along(.x)), starts_with('...')) |>
											mutate_at(vars(starts_with("maccs_")), as.integer)
# Select the subset of 15 the most distinguishing descriptors; 11 here - is just a number of descriptors, which is easy to handle in the future work; larger number -> in general better results;
# but to some extent, which is hard to define.
n_descriptors <- 15
data_select <- MRMR( data |> select(starts_with("maccs")), data |> pull(molregnos), n_descriptors)
data_select_names  <- names(data_select$selection)
data_selection <- data |> select(molregnos, canonical_smiles, entity, any_of(data_select_names)) # Here is the data along with the subset of original descriptors, which allows distinguishing between the compounds well according to MRMR results using limited resources, since this subset is kinda small
# Check if some of the cell-tested and mouse-tested compounds have the same descriptions now
data_2vis_raw <- data_selection |> unite(maccs_id, starts_with("maccs"), sep = "-") |>
								group_by(maccs_id) |>
								mutate(category = unique(entity) |> sort() |> str_c(collapse = "|")) |>
								mutate(molregnos = unique(molregnos) |> sort() |> str_c(collapse = "|")) |>
								ungroup() |>
								select(molregnos, category, maccs_id) |>
								distinct()
# Check categories
data_2vis_check <- data_2vis_raw |> select(category) |> distinct()
# cell__clinical_candidate																								-> cell, clinical development 				-> OK
# cell 																													-> cell 									-> OK
# cell|cell__clinical_candidate 																						-> cell, clinical development 				-> OK
# cell|mouse  																											-> cell & mouse
# cell__clinical_candidate|mouse__clinical_candidate 																	-> cell & mouse, clinical development 		-> OK
# cell|mouse 																											-> cell & mouse 							-> OK
# cell|cell__clinical_candidate|mouse 																					-> cell & mouse, clinical development 		-> OK
# cell__clinical_candidate|mouse 																						-> cell & mouse, clinical development 		-> OK
# cell|cell|mouse|mouse|mouse__clinical_candidate 																		-> cell & mouse, clinical development 		-> OK
# cell__clinical_candidate|mouse__clinical_candidate|mouse 																-> cell & mouse, clinical development 		-> OK
# cell|cell__clinical_candidate|mouse__clinical_candidate 																-> cell & mouse, clinical development 		-> OK
# cell|mouse__clinical_candidate 																						-> cell & mouse, clinical development 		-> OK
# cell|mouse|mouse__clinical_candidate																					-> cell & mouse, clinical development		-> OK
# cell|cell|mouse 																										-> cell & mouse 							-> OK
# mouse 																												-> mouse 									-> OK
# mouse__clinical_candidate 																							-> mouse, clinical development 				-> OK
# mouse|mouse__clinical_candidate 																						-> mouse, clinical development 				-> OK
# cell|cell__clinical_candidate|mouse|mouse__clinical_candidate       													-> cell & mouse, clinical development 		-> OK
# cell|cell|mouse|mouse 																								-> cell & mouse 							-> OK
# mouse|mouse__clinical_candidate 																						-> mouse, clinical development 				-> OK
# cell|cell__clinical_candidate|cell|mouse|mouse 																		-> cell & mouse, clinical development  		-> OK
# cell__clinical_candidate|mouse|mouse__clinical_candidate 																-> cell & mouse, clinical development  		-> OK
# cell|cell__clinical_candidate|cell__clinical_candidate|mouse__clinical_candidate|mouse|mouse__clinical_candidate		-> cell & mouse, clinical development  		-> OK
# cell|cell__clinical_candidate|cell__clinical_candidate|mouse__clinical_candidate|mouse|mouse__clinical_candidate		-> cell & mouse, clinical development  		-> OK
# cell|cell__clinical_candidate|cell|mouse|mouse|mouse__clinical_candidate												-> cell & mouse, clinical development  		-> OK
# cell|cell__clinical_candidate|mouse__clinical_candidate|mouse 														-> cell & mouse, clinical development  		-> OK
# cell__clinical_candidate|mouse__clinical_candidate|mouse|mouse__clinical_candidate 									-> cell & mouse, clinical development  		-> OK

data_2vis <- data_2vis_raw |> rowwise() |> mutate(category_old=category,
								category = case_when(
								category == "cell__clinical_candidate" | category == "cell|cell__clinical_candidate" ~ "cell, clinical development",
								category == "cell"  ~ "cell",
								category == "cell__clinical_candidate|mouse__clinical_candidate" |
									category == "cell|cell__clinical_candidate|mouse" |
									category == "cell__clinical_candidate|mouse" |
									category == "cell|cell|mouse|mouse|mouse__clinical_candidate" | 
									category == "cell__clinical_candidate|mouse__clinical_candidate|mouse" | 
									category == "cell|cell__clinical_candidate|mouse__clinical_candidate" |
									category == "cell__clinical_candidate|mouse|mouse__clinical_candidate" | 
									category == "cell|mouse__clinical_candidate" | 
									category == "cell|mouse|mouse__clinical_candidate" ~ "cell & mouse, clinical development",
									category == "cell__clinical_candidate|mouse__clinical_candidate|mouse|mouse__clinical_candidate" ~ "cell & mouse, clinical development",
									category == "cell|cell__clinical_candidate|mouse__clinical_candidate|mouse" ~ "cell & mouse, clinical development",
									category == "cell|mouse" ~ "cell & mouse",
									category == "cell|cell|mouse|mouse" ~ "cell & mouse",
									category == "cell|cell|mouse" ~ "cell & mouse",
									category == "cell|cell__clinical_candidate|cell__clinical_candidate|mouse__clinical_candidate|mouse|mouse__clinical_candidate" ~ "cell & mouse, clinical development",
									category == "cell|cell__clinical_candidate|cell|mouse|mouse|mouse__clinical_candidate" ~ "cell & mouse, clinical development",
									category == "cell|cell__clinical_candidate|mouse|mouse__clinical_candidate" ~ "cell & mouse, clinical development",
									category == "cell|cell__clinical_candidate|cell|mouse|mouse" ~ "cell & mouse, clinical development",
									category == "mouse" ~ "mouse",
									category == "mouse__clinical_candidate" ~ "mouse, clinical development",
									category == "mouse|mouse__clinical_candidate" ~ "mouse, clinical development",
									category == "mouse|mouse__clinical_candidate" ~ "mouse, clinical development",
									.default = as.character("What?")
								)) |>
								ungroup() |>
								mutate(category = as.factor(category)) |>
								mutate(maccs_vec = maccs_id) |>
										separate_wider_delim(maccs_vec, delim="-", names=c("maccs__1", "maccs__2", "maccs__3", "maccs__4", "maccs__5", "maccs__6",
										"maccs__7","maccs__8","maccs__9","maccs__10","maccs__11","maccs__12","maccs__13","maccs__14","maccs__15"))
# Generate background points to improve the preservation of global data structure during UMAP-based dimensionality reduction, see: "СОХРАНЕНИЕ ЛОКАЛЬНОЙ И ГЛОБАЛЬНОЙ СТРУКТУРЫ ДАННЫХ ПРИ СНИЖЕНИИ РАЗМЕРНОСТИ НА ПРИМЕРЕ АЛГОРИТМА UMAP" in "Новые горизонты прикладной математики", doi:10.20948/ngpm, https://keldysh.ru/ngpm/
# More detailed description of this idea arised during RSF 23-73-01058 is currently under review in Data Intelligence, https://www.sciengine.com/DI/home
# Associated code is provided, SEE: https://github.com/RSF-23-73-01058/test_UMAP_for_MAP
featureSpace <- thetaComb(theta = c(0,1), nfact = n_descriptors, intercept = FALSE) |> as_tibble(.name_repair = "universal") |> unite(maccs_id, starts_with("..."), sep = "-") |>
					mutate(molregnos = str_c(c("sample_", row_number()), collapse = ""),
							canonical_smiles = NA,
							entity = "yet_unknown") |>
					mutate(category = as.factor("yet_unknown")) |>
					select(molregnos, canonical_smiles, entity, maccs_id, category)
featureSpace_sample <- featureSpace |> anti_join(data_2vis, by = "maccs_id") |> slice_sample(prop = .05) |>
										mutate(maccs_vec = maccs_id, order = 7) |>
										separate_wider_delim(maccs_vec, delim="-", names=c("maccs__1", "maccs__2", "maccs__3", "maccs__4", "maccs__5", "maccs__6",
										"maccs__7","maccs__8","maccs__9","maccs__10","maccs__11","maccs__12","maccs__13","maccs__14","maccs__15"))
data_2vis_ready <- bind_rows(data_2vis, featureSpace_sample) |> mutate_at(vars(starts_with("maccs__")), as.integer)
# UMAP the data and plot the results
set.seed(87)
data_umap <- umap( data_2vis_ready |> select(starts_with("maccs__")) |> as.data.frame(), preserve.seed = TRUE, metric = "manhattan", n_neighbors = n_descriptors  )
umap_coordinates <- data_umap$layout |> bind_cols(data_2vis_ready |>
										select(category), data_2vis_ready |>
										select(molregnos)) |>
										rename(UMAP_1 = `...1`, UMAP_2 = `...2`) |>
										mutate(order = case_when(
											category == "cell & mouse, clinical development" ~ 7, 
											category == "cell, clinical development" ~ 6,
											category == "mouse, clinical development" ~ 5, 
											category == "cell & mouse" ~ 4, 
											category == "cell" ~ 3,
											category == "mouse" ~ 2,
											category == "yet_unknown" ~ 1, 
											.default = 0
										)) |>
										mutate(category = fct_reorder(category, order)) |>
										arrange(order)
umap_plot  <- ggplot(umap_coordinates, aes(x = UMAP_1, y = UMAP_2)) +
															 geom_point(aes(color = category, fill = category, alpha = category), shape=21) +
															 scale_alpha_manual(values = c(.7, .8, .8, .8, .8, .8, .8, 1)) +
															 scale_color_manual(values = c('grey80', '#765D69', '#8FB9A8', '#D9402AFF', 'black', 'black', 'black')) +
															 scale_fill_manual(values = c('grey80', '#765D69', '#8FB9A8', '#D9402AFF', '#765D69', '#8FB9A8', '#D9402AFF')) +
															 theme_classic() +
															 theme(legend.position = "right") +
															 coord_fixed()
umap_plot

# Export the results
ggsave(".../output/umap_plot.png", plot = umap_plot, width = 7, height = 5, units = "in", dpi = 300)