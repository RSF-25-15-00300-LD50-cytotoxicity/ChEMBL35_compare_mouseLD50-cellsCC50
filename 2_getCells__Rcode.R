library(tidyverse)
library(RMariaDB)
library(DBI)
library(networkD3)

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

# Input
hcl_assay_molregno <- read_tsv(".../output/hcl_tox_assays.tab")

# Get the available data on cells
# tid, description, tissue, CLO ID, EFO ID, Cellosaurus ID, LINCS ID, cell ontology ID
# IN ({tid_all})
assay_vec <- hcl_assay_molregno |> pull(assay) |> unique() |> str_c(collapse = "', '")
assay_vec <- str_c("'", assay_vec, "'")
assay_cell__query <- dbSendQuery(con, str_glue("SELECT a.chembl_id, a.tid, c.cell_name, c.clo_id, c.efo_id, c.cellosaurus_id, c.cl_lincs_id, c.cell_ontology_id, c.cell_description, c.cell_source_tissue
                                          FROM assays a
                                          JOIN cell_dictionary c
                                          WHERE a.chembl_id IN ({assay_vec}) AND
                                            a.cell_id = c.cell_id"))
assay_cell__result <- dbFetch(assay_cell__query) |> distinct() |> rename(assay = chembl_id)
dbClearResult(assay_cell__query)

# Summarize: molregno counts by cell
data <- hcl_assay_molregno |> inner_join(assay_cell__result)
# The data on cells are incomplete, -> save & adjust cell data manually using identifiers and data from DBs
celldata_with_gaps <- data |> select(tid, cell_name, clo_id, efo_id, cellosaurus_id, cl_lincs_id, cell_ontology_id, cell_description, cell_source_tissue) |> distinct()

# Export the results
write_tsv(data, ".../output/cell_data.tab")
write_tsv(celldata_with_gaps, ".../output/cells_raw.tab")

#### The data have to be completed with the following fields:

### Progeny of the cells
## Gather cell_source_tissue to the progenies from https://data.humancellatlas.org/:
# Adipose
# Breast
# Development (embryos, fetuses and all) -> distinct parameter
# Eye
# Genetic diversity (?) - human cells within complex tissues are highly variable, influenced by genetics, age, sex, and environment.
# Gut
# Heart and Vascular
# Immune
# Kidney
# Liver
# Lung
# Musculoskeletal
# Nervous system
# Oral and Craniofacial
# Organoid (?) - three-dimensional structures of cells that recapitulate organ development in vitro – hold tremendous potential for biomedical applications.
# Pancreas
# Reproduction
# Skin 
## -> more ore less in line

### Type of the cells corresponding to the degree of their deviation from normal human cells
## Normal human cells does not last in culture long, thus the studies are often conducted using established cultures of somehow transformed cells OR cells having limited number of divisions left.
## Cellosaurus is the basic source of this information, accessed on 10.09.2025 & 11.09.2025
# finite - normal cells having limited lifespan, including primary cells, which are finite cells, which have just been obtained from the source
# transformed - more or less normal cells, which obtained the prolonged lifespan during transformation in vitro (spontaneous, viral or induced by chemicals)
# cancer - cancer cell lines, not normal at all

### Developmental stage
## Cells obtained on the developmental stage may deviate from the normal somehow; applies only to the non-cancerous cell lines
# 1  - developmental stage
# NA - no clear indication on the developmental stage was found


### It should be that there is some degree of ambiguity in the cell's characteristics.
### For example, PRIMARY and FINITE cells derived from them may have the same `NAME` in different sources.
### PBMC is cummulative `NAME` 