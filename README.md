This repo is intended to be as data free as possible, since most of the data could be obtained programmaticaly, which is good.
However, some datasets may require manual adjustments. These datasets will be uploaded.
Namely:
- cells_completed.tab contains ChEMBL [1d] data and data from other sources. Thus, ChEMBL data license (https://creativecommons.org/licenses/by-sa/3.0/) along with the licenses of other sources mentioned in the file should be considered. Also, additional processing of the data in this file are in line with those described in [2d].



The code in this repo is highly dependent on R ecosystem [1t] in general and several readily available libraries [1t-4t].

Data References:

1d. Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of bioassay data. Nucleic Acids Research. 2019 Jan;47(D1):D930-D940. DOI: 10.1093/nar/gky1075. PMID: 30398643; PMCID: PMC6323927.

2d. Lagunin AA, Dubovskaja VI, Rudik AV, et al. CLC-Pred: A freely available web-service for in silico prediction of human cell line cytotoxicity for drug-like compounds. Plos one. 2018 ;13(1):e0191838. DOI: 10.1371/journal.pone.0191838. PMID: 29370280; PMCID: PMC5784992.


Tool References:

1t. Team, R. Core. "R language definition." Vienna, Austria: R foundation for statistical computing 3.1 (2000): 116.

2t. Wickham, Hadley, et al. "Welcome to the Tidyverse." Journal of open source software 4.43 (2019): 1686.

3t. Müller K, Ooms J, James D, DebRoy S, Wickham H, Horner J (2025). _RMariaDB: Database Interface and MariaDB Driver_. doi:10.32614/CRAN.package.RMariaDB <https://doi.org/10.32614/CRAN.package.RMariaDB>, R package version 1.3.4, <https://CRAN.R-project.org/package=RMariaDB>.

4t. R Special Interest Group on Databases (R-SIG-DB), Wickham H, Müller K (2024). _DBI: R Database Interface_. doi:10.32614/CRAN.package.DBI <https://doi.org/10.32614/CRAN.package.DBI>, R package version 1.2.3, <https://CRAN.R-project.org/package=DBI>.
