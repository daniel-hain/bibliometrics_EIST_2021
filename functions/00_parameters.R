
############################################################################
#  Main article selection
############################################################################

# Select variables to keep
#vars <- c("AU", "Author.s..ID", "TI", "PY", "SO", "VL", "IS", "PP", "TC", "DI", "Affiliations", "C1", "AB", "DE", "ID", "FU", "FX",
#          "CR", "RP", "LA", "JI", "DT", "DB", "UT", "J9", "AU_UN", "AU1_UN", "AU_UN_NR", "SR_FULL", "SR")

# Initial read parameters
TC_min <- 0 # Min number of citations
n_max <- 500 # Number of articles selected
PY_min <-1990 # Start year

############################################################################
#  Network Biblio
############################################################################

# Initial Filter
cutof_edge_bib <- 2
cutof_node_bib <- 5

cutof_edge_pct_bib <- 0.00
cutof_node_pct_bib <- 0.00

# community detection
com_size_bib <- 20

############################################################################
#  Network Co-Citation
############################################################################

# Initial Filter
cutof_edge_cit <- 2
cutof_node_cit <- 5

cutof_edge_pct_cit <- 0.05
cutof_node_pct_cit <- 0.05

# community detection
com_size_cit <- 40
