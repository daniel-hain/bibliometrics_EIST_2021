############################################################################
# Preamble
############################################################################
### Generic preamble
rm(list=ls())

options(scipen = 5) # To deactivate annoying scientific number notation
set.seed(1337) # To have a seed defined for reproducability

# own functions
source("functions/functions_basic.R")
source("functions/00_parameters.R")

### Load standard packages
library(tidyverse)
library(magrittr)

### Extra packages
library(bibliometrix)
library(tidygraph)

############################################################################
# Select Seed articels
############################################################################

# Load bibliographic data
M <- convert2df(file = 'data/EIST_scopus.csv', dbsource = "scopus", format = "csv")


# Restrict variables
M  %<>% select(-Molecular.Sequence.Numbers, -Chemicals.CAS, -Tradenames, -Manufacturers, -Sponsors, 
               -Conference.name, -Conference.date, -Conference.code, -Conference.location,
               -CODEN, -PubMed.ID, -LA, -Publication.Stage)

# Extract Meta Tags #TODO: Maybe more?
M %<>% metaTagExtraction(Field = "AU_CO", aff.disamb = TRUE, sep = ";")
M %<>% metaTagExtraction(Field = "AU1_CO", aff.disamb = TRUE, sep = ";")
M %<>% metaTagExtraction(Field = "AU1_UN", aff.disamb = TRUE, sep = ";")
M %<>% metaTagExtraction(Field = "SR", aff.disamb = TRUE, sep = ";")
M %<>% metaTagExtraction(Field = "CR_AU", aff.disamb = TRUE, sep = ";")
M %<>% metaTagExtraction(Field = "CR_SO", aff.disamb = TRUE, sep = ";")
# create label
M %<>% rownames_to_column('XX') %>% 
  mutate(XX = paste(str_extract(XX, pattern = ".*\\d{4}"), str_sub(TI, 1,25)) %>% str_replace_all("[^[:alnum:]]", " ") %>% str_squish() %>% str_replace_all(" ", "_") %>% make.unique(sep='_'))  

rownames(M) <- M$XX

# Save 
#€M %<>% column_to_rownames('XX')
M %>% saveRDS("../temp/M.RDS")

############################################################################
# Networks Bibliographic
############################################################################
# M <- readRDS("../temp/M.RDS")

mat_bib <- M  %>% biblioNetwork(analysis = "coupling", network = "references", sep = ";", shortlabel =  FALSE)
# mat_bib %>% saveRDS("../temp/mat_bib.RDS")
# mat_bib <- readRDS("../temp/mat_bib.RDS")

g_bib <- mat_bib %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) # %N>% left_join(M %>% select(XX, SR, PY, TC, J9), by = c("name" = "XX")) %>% mutate(id = 1:n()) 

g_bib_save <- g_bib

# Restrict the network
g_bib <- g_bib_save %E>% 
  filter(weight >= cutof_edge_bib) %N>%
  # filter(percent_rank(weight) >= cutof_edge_pct_bib) %N>%
  filter(!node_is_isolated()) %N>%
  mutate(dgr = centrality_degree(weights = weight)) 
  # %N>% filter(dgr >= cutof_node_bib) 
  # %N>% filter(percent_rank(dgr) >= cutof_node_pct_bib)

# # Inspect
g_bib %N>% mutate(dgr = centrality_degree(weights = weight)) %>% as_tibble() %>% skimr::skim()
g_bib %E>% as_tibble() %>% skimr::skim()

# Community Detection
g_bib <- g_bib %N>%
  mutate(com = group_louvain(weights = weight)) %>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight)) %N>%
  unmorph()


# Community size restriction
g_bib <- g_bib %N>%
  add_count(com, name = 'com_n') %>%
  mutate(com = ifelse(com_n >= com_size_bib, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_bib <- g_bib %N>%
  filter(!is.na(com))

# Update degree
g_bib <- g_bib %N>%
  mutate(dgr = centrality_degree(weights = weight))

# Save the objects we need lateron
g_bib %>% saveRDS("../temp/g_bib.RDS")

### Merge with main data
M_bib <- M %>% select(XX) %>% inner_join(g_bib %N>% as_tibble() %>% select(name, dgr, com, dgr_int), by = c('XX' = 'name')) %>%
  distinct(XX, .keep_all = TRUE) 

# Save and remove
M_bib %>% saveRDS("../temp/M_bib.RDS")

# ### Restrict original M
# M %<>%
#   semi_join(M_bib, by = 'XX')

### Aggregated Network
require(RNewsflow)
g_bib_agg <- g_bib %>%
  network_aggregate(by = "com", edge_attribute = "weight", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight) %>%
  select(from, to, weight)

# g_bib_agg <- g_bib_agg %E>%
#   rename(weight_count = weight) %>%
#   mutate(weight = weight_count / (.N()$N[from] * .N()$N[to]) ) %>%
#   mutate(weight = (weight * 100) %>% round(4)) %N>%
#   mutate(dgr = centrality_degree(weights = weight))

g_bib_agg %>% saveRDS("../temp/g_bib_agg.RDS")

# Delete all we dont need
rm(mat_bib, g_bib, com_size_bib, cutof_edge_bib, cutof_node_bib, g_bib_agg)

############################################################################
# Network Cocitation
############################################################################

mat_cit <- M %>%
  as.data.frame() %>% 
  biblioNetwork(analysis = "co-citation", network = "references", sep = ";", shortlabel = TRUE)

# mat_cit %>% saveRDS("../temp/mat_cit.RDS")
# mat_cit <- readRDS("../temp/mat_cit.RDS")

g_cit <- mat_cit %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) # %N>% left_join(M %>% select(XX, SR, PY, TC, J9), by = c("name" = "XX")) %>% mutate(id = 1:n()) 

# Restrict the network
g_cit <- g_cit %E>% 
  filter(weight >= cutof_edge_cit) %N>%
  # filter(percent_rank(weight) >= cutof_edge_pct_cit) %N>%
  filter(!node_is_isolated()) %N>%
  mutate(dgr = centrality_degree(weights = weight)) %N>%
  filter(dgr >= cutof_node_cit) #%>%
  #filter(percent_rank(dgr) >= cutof_node_pct_cit)

# Weighting # NOTE: Only in Cicit network
g_cit <- g_cit %E>%
  mutate(weight_jac = weight / (.N()$dgr[from] + .N()$dgr[to] - weight) ) %N>%
  mutate(dgr_jac = centrality_degree(weights = weight_jac))

# # Inspect
# g_cit %N>% mutate(dgr = centrality_degree(weights = weight)) %>% as_tibble() %>% skimr::skim()
# g_cit %E>% as_tibble() %>% skimr::skim()

# Community Detection
g_cit <- g_cit %N>%
  mutate(com = group_louvain(weights = weight)) %N>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight)) %>%
  unmorph()

# Community size restriction
g_cit <- g_cit %N>%
  add_count(com, name = 'com_n') %>%
  mutate(com = ifelse(com_n >= com_size_cit, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_cit <- g_cit %N>%
  filter(!is.na(com))

# Update degree
g_cit <- g_cit %N>%
  mutate(dgr = centrality_degree(weights = weight),
         dgr_jac = centrality_degree(weights = weight_jac))

# Save the objects we need lateron
g_cit %>% saveRDS("../temp/g_cit.RDS")

# generate citation report
C_nw <- g_cit %N>%
  as_tibble()

C_nw %>% saveRDS("../temp/C_nw.RDS")

###A ggregated Network
require(RNewsflow)
g_cit_agg <- g_cit %>%
  network_aggregate(by = "com", edge_attribute = "weight", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight) %>%
  select(from, to, weight)

saveRDS(g_cit_agg, "../temp/g_cit_agg.RDS")
rm(mat_cit, g_cit, g_cit_agg)

#### 2 mode network 
m_2m <- M %>% as.data.frame() %>% cocMatrix(Field = "CR", sep = ";")

g_2m <- m_2m %>% igraph::graph_from_incidence_matrix(directed = TRUE) %>% 
  igraph::simplify() 

el_2m <- g_2m %>%
  get.edgelist() %>%
  as_tibble() %>%
  rename(from = V1,
         to = V2)

el_2m %<>%
  inner_join(M_bib %>% select(XX, com), by = c('from' = 'XX')) %>%
  rename(com_bib = com) %>%
  inner_join(C_nw %>% select(name, com), by = c('to' = 'name')) %>%
  rename(com_cit = com)

saveRDS(el_2m, "../temp/el_2m.RDS")
rm(m_2m, g_2m, el_2m)

############################################################################
# Locan citations
############################################################################

CR <- M %>% citations(field = "article", sep = ";")
CR %>% saveRDS("../temp/CR.RDS")

CRL <- M %>% localCitations(sep = ";") # For some reason takes forever...
CRL %>% saveRDS("../temp/CRL.RDS")

rm(CR, CRL)

############################################################################
# Historical citation
############################################################################

# Create a historical citation network
histResults <- M %>% histNetwork(sep = ";")
histResults %>% histPlot(n = 50, size = 10, labelsize = 5)
histResults %>% saveRDS("../temp/histResults.RDS")

rm(histResults)

############################################################################
# Conceptual Structure
############################################################################

#CS <- M %>% conceptualStructure(field="ID", method="CA", minDegree=4, clust=5, stemming=FALSE, labelsize=10, documents=10)
#CS %>% saveRDS("../temp/CS.RDS")
#rm(CS)

############################################################################
# Other stuff
############################################################################
M_threefield <- M %>% as.data.frame() %>% threeFieldsPlot(fields = c("AU", "DE", "CR_SO"), n = c(20, 20, 10))
M_threefield

M_threefield %>% saveRDS("../temp/M_threefield.RDS")
rm(M_threefield)

# M %>% authorProdOverTime(k = 10, graph = TRUE)
# M %>% rpys(sep = ";", graph = T)
# M %>% thematicMap()
# M_them_evo <- M %>% thematicEvolution(years = c(2000, 2019))

############################################################################
# Topicmodel
############################################################################

library(tidytext)
library(topicmodels)
library(textstem)

# Extract text to work with
text_tidy <- M %>% 
  as_tibble() %>%
  select(XX, AB) %>%
  rename(document = XX) 

# Some initial cleaning
text_tidy %<>% 
  mutate(AB = AB %>% 
           str_to_lower() %>%
           str_replace_all("&", "-and-") %>%
           str_remove_all("/(&trade;|&reg;|&copy;|&#8482;|&#174;|&#169;)/.*") %>%
           iconv(to = "UTF-8", sub = "byte") %>%
           str_remove_all("�.*") %>%
           str_remove_all('[:digit:]') %>%
           str_squish() 
  )  %>%
  drop_na() 

# Unnesting
#text_tidy %<>% 
#  unnest_tokens(term, AB) 

text_tidy %<>% 
  unnest_ngrams(term, AB, n_min = 1, n = 3)

# filtering
#text_tidy %<>%
#  filter(str_length(term) > 2)

text_tidy %<>% 
  separate(term, c("word1", "word2", "word3"), sep = " ")

# Stopwords
stop_words_own <- tibble(
  word =c("rights","reserved" , "study", "studies", "this", "paper", "result", "model", "approach", "article", "author", "method", "understand", "focus", "examine", "aim", "argue", "identify",
                   "increase", "datum", "potential", "explore", "include", "issue", "propose", "address", "apply", "require", "analyse", "relate", "finding",
                   "analyze", "discuss", "contribute", "publish", "involve", "draw", "lead", "exist", "set", "reduce", "create", "form", "explain", "play",
                   "affect", "regard", "associate", "establish", "follow", "conclude", "define", "strong", "attempt", "finally", "elsevier", "offer",
                   "taylor", "francis", "copyright", "springer", "wiley", "emerald", "copyright", "b.v"),
  lexicon = 'own') %>% 
  bind_rows(stop_words)

text_tidy %<>%
  filter(!word1 %in% stop_words_own$word,
         !word2 %in% stop_words_own$word,
         !word3 %in% stop_words_own$word,
         is.na(word1) | str_length(word1) > 2,
         is.na(word2) | str_length(word2) > 2,
         is.na(word3) | str_length(word3) > 2)

# Lemmatizing 
lemma_own <- tibble( # WORK IN THAT !!!!!!!!!!
  token = c("institutional", "technological", "national", "regional", "sustainable",    "environmental", "political"),
  lemma = c("institution",   "technology",    "nation",   "region",   "sustainability", "environment",   "policy"))
    
lemma_new <- lexicon::hash_lemmas %>% 
  filter(token != 'data') %>%
  anti_join(lemma_own, by = 'token') %>%
  bind_rows(lemma_own)

text_tidy %<>%
  mutate(word1 = word1 %>% lemmatize_words(dictionary = lemma_new),
         word2 = word2 %>% lemmatize_words(dictionary = lemma_new),
         word3 = word3 %>% lemmatize_words(dictionary = lemma_new))

text_tidy %<>%
  unite(term, word1, word2, word3, na.rm = TRUE, sep = " ")

# TFIDF weighting
text_tidy %<>%
  count(document, term) %>%
  bind_tf_idf(term, document, n)

text_tidy %>% saveRDS("../temp/text_tidy.RDS")

rm(stop_words_own, lemma_own, lemma_new)

# TTM
text_dtm <- text_tidy %>%
  cast_dtm(document, term, n) 
# %>% tm::removeSparseTerms(sparse = .99)

# Finding nummer of topics
library("ldatuning")

find_topics <- text_dtm %>% 
  FindTopicsNumber(
    topics = seq(from = 2, to = 15, by = 1),
    metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
    method = "Gibbs",
    control = list(seed = 77),
    mc.cores = 4L,
    verbose = TRUE
)

find_topics %>% FindTopicsNumber_plot() # Taking 6 topics

# LDA
text_lda <- text_dtm %>% LDA(k = 6, method= "Gibbs", control = list(seed = 1337))
text_lda %>% saveRDS("../temp/text_lda.RDS")

### LDA Viz
library(LDAvis)
json_lda <- topicmodels_json_ldavis(fitted = text_lda, 
                                    doc_dtm = text_dtm, 
                                    method = "TSNE")
json_lda %>% serVis()
json_lda %>% serVis(out.dir = 'output/LDAviz')

# clean up
rm(text_tidy, text_dtm, text_lda)


############################################################################
# Other network levels
############################################################################
# mat_bib %<>% normalizeSimilarity(type = "association") # NOTE: We do not normalize on the biblio-network publication level anymore.

