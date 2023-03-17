###########################################################################################
########################### Preamble
###########################################################################################

### Generic preamble
rm(list=ls())
set.seed(1337)

### Load packages  
library(tidyverse)
library(magrittr)

### Extra packages
# Biblis & NWs
library(bibliometrix)
library(tidygraph)
require(RNewsflow)
#NLP
library(tidytext)
library(topicmodels)
library(textstem)
library(ldatuning)

###########################################################################################
########################### Variable definitions
###########################################################################################

rm(list=ls())
# own Parameters
source("functions/00_parameters.R")
source("functions/functions_basic.R")

# institute and department
var_inst <- 'EIST'
var_dept <- 'v2'

###########################################################################################
########################### Load & preprocessing articles
###########################################################################################

print('Starting: Loading Files')

files <- list.files(path = 'data', pattern = paste0('scopus_', str_to_lower(var_inst), '_', str_to_lower(var_dept)), full.names = TRUE)

# Load bibliographic data
M <- convert2df(file = files, dbsource = "scopus", format = "csv") %>% 
  # Delete duplicates 
  distinct(UT, .keep_all = TRUE) 

# Filter 
M %<>% 
  # Filter number references
  mutate(CR_n = CR %>% str_count(';')) %>%
  # Abstract
  filter(AB != '') %>%
  filter(AB %>% str_length() >= 25) %>%
  # Number of cited references and citations
  mutate(TC_year = TC / (PY_max + 1 - PY)) 

# create label & Rownames
M %<>% rownames_to_column('XX') %>% 
  mutate(XX = paste(str_extract(XX, pattern = ".*\\d{4}"), str_sub(TI, 1,25)) %>% str_replace_all("[^[:alnum:]]", " ") %>% str_squish() %>% str_replace_all(" ", "_") %>% make.unique(sep='_'))
  
rownames(M) <- M$XX

# Save whole compilation
M %>% saveRDS(paste0('../temp/M_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))


# M <- read_rds(paste0('../temp/M_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
M %>% biblioAnalysis(sep = ";") %>% saveRDS(paste0('../temp/M_res_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

###########################################################################################
########################### Networks Bibliographic
###########################################################################################

print('Starting: Bibliographic Coupling network')

mat_bib <- M  %>% biblioNetwork(analysis = "coupling", network = "references", sep = ";", shortlabel =  FALSE)
mat_bib %>% saveRDS(paste0('../temp/mat_bib__', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
# mat_bib <- readRDS(paste0('../temp7mat_bib__', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

g_bib <- mat_bib %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) %N>% 
  left_join(M %>% select(XX, UT, PY, CR_n, TC_year), by = c("name" = "XX"))

## Restrict the network
g_bib <- g_bib %E>% 
  filter(weight >= cutof_edge_bib)

g_bib <- g_bib %N>%
  filter(!node_is_isolated()) %N>% # Could add:  | int_dept == TRUE
  mutate(dgr = centrality_degree(weights = weight)) %N>% 
  filter(dgr >= cutof_node_bib)

# Jaccard weighting
g_bib <- g_bib %E>% 
  mutate(weight_jac = weight / (.N()$CR_n[from] + .N()$CR_n[to] - weight) ) %E>%
  mutate(weight_jac = if_else(weight_jac > 1, 1, weight_jac) ) %N>%
  mutate(dgr_jac = centrality_degree(weights = weight_jac)) 

# # Further restrictions
# g_bib <- g_bib  %N>%
#   filter(percent_rank(dgr_jac) >= cutof_node_pct_bib) %E>% 
#   filter(percent_rank(weight_jac) >= cutof_edge_pct_bib) %N>%
#   filter(!node_is_isolated())

## Community Detection
g_bib <- g_bib %N>%
  mutate(com = group_louvain(weights = weight_jac)) %>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight_jac)) %N>%
  unmorph()

# Community size restriction
com_size_bib <- (g_bib %N>% as_tibble() %>% nrow()) * 0.05
g_bib %N>% as_tibble() %>% count(com, sort = TRUE)

g_bib <- g_bib %N>%
  group_by(com) %>%
    mutate(com_n = n()) %>%
  ungroup() %>%
  mutate(com = ifelse(com_n >= com_size_bib, com, NA) ) %>%
  mutate(com = ifelse(com <= com_max_bib, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_bib <- g_bib %N>%
  filter(!is.na(com))

# Update degree
g_bib <- g_bib %N>%
  mutate(dgr = centrality_degree(weights = weight),
         dgr_jac = centrality_degree(weights = weight_jac))


# Save the objects we need lateron
g_bib %>% saveRDS(paste0('../temp/g_bib_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

## Merge with main data
M_bib <- M %>% select(UT, XX) %>% inner_join(g_bib %N>% as_tibble() %>% select(UT, dgr, dgr_jac, com, dgr_int), by = 'UT') %>%
  distinct(UT, .keep_all = TRUE) 

M_bib %>% saveRDS(paste0('../temp/M_bib_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))


## Aggregated Network
g_bib_agg <- g_bib %N>%
  filter(!is.na(com)) %>%
  network_aggregate(by = "com", edge_attribute = "weight_jac", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight_jac) %>%
  select(from, to, weight)

## Weight edges
# g_bib_agg <- g_bib_agg %E>%
#   rename(weight_count = weight) %>%
#   mutate(weight = weight_count / (.N()$N[from] * .N()$N[to]) ) %>%
#   mutate(weight = (weight * 100) %>% round(4)) %N>%
#   mutate(dgr = centrality_degree(weights = weight))

# Save the objects we need lateron
g_bib_agg %>% saveRDS(paste0('../temp/g_bib_agg_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# Delete all we dont need
rm(mat_bib, g_bib, com_size_bib, cutof_edge_bib, cutof_node_bib, g_bib_agg)




###########################################################################################
########################### Collaboration NW
###########################################################################################

pub_inst <- M %>% as_tibble() %>% metaTagExtraction(Field = "AU_UN") %>% 
  select(XX, UT, PY, AU_UN) %>% 
  # Sepperate AU_UN for 2_m edgelist
  separate_rows(AU_UN, sep = ';') %>%
  # filter
  drop_na(AU_UN) %>%
  filter(!(AU_UN %in% c('', ' ', 'NA', 'NOTREPORTED', 'NOTDECLARED'))) %>%
  # Only 1 link per paper, independent of author number
  distinct(UT, AU_UN, .keep_all = TRUE) 

el_inst <- pub_inst %>% 
  left_join(pub_inst %>% select(UT, AU_UN), by = 'UT') %>% 
  rename(from = AU_UN.x, to = AU_UN.y) %>%
  filter(from != to) %>%
  group_by(from, to) %>%
  summarise(weight = n()) %>%
  ungroup()

# Save it 
pub_inst %>% saveRDS(paste0('../temp/pub_inst_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
el_inst %>% saveRDS(paste0('../temp/el_inst_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# clean 
# rm(pub_inst, el_inst)

###########################################################################################
########################### Network Cocitation 
###########################################################################################

print('Starting: Co-Citation Network')

mat_cit <- M %>%
  semi_join(M_bib, by = 'UT') %>%
  as.data.frame() %>% 
  biblioNetwork(analysis = "co-citation", network = "references", sep = ";", shortlabel = FALSE)

mat_cit %>% saveRDS(paste0('../temp/mat_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
# mat_cit <- readRDS(paste0('../temp/mat_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

g_cit <- mat_cit %>% igraph::graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE, diag = FALSE) %>% 
  igraph::simplify() %>%
  as_tbl_graph(directed = FALSE) # %N>% left_join(M %>% select(XX, SR, PY, TC, J9), by = c("name" = "XX")) %>% mutate(id = 1:n())

# Restrict the network
g_cit <- g_cit %E>% 
  filter(weight >= cutof_edge_cit) %N>%
  filter(!node_is_isolated())

g_cit <- g_cit %N>%
  mutate(dgr = centrality_degree(weights = weight)) %N>%
  filter(dgr >= cutof_node_cit) 

# Further restrictions
g_cit <- g_cit %N>% 
  filter(percent_rank(dgr) >= cutof_node_pct_cit) %E>%
  filter(percent_rank(weight) >= cutof_edge_pct_cit) %N>%
  filter(!node_is_isolated())

## Community Detection
g_cit <- g_cit %N>%
  mutate(com = group_louvain(weights = weight)) %N>%
  morph(to_split, com) %>% 
  mutate(dgr_int = centrality_degree(weights = weight)) %>%
  unmorph()

g_cit %N>% as_tibble() %>% count(com)

# community detection
com_size_cit <- (g_cit %N>% as_tibble() %>% nrow()) * 0.05

# Community size restriction
g_cit <- g_cit %N>%
  group_by(com) %>%
  mutate(com_n = n()) %>%
  ungroup() %>%
  mutate(com = ifelse(com_n >= com_size_cit, com, NA) ) %>%
  mutate(com = ifelse(com <= com_max_cit, com, NA) ) %>%
  select(-com_n)  

# Delete nodes withou community
g_cit <- g_cit %N>%
  filter(!is.na(com))

# Update degree
g_cit <- g_cit %N>%
  mutate(dgr = centrality_degree(weights = weight))

# Save the objects we need lateron
g_cit %>% saveRDS(paste0('../temp/g_cit_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# generate citation report
C_nw <- g_cit %N>% as_tibble() 
C_nw %>%  saveRDS(paste0('../temp/C_nw_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

## A ggregated Network
require(RNewsflow)
g_cit_agg <- g_cit %>%
  network_aggregate(by = "com", edge_attribute = "weight", agg_FUN = sum)  %>%
  as.undirected(mode = "collapse", edge.attr.comb = "sum") %>%
  as_tbl_graph(directed = FALSE) %N>%
  select(-name) %>%
  mutate(id = 1:n()) %E>%
  rename(weight = agg.weight) %>%
  select(from, to, weight)

g_cit_agg %>% saveRDS(paste0('../temp/g_cit_agg_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(mat_cit, g_cit, g_cit_agg)

###########################################################################################
########################### 2 mode network 
###########################################################################################

print('Starting: 2 Mode Network')

rownames(M) <- M %>% pull(UT)

m_2m <- M %>% 
  semi_join(M_bib) %>%
  as.data.frame() %>% cocMatrix(Field = "CR", sep = ";", short = FALSE)

g_2m <- m_2m %>% igraph::graph_from_incidence_matrix(directed = TRUE, mode = 'out', multiple = FALSE) %>% 
  igraph::simplify() 

el_2m <- g_2m %>%
  get.edgelist() %>%
  as_tibble() %>%
  rename(from = V1,
         to = V2)

el_2m %<>%
  left_join(M_bib %>% select(UT, com), by = c('from' = 'UT')) %>%
  rename(com_bib = com) %>%
  left_join(M %>% select(UT, PY), by = c('from' = 'UT')) %>%
  left_join(C_nw %>% select(name, com), by = c('to' = 'name')) %>%
  rename(com_cit = com) %>% 
  drop_na(PY, com_bib, com_cit)

# save
el_2m %>% saveRDS(paste0('../temp/el_2m_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(m_2m, g_2m, el_2m, C_nw)

###########################################################################################
########################### Topicmodel
########################################################################################### 

# Extract all for externbal BERTopic Modelling
M %>% 
  as_tibble() %>%
  select(UT, PY, TI, AB) %>%
  mutate(text = paste(TI, AB, sep = '. ')  %>% 
           str_to_lower() %>% 
           str_remove_all("©.*") %>%          
           str_remove_all("/(&trade;|&reg;|&copy;|&#8482;|&#174;|&#169;)/.*") %>%
           str_remove_all('(&elsevier;|&springer;|&rights reserved)/.*') %>%
           str_squish())  %>%
  select(UT, PY, text) %>%
  write_csv('data/data_text.csv')


print('Starting: Topic Modelling')

# Extract text to work with
text_tidy <- M %>% 
  as_tibble() %>%
  select(UT, AB) %>%
  rename(document = UT,
         text = AB) 



# Some initial cleaning
text_tidy %<>% 
  mutate(text = text %>% 
           str_to_lower() %>%
           str_replace_all("&", "-and-") %>%
           str_remove_all("/(&trade;|&reg;|&copy;|&#8482;|&#174;|&#169;)/.*") %>%
           iconv(to = "UTF-8", sub = "byte") %>%
           str_remove_all("�.*") %>%
           str_remove_all('[:digit:]') %>%
           str_squish() 
  )  %>%
  drop_na() 

# n grams
text_tidy %<>% 
  unnest_ngrams(term, text, ngram_delim = ' ', n_min = 1, n = 3) %>% 
  separate(term, c("word1", "word2", "word3"), sep = " ")

# Stopwords
stop_words_own <- tibble(
  word =c("the", "rights","reserved" , "study", "studies", "these", "this", "paper", "result", "model", "approach", "article", "author", "method", "understand", "focus", "examine", "aim", "argue", "identify",
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
  token = c("systems", "institutional", "technological", "national", "regional", "sustainable",    "environmental", "political", "politic", "politics"),
  lemma = c("system", "institution",   "technology",    "nation",   "region",   "sustainability", "environment", "policy", "policy", "policy"))

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

# TTM
text_dtm <- text_tidy %>%
  cast_dtm(document, term, n) %>% tm::removeSparseTerms(sparse = .99)

# # Finding nummer of topics
# find_topics <- text_dtm %>%
#   FindTopicsNumber(
#     topics = seq(from = 5, to = 15, by = 1),
#     metrics = c("CaoJuan2009", "Arun2010", "Deveaud2014", "Griffiths2004"), # NOTE: Leaving out for now matric "Griffiths2004", due to problems installing the "gmp" package
#     method = "Gibbs",
#     control = list(seed = 1337),
#     mc.cores = 3L,
#     verbose = TRUE
#   )
# 
# find_topics %>% FindTopicsNumber_plot()
n_topic = 7

text_lda <- text_dtm %>% LDA(k = n_topic, method= "Gibbs", control = list(seed = 1337))

### LDA Viz
library(LDAvis)
json_lda <- topicmodels_json_ldavis(fitted = text_lda, 
                                    doc_dtm = text_dtm, 
                                    method = "TSNE")
#json_lda %>% serVis()

# Save
text_tidy %>% saveRDS(paste0('../temp/text_tidy_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
text_lda %>% saveRDS(paste0('../temp/text_lda_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
json_lda %>% serVis(out.dir = paste0('output/topic_modelling/LDAviz_', str_to_lower(var_inst), '_', str_to_lower(var_dept)))

###########################################################################################
########################### Similarity to past & future
###########################################################################################


# # old version with topic
# text_lda_gamma <- text_lda %>% 
#   tidy(matrix = "gamma")
# 
# el_sim_topic <- text_lda_gamma %>%
#   widyr::pairwise_similarity(document, topic, gamma, diag = FALSE, upper = TRUE)

embeddings <- read_csv('data/data_text.csv') %>%
  select(UT) %>%
  bind_cols(read_csv('data/embeddings_specter.csv')) %>%
  pivot_longer(-UT, names_to = 'dim')

el_sim <- embeddings %>%
  widyr::pairwise_similarity(UT, dim, value, diag = FALSE, upper = TRUE)

# join with year = uni
el_sim %<>%
  inner_join(M %>% select(UT, PY), by = c('item1' = 'UT')) %>%
  inner_join(M %>% select(UT, PY), by = c('item2' = 'UT')) %>%
  rename(PY_from = PY.x, PY_to = PY.y)

# decide for similarity past or future
el_sim %<>%
  mutate(delta = PY_to - PY_from) %>%
  mutate(sim_type = case_when(
    delta == 0 ~ "present",
    delta >= 1 ~ "future",
    delta <= -1 ~ "past") 
  )

# aggregate on document level
pub_sim <- el_sim %>%
  group_by(item1, sim_type) %>%
  summarise(sim = mean(similarity)) %>%
  pivot_wider(names_from = sim_type, names_prefix = 'sim_', values_from = sim) %>%
  drop_na()

uni_sim <- pub_sim %>%
  inner_join(pub_inst %>% select(UT, AU_UN), by = c('item1' = 'UT')) %>%
  group_by(AU_UN) %>%
  summarise(sim_past = mean(sim_past),
            sim_present = mean(sim_present),
            sim_future = mean(sim_future),
            n = n()) %>%
  mutate(future_trend = sim_future - sim_past)

# Save it 
uni_sim %>% saveRDS(paste0('../temp/uni_sim_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

# clean up
rm(text_tidy, text_dtm, text_lda, json_lda, uni_sim, pub_sim)

###########################################################################################
########################### Local citations
########################################################################################### 

print('Starting: Further Analysis')

CR <- M %>% citations(sep = ";")
CR %>% saveRDS(paste0('../temp/CR_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

#CRL <- M %>% localCitations(sep = ";") # For some reason takes forever...
#CRL %>% saveRDS(paste0('../temp/CRL_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))

rm(CR)

###########################################################################################
############################ Threefield Plot
########################################################################################### 

M_threefield <- M %>% as.data.frame() %>% threeFieldsPlot(fields = c("AU", "DE", "CR_SO"), n = c(20, 20, 10))
M_threefield %>% saveRDS(paste0('../temp/threefield_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
rm(M_threefield)

###########################################################################################
########################### Historical citation
########################################################################################### 

#histResults <- M %>% histNetwork(sep = ";")
#histResults %>% saveRDS(paste0('../temp/histResult_', str_to_lower(var_inst), '_', str_to_lower(var_dept), '.rds'))
#rm(histResults)

############################################################################
# Conceptual Structure
############################################################################

#CS <- M %>% conceptualStructure(field="ID", method="CA", minDegree=4, clust=5, stemming=FALSE, labelsize=10, documents=10)
#CS %>% saveRDS("../temp/CS.RDS")
#rm(CS)

############################################################################
# Other Stuff
############################################################################

#M %>% authorProdOverTime(k = 10, graph = TRUE)
#M %>% rpys(sep = ";", graph = T)
#M %>% thematicMap()
#M_them_evo <- M %>% thematicEvolution(years = c(2016, 2018,2021))

############################################################################
# Other network levels
############################################################################
# mat_bib %<>% normalizeSimilarity(type = "association") # NOTE: We do not normalize on the biblio-network publication level anymore.