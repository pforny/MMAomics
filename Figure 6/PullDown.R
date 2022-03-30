# libraries

library(ggsci)
library(data.table)
require(readxl)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(igraph)
library(ggplotify)
library(ggVennDiagram)
library(UpSetR)






#colour_palette
##blue 
#82
#129
#163
myblue=rgb(60, 107, 131, maxColorValue = 255)
#orange
#231
#129
#48
myorange=rgb(231, 129, 48, maxColorValue = 255)
#brown
#180
#100
#36
mybrown=rgb(180, 100, 36, maxColorValue = 255)
#mygrey
#212
#212
#212
mygrey=rgb(180, 180, 180, maxColorValue = 255)


mypal <- pal_aaas("default", alpha = 1)(9)
mypal2 <- pal_npg("nrc", alpha = 1)(9)



# create figures path
system("mkdir Figs/pullDown/")
fig_path <- c("Figs/pullDown/")
fig_path_pdf <- c("Figs/v8/pdf/Fig6/")
fig_path_pdf_supp <- c("Figs/v8/pdf/SuppFig7/")


##################################################
# pull down data analysis, prepare data

# values are total spectrum counts (exported from Scaffold v5)

tbl <- data.table(read_excel("Data/pullDown/Samples View Report for EV MCEE MMAA MMAB MUT VLCAD.xlsx"))
setnames(tbl, c("number", "visible", "starred", "protein_description", "name", "alternate_id", "weight", "protein_group", "aov_pval", "variance", "profile", "pulled", "taxonomy", "EV1", "EV2", "EV3", "MCEE1", "MCEE2", "MCEE3", "MMAA1", "MMAA2", "MMAA3", "MMAB1", "MMAB2", "MMAB3", "MMUT1", "MMUT2", "MMUT3", "VLCAD1", "VLCAD2", "VLCAD3", "aov_prot_x")) # aov_prot_x indicates the protein against which the ANOVA was run, e.g. MMUT against EV and VLCAD, # pulled indicates whether the aov_prot_x was up while EV and VLCAD were down
tbl <- tbl[-c(1:3), ]
tbl$aov_pval <- as.numeric(tbl$aov_pval)

# exclude non-human proteins
tbl1 <- tbl[grep("HUMAN", tbl$name), ]
tbl1[, mod_name := sub("_.+", "", tbl1$name)]
tbl1[, Entry.name := sub(" .*", "", tbl1$name)]


#exclude spectrum count columns
tbl2 <- tbl1[, -c("EV1", "EV2", "EV3", "MCEE1", "MCEE2", "MCEE3", "MMAA1", "MMAA2", "MMAA3", "MMAB1", "MMAB2", "MMAB3", "MMUT1", "MMUT2", "MMUT3", "VLCAD1", "VLCAD2", "VLCAD3")]

tbl3 <- tbl2[aov_prot_x == "MMUT", ]




#SCRAPING UniProt (to assign locations to each protein (mito vs non-mito))
####################################################################################

# tblx <- "https://www.uniprot.org/uniprot/?query=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+OGDH&sort=score" %>% 
#   read_html() %>% 
#   html_table() %>% 
#   extract2(1)
# wtbl <- read.tab("Data/pullDown/uniprot-human-filtered-organism Homo+sapiens+(Human)+[9606] .tab")
# names(tbl3)
# names(wtbl)

# #add entry to scaffold table
# wjoined_scaffold <- wtbl %>% 
#   filter(Status == "reviewed") %>%
#   mutate(gname = str_extract(`Gene.names`, "[:graph:]+")) %>%
#   select(Entry, Entry.name) %>%
#   right_join(tbl3)

# wjoined_scaffold$bioprocess <- NA
# wjoined_scaffold$location <- NA
# wjoined_scaffold$link <- NA


# base_url <- "https://www.uniprot.org/uniprot/"
# #first round
# start_time1 <- Sys.time()
# for(i in 1:nrow(wjoined_scaffold)) {
#   tryCatch({
#     #nimm base url und add identifier for gname und speichere get html > save unter page als xml
#     page <- paste0(base_url, wjoined_scaffold$Entry[i]) %>% read_html()
#     link1 <- paste0(base_url, wjoined_scaffold$Entry[i])
#     wjoined_scaffold$link[i] <- link1
#     # get location
#     wjoined_scaffold$bioprocess[i] <- page %>%
#       html_node(".databaseTable") %>%
#       html_table() %>%
#       filter(X1 == "Biological process") %>%
#       .$X2
    
#     # }, error = function(e) {
#     #   cat(glue::glue("Error in gname {wjoined_essig$gname[i]}"))
#     # })  
#     # tryCatch({  
#     wjoined_scaffold$location[i] <- page %>%
#       html_nodes(".subcell_name") %>%
#       html_text() %>%
#       paste(collapse = ", ")
#   }, error = function(e) {
#     cat(glue::glue("Error in Entry.name {wjoined_scaffold$Entry.name[i]}"))
#   })
# }
# End_time1<- Sys.time()
# Run_time1 <- End_time1-start_time1 #4.146655 mins
# Run_time1


# #same amount missing in bioprocess and location > problem in reaching site?!
# table(is.na(wjoined_scaffold$bioprocess))#145 entries for bio process missing (36%) - final
# table(is.na(wjoined_scaffold$location)) #145 entries for location missing (36%) - should decrease in next round
# table(is.na(wjoined_scaffold$link)) # all links are there


# start_time2 <- Sys.time()
# for (i in which(is.na(wjoined_scaffold$location))) {
#   tryCatch({
#     #nimm base url und add identifier for gname und speichere get html > save unter page als xml
#     page <- paste0(base_url, wjoined_scaffold$Entry[i]) %>% read_html()
    
#     # get location
#     wjoined_scaffold$location[i] <- page %>%
#       html_nodes(".subcell_name") %>%
#       html_text() %>%
#       paste(collapse = ", ")
#   }, error = function(e) {
#     cat(glue::glue("Error in Entry.name {wjoined_scaffold$Entry.name[i]}"))
#   })
# }
# end_time2 <- Sys.time()
# Run_time2 <- start_time2 -end_time2 #1.149099 hours
# Run_time2

# table(is.na(wjoined_scaffold$location)) #no locations missing - not true ~35 are missing, not all on uniprot annotation but some have GO annotation
# write.csv(wjoined_scaffold, "interimData/pullDown/scrapped_scaffold_uniprot.csv")

# file of uniprot localization generated on 5 Aug 4.15 PM
wjoined_scaffold_import <- read.csv("interimData/pullDown/scrapped_scaffold_uniprot.csv", header = T, sep=",")
class(wjoined_scaffold_import)
colnames(wjoined_scaffold_import)


#rewrite location in new column as binary (1= mitochon is contained in location)
wjoined_scaffold_dt <- as.data.table(wjoined_scaffold_import)
colnames(wjoined_scaffold_dt)
wjoined_scaffold_dt[, mitochondrial := 0]
wjoined_scaffold_dt[grep("mitoch", wjoined_scaffold_import$location, ignore.case = TRUE), mitochondrial := 1]
sum(wjoined_scaffold_dt$mitochondrial)



# prep lists for venn diagram

count_tbl <- tbl1[, c("EV1", "EV2", "EV3", "MCEE1", "MCEE2", "MCEE3", "MMAA1", "MMAA2", "MMAA3", "MMAB1", "MMAB2", "MMAB3", "MMUT1", "MMUT2", "MMUT3", "VLCAD1", "VLCAD2", "VLCAD3")]
count_tbl <- as.data.table(lapply(count_tbl, as.numeric))

sum_tbl <- count_tbl %>% 
  mutate(EV = (EV1 + EV2 + EV3)) %>% 
  mutate(MCEE = MCEE1 + MCEE2 + MCEE3) %>% 
  mutate(MMAA = MMAA1 + MMAA2 + MMAA3) %>% 
  mutate(MMAB = MMAB1 + MMAB2 + MMAB3) %>% 
  mutate(MMUT = (MMUT1 + MMUT2 + MMUT3)) %>% 
  mutate(VLCAD = (VLCAD1 + VLCAD2 + VLCAD3))
sum_tbl1 <- sum_tbl %>% 
  select(EV, MCEE, MMAA, MMAB, MMUT, VLCAD)
sum_tbl1 <- cbind(tbl1$mod_name, sum_tbl1)
names(sum_tbl1)[1] <- "mod_name"


venn_list <- list(EV = sum_tbl1[sum_tbl1$EV != 0, ]$mod_name,
	MCEE = sum_tbl1[sum_tbl1$MCEE != 0, ]$mod_name,
	MMAB = sum_tbl1[sum_tbl1$MMAB != 0, ]$mod_name,
	MMUT = sum_tbl1[sum_tbl1$MMUT != 0, ]$mod_name,
	VLCAD = sum_tbl1[sum_tbl1$VLCAD != 0, ]$mod_name
	)


# prepare data/list for upset plot
# upset_list <- list(EV = sum_tbl1[sum_tbl1$EV != 0, ]$mod_name,
# 	MCEE = sum_tbl1[sum_tbl1$MCEE != 0, ]$mod_name,
# 	MMAA = sum_tbl1[sum_tbl1$MMAA != 0, ]$mod_name,
# 	MMAB = sum_tbl1[sum_tbl1$MMAB != 0, ]$mod_name,
# 	MMUT = sum_tbl1[sum_tbl1$MMUT != 0, ]$mod_name,
# 	VLCAD = sum_tbl1[sum_tbl1$VLCAD != 0, ]$mod_name
# 	)



# using gplots package to generate intersection list # doesn't work
# itemsList <- venn(venn_list, show.plot = FALSE)
# lengths(attributes(itemsList)$intersections)

# p <- calculate.overlap(venn_list)
# p

# combs <- 
#   unlist(lapply(1:length(venn_list), 
#                 function(j) combn(names(venn_list), j, simplify = FALSE)),
#          recursive = FALSE)
# names(combs) <- sapply(combs, function(i) paste0(i, collapse = ""))
# str(combs)
# combs$1



#########################################
#########################################
#########################################
# PLOTS




#########################################
# analyse ANOVA for MMUT
# select table for ANOVA with MMUT
# tbl_mmut <- tbl2[aov_prot_x == "MMUT"]
tbl_mmut <- copy(wjoined_scaffold_dt)
tbl_mmut <- tbl_mmut[!duplicated(tbl_mmut$mod_name), ]
tbl_mmut$mod_name <- factor(tbl_mmut$mod_name, levels = tbl_mmut[order(tbl_mmut$aov_pval), ]$mod_name)

prots_interest <- tbl_mmut[mod_name == "MUTA" | mod_name == "ODO2" | mod_name == "DHE3" | mod_name == "AATM", ]
prots_interest[mod_name == "MUTA", new_name := "MMUT"]
prots_interest[mod_name == "DHE3", new_name := "GLUD1"]
prots_interest[mod_name == "ODO2", new_name := "DLST"]
prots_interest[mod_name == "AATM", new_name := "GOT2"]



pval_hist <- 
ggplot(tbl_mmut[aov_pval != 1, ], aes(x = aov_pval)) +
	geom_histogram() +
	xlab("p-value") +
	ylab("Count") +
	theme_pubr()



# how many proteins are significantly enriched in the pull-down?
dim(tbl_mmut[aov_pval < 0.05, ])[1]
# what is the localization of these significant proteins?
tbl_mmut[aov_pval < 0.05, ]$mitochondrial



pval_ranks <-
ggplot(tbl_mmut[aov_pval != 1, ], aes(x = mod_name, y = aov_pval)) +
	geom_point(alpha = 0.6, shape = 1) +
	geom_point(data = tbl_mmut[aov_pval != 1 & mitochondrial == 1, ], inherit.aes = FALSE, aes(x = mod_name, y = aov_pval, color = mypal2[4]), size = 4, alpha = 0.6) +
  geom_point(data = prots_interest, inherit.aes = FALSE, aes(x = mod_name, y = aov_pval), color = mypal2[1], size = 3) +
  geom_text_repel(data = prots_interest, aes(label = new_name), color = "black", nudge_x = 150, nudge_y = 0.08, size = 3) +
  xlab("Proteins") +
	ylab("p-value") +
	ggtitle("MMUT pull-down") +
	expand_limits(x = c(-10, length(levels(tbl_mmut[aov_pval != 1, ]$mod_name)) + 5)) +
    scale_color_identity(guide = "legend", labels = c("Mitochondrial location")) +
	theme_pubr() +
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.4, 0.9), legend.title = element_blank(), plot.title = element_text(size = 12))

# pull_pvalues <- 
# ggarrange(pval_hist, pval_ranks)

# ggsave(paste(fig_path,"pvalsPullDown.png", sep = ""), pull_pvalues, width = 6.5, height = 3)
# dev.off()

ggsave(paste(fig_path,"pvalsRankedPullDown.png", sep = ""), pval_ranks, width = 3.5, height = 3)
dev.off()
ggsave(paste(fig_path_pdf,"pvalsRankedPullDown.pdf", sep = ""), device = "pdf", pval_ranks, width = 3.5, height = 3)
dev.off()



#########################################################
# illustrate relationships with network plot

edges_tbl <- tbl2[pulled == "yes" & aov_prot_x != "MMAA", c("aov_prot_x", "mod_name", "aov_pval")]
edges_tbl[mod_name == "MUTA", ]$mod_name <- "MMUT"
edges_tbl[mod_name == "ODO2", ]$mod_name <- "DLST"
edges_tbl[mod_name == "DHE3", ]$mod_name <- "GLUD1"
edges_tbl[mod_name == "AATM", ]$mod_name <- "GOT2"

edges_tbl$aov_pval <- log(1/(edges_tbl$aov_pval*10))

vertices <- unique(c(edges_tbl$mod_name, "MMUT", "MMAB", "MCEE"))
vert_tbl <- data.table(name = vertices, flagged = "No", target = "No")
vert_tbl[name == "MMUT" | name == "MMAB" | name == "MCEE", flagged := "Yes"]
vert_tbl[name == "DLST" | name == "GLUD1" | name == "GOT2", target := "Yes"]

graph <- graph_from_data_frame(d = edges_tbl, vertices = vert_tbl, directed = FALSE)

V(graph)$color <- ifelse(V(graph)$flagged == "Yes", mypal2[6], ifelse(V(graph)$target == "Yes", mypal2[5], "white"))

degree(graph)
betweenness(graph)
get.diameter(graph)


graph1 <- simplify(graph, remove.loops = TRUE) 

network <- as.ggplot(base2grob(~plot(graph1, vertex.label.color = "black", layout = layout_with_fr, edge.width = E(graph)$aov_pval, edge.color = "black", edge.curved = 0.1, vertex.label.family = "Helvetica", vertex.label.cex = 0.6, vertex.size = 17)))
network

# vertex.size = log(degree(graph))*10

ggsave(paste(fig_path,"PullDownNetwork.png", sep = ""), network, width = 7, height = 7)
dev.off()
ggsave(paste0(fig_path_pdf,"PullDownNetwork.pdf"), network, device = "pdf", width = 7, height = 7)
dev.off()




# plot venn diagram

ggVennDiagram(venn_list, label = "count")

venn <- Venn(venn_list)
data <- process_data(venn)

venn_plt <- 
ggplot() +
	geom_sf(aes(fill = id), data = venn_region(data), alpha = 0.7, show.legend = FALSE) +  
	geom_sf(aes(color = id), size = 1, data = venn_setedge(data), alpha = 0.7, show.legend = FALSE) + 
	geom_sf_text(aes(label = name), data = venn_setlabel(data)) +  
	geom_sf_text(aes(label = count), data = venn_region(data), alpha = 1) +  
	scale_fill_manual(values = c(rep("white", 18), "grey", rep("white", 13))) +
	scale_color_brewer(palette = "Dark2") +
	theme_void()

ggsave(paste(fig_path,"VennDiag.png", sep = ""), venn_plt, width = 4, height = 4)
dev.off()
ggsave(paste0(fig_path_pdf,"VennDiag.pdf"), venn_plt, device = "pdf", width = 3, height = 3)
dev.off()




# upset plot
upset(fromList(venn_list), nsets = 5)
pdf(paste0(fig_path_pdf,"upsetPullDown.pdf"), width = 6, height = 5)
upset(fromList(venn_list), nsets = 5)
dev.off()














