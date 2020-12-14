# Title and author information --------------------------------------------
#!/usr/bin/R

#########################################
#                                       #
#     2020_11_23_adenoma_analysis.R     #
#                                       #
#########################################

#Title: Microbial signatures of colonic adenomas
#
#Copyright (C) 2020-2021  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

#Purpose: To examine differences in microbiome diversity between adenoma formers
# and non-formers

# SET ENVIRONMENT -------------------------------------------------------------
library(ggplot2)
library(vegan)
library(reshape2)
library(randomForest)
library(qvalue)
library(gplots)
library(pROC)
library(MASS)
library(glmmTMB)

options("stringsAsFactors" = F)

# FUNCTIONS ---------------------------------------------------------------

###
#        Function diversity_analysis             #
###

#this function performs alpha and beta diversity analysis on the
#user provided data frame

#depends on vegan

diversity_analysis <- function(obj) {
  obj$shannon <-
    vegan::diversity(obj$data, index = "shannon", MARGIN = 1)
  obj$simpson <-
    vegan::diversity(obj$data, index = "simpson", MARGIN = 1)
  obj$invsimpson <-
    vegan::diversity(obj$data, index = "invsimpson", MARGIN = 1)
  return(obj)
}

###
#              Function ordinate                 #
###

#This function creates ordination objects for plotting later

ordinate <- function(obj) {
  obj_mmds <- metaMDS(obj$data,
                      k = 5,
                      distance = "bray")


  obj_prcomp <- prcomp(obj$data,
                       scale = T,
                       center = T)

  obj$mds  <- obj_mmds
  obj_mmds <- as.data.frame(obj_mmds$points)
  obj_mmds <- obj_mmds[rownames(obj$meta), ]
  obj$mds.df <- obj_mmds


  obj$prcomp <- obj_prcomp
  obj_prcomp    <- as.data.frame(obj_prcomp$x[, 1:5])
  obj_prcomp    <- obj_prcomp[rownames(obj$meta), ]
  obj$prcomp.df <- obj_prcomp


  return(obj)
}

###
#             Function run_beta                  #
###

#calculates beta-div between all samples

run_beta <- function(obj) {
  bdiv <- as.matrix(vegdist(obj$data))
  obj$bdiv <- bdiv
  return(obj)
}

###
#        Function phylotype_analysis             #
###

phylotype_analysis <- function(obj, tax) {
  #obj: microbiome object with at least 1 slot (data)
  #tax: a tax object (named list taxa as names values in the list are seq ids)
  obj.out <- NULL
  for (h in 1:length(tax)) {
    df <- NULL
    for (i in 1:length(tax[[h]])) {
      v1       <- obj$data[, unlist(tax[[h]][[i]])]
      v2       <- names(tax[[h]])[i]
      if (is.null(dim(v1))) {
        df[[v2]] <- v1
      } else{
        df[[v2]] <- rowSums(v1)
      }
    }
    obj.out[[names(tax)[h]]] <- as.data.frame(df)
  }
  return(obj.out)
}

###
#            Function group_means                #
###

#This function will calculate the group means across all columns in a
#dataframe

group_means <- function(df, mapping, t = F) {
  #df : data frame with rows as sample IDs and columns as objects (e.g. OTUs)
  #mapping : a mapping df with rownames df as rownames and group id in col2
  #t : boolean, defaults to FALSE. Use if df has sample ids as colnames
  if (t == T) {
    df <- t(df)
  }
  groups <- base::unique(x = mapping[, 2])
  my_df <-
    data.frame(matrix(nrow = length(groups), ncol = ncol(df)))
  for (i in 1:ncol(df)) {
    tgvec <- NULL
    for (j in 1:length(groups)) {
      s <- base::rownames(mapping[base::which(mapping[, 2] %in% groups[j]),
                                  , drop = F])
      m <- base::mean(df[s, i])
      tgvec <- c(tgvec, m)
    }
    my_df[, i] <- tgvec
  }
  rownames(my_df) <- groups
  colnames(my_df) <- colnames(df)
  return(my_df)
}

#              Function filter_df                #

#filter data frame to remove zero sum columns and rows
filter <- function(df) {
  df <- df[which(rowSums(df) > 0), which(colSums(df) > 0)]
  return(df)
}

# IMPORT DATA -------------------------------------------------------------

polyp2_asv.df <-
  read.table(
    "/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/rarefied_10000_adenoma_asv.txt",
    sep = "\t",
    row.names = 1,
    header = T
  )
polyp2_asv_rel.df <-
  read.table(
    "/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/rel_abd_adenoma_asv.txt",
    sep = "\t",
    row.names = 1,
    header = T
  )
polyp2_tax.df <-
  read.table(
    "/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_adenoma_tax.txt",
    sep = "\t",
    row.names = 1,
    header = T
  )
polyp2_metadata.df <-
  read.table(
    "/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_combined_metadata.txt",
    sep = "\t",
    row.names = 1,
    header = T
  )

polyp2_sequence_stats.df <-
  read.table(
    "/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_adenoma_sequence_stats.txt",
    sep = "\t",
    row.names = 1,
    header = T
  )

# ANALYSIS: READ COUNT SUMMARY STATS --------------------------------------

#just some quick summary stats

sum(polyp2_sequence_stats.df$reads.in)
#[1] 37183608
mean(polyp2_sequence_stats.df$reads.in)
#[1] 70290.37
median(polyp2_sequence_stats.df$reads.in)
#[1] 68913

sum(polyp2_sequence_stats.df$reads.out)
#[1] 27523989
mean(polyp2_sequence_stats.df$reads.out)
#[1] 52030.22
median(polyp2_sequence_stats.df$reads.out)
#[1] 50882

# ANALYSIS: FILTERED READ COUNT SUMMARY STATS ---------------------------------

# remove samples filtered for rarefaction
filtered_polyp2_sequence_stats.df <-
  polyp2_sequence_stats.df[which(rownames(polyp2_sequence_stats.df) %in% rownames(polyp2_metadata.df)), ]

filtered_polyp2_sequence_stats.df$type <- polyp2_metadata.df$type

mfiltered_polyp2_sequence_stats.df <-
  melt(filtered_polyp2_sequence_stats.df)

reads_by_type.hist <-
  ggplot(data = mfiltered_polyp2_sequence_stats.df[which(mfiltered_polyp2_sequence_stats.df$variable == "reads.out"), ],
         aes(x = value,
             fill = type))

pdf("figs/qcreads_by_tissue.pdf")
reads_by_type.hist +
  geom_histogram(alpha = .5,
                 position = "identity",
                 color = "black") +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black")
  ) +
  scale_x_continuous(expand = c(.02, 0)) +
  scale_y_continuous(expand = c(.01, 0)) +
  ylab("Count") +
  xlab("Read Depth") +
  scale_fill_brewer("Tissue", palette = "Set1")
dev.off()

# Analysis: Identify Complete cases  -----------------------------------------

# Some of our planned analyses involve identifying location level shifts in the
# microbiome. However, there are several tissue sites that have few samples
# which will likely reduce the power of the analyses we are conducting. Thus,
# I will first identify the complete cases that are present in our data. For our
# purposes we will define complete cases as those that have RC, LC, and REC
# tissue. Ideally we would also include transverse but few patients actually
# had this tissue. I think that we should work with these cases for the longitudinal
# analysis.

lc.vec    <-
  polyp2_metadata.df[polyp2_metadata.df$location == "LC", "id"]
rc.vec    <-
  polyp2_metadata.df[polyp2_metadata.df$location == "RC", "id"]
rec.vec   <-
  polyp2_metadata.df[polyp2_metadata.df$location == "REC", "id"]
fecal.vec <-
  polyp2_metadata.df[polyp2_metadata.df$location == "Fecal", "id"]
oral.vec  <-
  polyp2_metadata.df[polyp2_metadata.df$location == "Oral", "id"]

#only the three gut tissues
complete_gut.cases <-
  unique(rec.vec[which(rec.vec %in% lc.vec & rec.vec %in% rc.vec)])

#three gut tissues plus fecal and oral
complete.cases <- unique(rec.vec[which(rec.vec %in% lc.vec &
                                         rec.vec %in% rc.vec &
                                         rec.vec %in% fecal.vec &
                                         rec.vec %in% oral.vec)])

# DATA: MAKE OBJECT -------------------------------------------------------

polyp2_obj <- NULL
polyp2_obj$data <- polyp2_asv.df
polyp2_obj$meta <- polyp2_metadata.df
polyp2_obj$rel_abd <- polyp2_asv_rel.df

# DATA: AGGREGATE PHYLOTYPES -----------------------------------------------------

polyp2_taxadf <-
  polyp2_tax.df[which(rownames(polyp2_tax.df) %in% colnames(polyp2_obj$data)),
                , drop = F]

kingdom.df <- replicate(length(unique(polyp2_taxadf[, 2])), c())
names(kingdom.df) <- unique(polyp2_taxadf[, 2])
phylum.df  <- replicate(length(unique(polyp2_taxadf[, 3])), c())
names(phylum.df) <- unique(polyp2_taxadf[, 3])
class.df   <- replicate(length(unique(polyp2_taxadf[, 4])), c())
names(class.df) <- unique(polyp2_taxadf[, 4])
order.df   <- replicate(length(unique(polyp2_taxadf[, 5])), c())
names(order.df) <- unique(polyp2_taxadf[, 5])
family.df  <- replicate(length(unique(polyp2_taxadf[, 6])), c())
names(family.df) <- unique(polyp2_taxadf[, 6])
genus.df   <- replicate(length(unique(polyp2_taxadf[, 7])), c())
names(genus.df) <- unique(polyp2_taxadf[, 7])

for (i in 1:nrow(polyp2_taxadf)) {
  kingdom.df[[polyp2_taxadf[i, 2]]] <-
    c(kingdom.df[[polyp2_taxadf[i, 2]]], polyp2_taxadf[i, 1])
  phylum.df[[polyp2_taxadf[i, 3]]]  <-
    c(phylum.df[[polyp2_taxadf[i, 3]]], polyp2_taxadf[i, 1])
  class.df[[polyp2_taxadf[i, 4]]]   <-
    c(class.df[[polyp2_taxadf[i, 4]]], polyp2_taxadf[i, 1])
  order.df[[polyp2_taxadf[i, 5]]]   <-
    c(order.df[[polyp2_taxadf[i, 5]]], polyp2_taxadf[i, 1])
  family.df[[polyp2_taxadf[i, 6]]]  <-
    c(family.df[[polyp2_taxadf[i, 6]]], polyp2_taxadf[i, 1])
  genus.df[[polyp2_taxadf[i, 7]]]   <-
    c(genus.df[[polyp2_taxadf[i, 7]]], polyp2_taxadf[i, 1])

}

tax.obj <- NULL
tax.obj$kingdom <- kingdom.df
tax.obj$phylum  <- phylum.df
tax.obj$class   <- class.df
tax.obj$order   <- order.df
tax.obj$family  <- family.df
tax.obj$genus   <- genus.df

#aggregate phylotype counts (not really a df, actually an obj)
polyp2_obj$phylotype <- phylotype_analysis(polyp2_obj, tax = tax.obj)

# ANALYSIS: COMPUTE DIVERSITY ---------------------------------------------

set.seed(731)
polyp2_obj <- diversity_analysis(polyp2_obj)
polyp2_obj <- run_beta(polyp2_obj)
polyp2_obj <- ordinate(polyp2_obj)

# ANALYSIS: SHANNON ------------------------------------------------------

polyp2_obj$shannon.df <- cbind(
  polyp2_obj$shannon,
  polyp2_obj$meta[names(polyp2_obj$shannon), "type"],
  polyp2_obj$meta[names(polyp2_obj$shannon), "location"],
  polyp2_obj$meta[names(polyp2_obj$shannon), "polyp"],
  polyp2_obj$meta[names(polyp2_obj$shannon), "Adenoma"],
  polyp2_obj$meta[names(polyp2_obj$shannon), "polyp.tissue"],
  polyp2_obj$meta[names(polyp2_obj$shannon), "id"]
)

colnames(polyp2_obj$shannon.df) <-
  c("shannon",
    "tissue",
    "location",
    "former",
    "npolyp",
    "polyptissue",
    "id")

polyp2_obj$shannon.df <- as.data.frame(polyp2_obj$shannon.df)
polyp2_obj$shannon.df$shannon <-
  as.numeric(polyp2_obj$shannon.df$shannon)
polyp2_obj$shannon.df$npolyp <-
  as.numeric(polyp2_obj$shannon.df$npolyp)

polyp2_shannon.boxplot <- ggplot(na.omit(polyp2_obj$shannon.df),
                                 aes(x = tissue,
                                     y = shannon,
                                     fill = former))

pdf("figs/shannon_tissue_polyp.pdf")
polyp2_shannon.boxplot +
  geom_boxplot() +
  geom_point(
    position = position_dodge(width = .75),
    color = "black",
    shape = 21,
    alpha = .5
  ) +
  theme(
    text = element_text(size = 24, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1.5,
    legend.position = "top"
  ) +
  ylab("Shannon") +
  xlab("") +
  scale_fill_brewer("",
                    palette = "Dark2",
                    labels = c("Non-Former", "Former")) +
  scale_color_brewer("", palette = "Dark2")
dev.off()

polyp2_shannon.boxplot <- ggplot(na.omit(polyp2_obj$shannon.df),
                                 aes(x = location,
                                     y = shannon,
                                     fill = former))

pdf("figs/shannon_location_boxplot.pdf")
polyp2_shannon.boxplot +
  geom_boxplot() +
  geom_point(
    position = position_dodge(width = .75),
    color = "black",
    shape = 21,
    alpha = .5
  ) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black")
  ) +
  ylab("Shannon") +
  xlab("") +
  scale_fill_brewer("Former", palette = "Dark2") +
  scale_color_brewer("Former", palette = "Dark2")
dev.off()

#this is the whole enchilada here. It might just be me but I find this really
#difficult to determine if this actually does what I want here.

lmefit1 <- glmmTMB(shannon ~ npolyp *
                     factor(tissue) + (1 | id) ,
                   na.omit(polyp2_obj$shannon.df))

summary(lmefit1) #sig

#mucosal former effects


lmefit2 <- glmmTMB(shannon ~
                     npolyp + (1 | id),
                   na.omit(subset(
                     polyp2_obj$shannon.df, subset = tissue == "Mucosal"
                   )))
summary(lmefit2) # no sig

#mucosal adenomas with location effects

lmefit3 <- glmmTMB(shannon ~
                     npolyp + location + (1 | id),
                   na.omit(subset(
                     polyp2_obj$shannon.df, subset = tissue == "Mucosal"
                   )))
summary(lmefit3) # no sig

# Note:  fecal and oral models only include main effects bc there aren't
# multiple samples per individual. We use stats::glm

#fecal npolyp
lmfit1 <- glm(shannon ~
                npolyp ,
              data = na.omit(subset(polyp2_obj$shannon.df, tissue == "Fecal")))

summary(lmfit1) #sig

#oral npolyp

lmfit2 <- glm(shannon ~
                npolyp ,
              data = na.omit(subset(polyp2_obj$shannon.df, tissue == "Oral")))
summary(lmfit2) #no sig


adenoma_shannon.point <- ggplot(na.omit(polyp2_obj$shannon.df),
                                aes(x = npolyp,
                                    y = shannon,
                                    color = tissue))

pdf("figs/shannon_regression.pdf",
    width = 14,
    height = 7)

adenoma_shannon.point +
  geom_point(size = 3, alpha = .7) +
  geom_smooth(method = "glm", aes(fill = tissue), alpha = .1) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1
  ) +
  ylab("Shannon") +
  xlab("Adenomas") +
  facet_wrap(. ~ tissue) +
  scale_color_brewer("Tissue", palette = "Dark2") +
  scale_fill_brewer("Tissue", palette = "Dark2")

dev.off()


# ANALYSIS: RICHNESS  -----------------------------------------------------


#Richness
polyp2_obj$richness <- as.data.frame(specnumber(polyp2_obj$data))
colnames(polyp2_obj$richness)[1] <- "richness"
polyp2_obj$richness$tissue <- polyp2_obj$meta$type
polyp2_obj$richness$location <- polyp2_obj$meta$location
polyp2_obj$richness$former <- factor(polyp2_obj$meta$polyp)
polyp2_obj$richness$npolyp <- polyp2_obj$meta$Adenoma
polyp2_obj$richness$polyptissue <- polyp2_obj$meta$polyp.tissue
polyp2_obj$richness$id <- polyp2_obj$meta$id

polyp2_richness.plot <- ggplot(na.omit(polyp2_obj$richness),
                               aes(x = tissue,
                                   y = richness,
                                   fill = former))

pdf("figs/richness_tissue_polyp.pdf")
polyp2_richness.plot +
  geom_boxplot() +
  geom_point(
    position = position_dodge(width = .75),
    color = "black",
    shape = 21,
    alpha = .5
  ) +
  theme(
    text = element_text(size = 24, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1.5,
    legend.position = "top"
  ) +
  ylab("Richness") +
  xlab("") +
  scale_fill_brewer("",
                    palette = "Dark2",
                    labels = c("Non-Former", "Former")) +
  scale_color_brewer("", palette = "Dark2")
dev.off()


polyp2_richness.boxplot <- ggplot(na.omit(polyp2_obj$richness),
                                  aes(x = location,
                                      y = richness,
                                      fill = former))

pdf("figs/richness_location_boxplot.pdf")
polyp2_richness.boxplot +
  geom_boxplot() +
  geom_point(
    position = position_dodge(width = .75),
    color = "black",
    shape = 21,
    alpha = .5
  ) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black")
  ) +
  ylab("Richness") +
  xlab("") +
  scale_fill_brewer("Former", palette = "Dark2") +
  scale_color_brewer("Former", palette = "Dark2")
dev.off()


adenoma_richness.point <- ggplot(na.omit(polyp2_obj$richness),
                                 aes(x = npolyp,
                                     y = richness,
                                     color = tissue))

# So I think that what we can say here is that the there are tissue level
# differences in richness (no surprise here). The relationship with number
# of polyps is interesting.

lmerich_fit1 <- glmmTMB(richness ~ npolyp * factor(tissue) +
                          (1 | id) ,
                        na.omit(polyp2_obj$richness))

summary(lmerich_fit1)#sig

#mucosal former effects
lmerich_fit2 <- glmmTMB(richness ~
                          npolyp + (1 | id),
                        na.omit(subset(polyp2_obj$richness, subset = tissue ==
                                         "Mucosal")))
summary(lmerich_fit2)#no sig

#mucosal adenomas with location effects

lmerichfit3 <- glmmTMB(richness ~
                         npolyp + location + (1 | id),
                       na.omit(subset(polyp2_obj$richness,  tissue ==
                                        "Mucosal")))
lmerichfit3.null <- glmmTMB(richness ~
                              (1 | id),
                            na.omit(subset(polyp2_obj$richness,  tissue ==
                                             "Mucosal")))

anova(lmerichfit3.null, lmerichfit3) #marginally sig
summary(lmerichfit3) # sig only for rectum

#fecal former
# Note: These models only include main effects because there aren't multiple
# samples for each individual
lmrich_fit1 <- glm(richness ~
                     npolyp ,
                   data = na.omit(subset(polyp2_obj$richness, tissue == "Fecal")))
summary(lmrich_fit1)#sig

#oral former
# Note: These models only include main effects because there aren't multiple
# samples for each individual

lmrich_fit2 <- glm(richness ~
                     npolyp ,
                   data = na.omit(subset(polyp2_obj$richness, tissue == "Oral")))
summary(lmrich_fit2) #no sig


# These results mirror those from the original (separate) analyses which indicate
# that only fecal alpha diversity associates with number of adenomas


pdf("figs/richness_regression.pdf",
    width = 14,
    height = 7)

adenoma_richness.point +
  geom_point(size = 3, alpha = .7) +
  geom_smooth(method = "glm", aes(fill = tissue), alpha = .1) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1
  ) +
  ylab("Richness") +
  xlab("Adenomas") +
  facet_wrap(. ~ tissue) +
  scale_color_brewer("Tissue", palette = "Dark2") +
  scale_fill_brewer("Tissue", palette = "Dark2")

dev.off()


# ANALYSIS: BETA DIVERSITY NMDS ---------------------------------------------

#Start with a simple all tissue ordination

#add metadata
polyp2_obj$mds.df$tissue       <- polyp2_obj$meta$type
polyp2_obj$mds.df$former       <- polyp2_obj$meta$polyp
polyp2_obj$mds.df$npolyp       <- polyp2_obj$meta$Adenoma
polyp2_obj$mds.df$polyp.tissue <- polyp2_obj$meta$polyp.tissue
polyp2_obj$mds.df$location     <- polyp2_obj$meta$location
polyp2_obj$mds.df$nbin         <-
  cut(
    polyp2_obj$mds.df$npolyp,
    # for plotting
    breaks = c(-1, 0.1, 2.1, 5, 10, 20),
    labels = c("0", "1-2", "3-5", "6-10", "10+")
  )

tissue_mds.ord <- ggplot(na.omit(polyp2_obj$mds.df),
                         aes(x = MDS1,
                             y = MDS2,
                             color = tissue))

pdf("figs/tissue_nmds.pdf")
tissue_mds.ord +
  geom_point(size = 3, alpha = .4) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(fill = NA),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1,
    legend.position = "top"
  ) +
  scale_color_brewer("", palette = "Dark2")
dev.off()



mucosal_mds.ord <-
  ggplot(na.omit(subset(polyp2_obj$mds.df, subset = tissue == "Mucosal")),
         aes(x = MDS1,
             y = MDS2,
             color = nbin))

pdf("figs/mucosal_nmds.pdf")
mucosal_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(fill = NA),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1
  ) +
  scale_alpha_manual(
    "",
    breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
    values =  c(0.1, 0.7, 0.7, 0.7, .7)
  ) +
  scale_color_brewer("Adenomas", palette = "Set1") +
  guides(alpha = FALSE)
dev.off()

#fecal

fecal_mds.ord <-
  ggplot(na.omit(subset(polyp2_obj$mds.df, subset = tissue == "Fecal")),
         aes(x = MDS1,
             y = MDS2,
             color = nbin))

pdf("figs/fecal_nmds.pdf")
fecal_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(fill = NA),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1
  ) +
  scale_alpha_manual(
    "",
    breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
    values =  c(0.1, 0.7, 0.7, 0.7, .7)
  ) +
  scale_color_brewer("Adenomas", palette = "Set1") +
  guides(alpha = FALSE)
dev.off()



#oral

oral_mds.ord <-
  ggplot(na.omit(subset(polyp2_obj$mds.df, subset = tissue == "Oral")),
         aes(x = MDS1,
             y = MDS2,
             color = nbin))

pdf("figs/oral_nmds.pdf")
oral_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_rect(fill = NA),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    aspect.ratio = 1
  ) +
  scale_alpha_manual(
    "",
    breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
    values =  c(0.1, 0.7, 0.7, 0.7, .7)
  ) +
  scale_color_brewer("Adenomas", palette = "Set1") +
  guides(alpha = FALSE)
dev.off()


# ANALYSIS: BETA DIVERSITY ADONIS (ASV) -------------------------------------

#start with one big model that includes all the data
set.seed(731)
tissue_npolyp.adonis <-
  adonis(
    polyp2_obj$data ~  polyp2_obj$meta$Adenoma * polyp2_obj$meta$type,
    permutations = 5000,
    method = "bray"
  )
tissue_npolyp.adonis

#this model indicates that there are significant impacts of tissue and number of
#adenomas on the microbiome as well as an interaction between these two effects


#lets look at mucosa separately
set.seed(731)
mucosal_npolyp.adonis <-
  adonis(
    polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),]), ] ~
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), "Adenoma"] +
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), "location"],
    permutations = 5000,
    method = "bray"
  )

mucosal_npolyp.adonis

#significant effect of adenoma # but not colon location.


#now fecal
set.seed(731)

fecal_npolyp.adonis <-
  adonis(polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),]), ] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

fecal_npolyp.adonis #not significant

#now oral
set.seed(731)

oral_npolyp.adonis <-
  adonis(polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"),]), ] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

oral_npolyp.adonis #not significant

# ANALYSIS: BETA DIVERSITY ADONIS (Genus) -------------------------------------

#start with one big model that includes all the data
set.seed(731)
gtissue_npolyp.adonis <-
  adonis(
    polyp2_obj$phylotype$genus ~  polyp2_obj$meta$Adenoma * polyp2_obj$meta$type,
    permutations = 5000,
    method = "bray"
  )
gtissue_npolyp.adonis

#this model indicates that there are significant impacts of tissue and number of
#adenomas on the microbiome as well as an interaction between these two effects


#lets look at mucosa separately
set.seed(731)
gmucosal_npolyp.adonis <-
  adonis(
    polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),]), ] ~
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), "Adenoma"] +
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), "location"],
    permutations = 5000,
    method = "bray"
  )

gmucosal_npolyp.adonis

#significant effect of adenoma # but not colon location.


#now fecal
set.seed(731)

gfecal_npolyp.adonis <-
  adonis(
    polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),]), ] ~
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), "Adenoma"] ,
    permutations = 5000,
    method = "bray"
  )

gfecal_npolyp.adonis #not significant

#now oral
set.seed(731)

goral_npolyp.adonis <-
  adonis(
    polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"),]), ] ~
      polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), "Adenoma"] ,
    permutations = 5000,
    method = "bray"
  )

goral_npolyp.adonis #not significant


#so only the mucosal adonis and the big model are significant

# ANALYSIS: TAXA MODELS SETUP ------------------------------------------

# For this analysis the data needs to be split and filtered to ensure quality
# results.

# Sorry for the mess here... Start by subsetting based on type of tissue
# and filtering to get taxa that are present in at least 20% of samples.
# This threshold  was chosen because ~22% of mucosal tissue samples are from
# adenomas.By selecting a threshold that is slightly lower than this we ensure
# that taxa that are only present on adenoma tissue will not be filtered. For
# consistency we will apply this same threshold to other tissues.


#Fecal ASV
fecal_adenoma_asv.df <-
  polyp2_obj$data[which(rownames(polyp2_obj$data) %in%
                          rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),])), ]

#Filter
fecal_adenoma_asv.df <-
  fecal_adenoma_asv.df[, apply(fecal_adenoma_asv.df, 2,
                               function(x) {
                                 sum(as.logical(x)) >
                                   round(nrow(fecal_adenoma_asv.df) * .2)
                               })]

#Fecal Genus
fecal_adenoma_genus.df <-
  polyp2_obj$phylotype$genus[which(rownames(polyp2_obj$phylotype$genus) %in%
                                     rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),])), ]

#Filter
fecal_adenoma_genus.df <-
  fecal_adenoma_genus.df[, apply(fecal_adenoma_genus.df,
                                 2,
                                 function(x) {
                                   sum(as.logical(x)) >
                                     round(nrow(fecal_adenoma_genus.df) * .2)
                                 })]



#Mucosal ASV
mucosal_adenoma_asv.df <-
  polyp2_obj$data[which(rownames(polyp2_obj$data) %in%
                          rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),])), ]

#Filter
mucosal_adenoma_asv.df <-
  mucosal_adenoma_asv.df[, apply(mucosal_adenoma_asv.df, 2,
                                 function(x) {
                                   sum(as.logical(x)) >
                                     round(nrow(mucosal_adenoma_asv.df) * .2)
                                 })]

#MucosalGenus
mucosal_adenoma_genus.df <-
  polyp2_obj$phylotype$genus[which(rownames(polyp2_obj$phylotype$genus) %in%
                                     rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),])), ]

#Filter
mucosal_adenoma_genus.df <-
  mucosal_adenoma_genus.df[, apply(mucosal_adenoma_genus.df, 2,
                                   function(x) {
                                     sum(as.logical(x)) >
                                       round(nrow(mucosal_adenoma_genus.df) * .2)
                                   })]


#Oral ASV
oral_adenoma_asv.df <-
  polyp2_obj$data[which(rownames(polyp2_obj$data) %in%
                          rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), ])),]

#Filter
oral_adenoma_asv.df <-
  oral_adenoma_asv.df[, apply(oral_adenoma_asv.df, 2,
                              function(x) {
                                sum(as.logical(x)) >
                                  round(nrow(oral_adenoma_asv.df) * .2)
                              })]

#Oral Genus
oral_adenoma_genus.df <-
  polyp2_obj$phylotype$genus[which(rownames(polyp2_obj$phylotype$genus) %in%
                                     rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"),])), ]

#Filter

oral_adenoma_genus.df <-
  oral_adenoma_genus.df[, apply(oral_adenoma_genus.df, 2,
                                function(x) {
                                  sum(as.logical(x)) >
                                    round(nrow(oral_adenoma_genus.df) * .2)
                                })]

#make master objects
asv.models.obj   <- NULL
genus.models.obj <- NULL

# ANALYSIS: TAXA MODELS MUCOSAL GENUS -------------------------------------

#nadenoma models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_genus.df[, i] ~ Adenoma + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )

    mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <-
      fit1
  }, error = function(e)
    e)
}

#assign to proper slot so we can recycle the object
genus.models.obj$nadenoma <- mucosal_adenoma_models

#former models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_genus.df[, i] ~ factor(polyp) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )


    mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <-
      fit1
  }, error = function(e)
    e)
}

#assign to proper slot so we can recycle the object
genus.models.obj$former <- mucosal_adenoma_models

#polyp.tissue models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_genus.df[, i] ~ factor(polyp.tissue) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )


    mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <-
      fit1
  }, error = function(e)
    e)
}

#assign to proper slot so we can recycle the object
genus.models.obj$polyp.tissue <- mucosal_adenoma_models

#location models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_genus.df[, i] ~ factor(location) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )


    mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <-
      fit1
  }, error = function(e)
    e)
}

#assign to proper slot so we can recycle the object
genus.models.obj$location <- mucosal_adenoma_models

#basic stats on fdr control genus
sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$nadenoma,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$former,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$polyp.tissue,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$location,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

# ANALYSIS: TAXA MODELS MUCOSAL ASV ---------------------------------------

#nadenoma models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_asv.df[, i] ~ Adenoma + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )

    mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$nadenoma <- mucosal_adenoma_models

#former models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_asv.df[, i] ~ factor(polyp) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )

    mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$former <- mucosal_adenoma_models

#polyp.tissue models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_asv.df[, i] ~ factor(polyp.tissue) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )

    mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$polyp.tissue <- mucosal_adenoma_models

#location models
set.seed(731)
mucosal_adenoma_models <- NULL

for (i in 1:ncol(mucosal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glmmTMB(
        mucosal_adenoma_asv.df[, i] ~ factor(location) + (1 | id),
        data = subset(polyp2_obj$meta, type == "Mucosal"),
        family = gaussian
      )

    mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$location <- mucosal_adenoma_models

#basic stats on fdr control ASV
sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$nadenoma,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$former,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$polyp.tissue,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$location,
    FUN = function(x) {
      coefficients(summary(x))$cond[2, 4]
    }
  )
)$qvalue < .1))

# ANALYSIS: TAXA MODELS FECAL GENUS -------------------------------------

# Fecal and oral models will be a little different. This is because there is
# only 1 samples / ID and thus using ID as a random effect is probably not a
# good idea. I will just use a standard neg binomial glm for these models.

#fecal #adenoma

set.seed(731)
fecal_adenoma_models <- NULL

for (i in 1:ncol(fecal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(fecal_adenoma_genus.df[, i] ~ Adenoma ,
             data = subset(polyp2_obj$meta, type == "Fecal"))

    fecal_adenoma_models[[colnames(fecal_adenoma_genus.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

genus.models.obj$fecal_nadenoma <- fecal_adenoma_models

#fecal former
set.seed(731)
fecal_adenoma_models <- NULL

for (i in 1:ncol(fecal_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(fecal_adenoma_genus.df[, i] ~ factor(polyp) ,
             data = subset(polyp2_obj$meta, type == "Fecal"))

    fecal_adenoma_models[[colnames(fecal_adenoma_genus.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

genus.models.obj$fecal_former <- fecal_adenoma_models

# ANALYSIS: TAXA MODELS FECAL ASV ---------------------------------------

#fecal asv #adenoma

set.seed(731)
fecal_adenoma_models <- NULL

for (i in 1:ncol(fecal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(fecal_adenoma_asv.df[, i] ~ Adenoma ,
             data = subset(polyp2_obj$meta, type == "Fecal"))

    fecal_adenoma_models[[colnames(fecal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$fecal_nadenoma <- fecal_adenoma_models


#fecal former
set.seed(731)
fecal_adenoma_models <- NULL

for (i in 1:ncol(fecal_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(fecal_adenoma_asv.df[, i] ~ factor(polyp) ,
             data = subset(polyp2_obj$meta, type == "Fecal"))

    fecal_adenoma_models[[colnames(fecal_adenoma_asv.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

asv.models.obj$fecal_former <- fecal_adenoma_models
fecal_adenoma_models <- NULL

# ANALYSIS: TAXA MODELS ORAL GENUS --------------------------------------

#oral #adenoma

set.seed(731)
oral_adenoma_models <- NULL

for (i in 1:ncol(oral_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(oral_adenoma_genus.df[, i] ~ Adenoma ,
             data = subset(polyp2_obj$meta, type == "Oral"))

    oral_adenoma_models[[colnames(oral_adenoma_genus.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

genus.models.obj$oral_nadenoma <- oral_adenoma_models

#oral former
set.seed(731)
oral_adenoma_models <- NULL

for (i in 1:ncol(oral_adenoma_genus.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(oral_adenoma_genus.df[, i] ~ factor(polyp) ,
             data = subset(polyp2_obj$meta, type == "Oral"))

    oral_adenoma_models[[colnames(oral_adenoma_genus.df)[i]]] <-
      fit1
  },
  error = function(e)
    e)
}

genus.models.obj$oral_former <- oral_adenoma_models

#basic stats on fdr control genus
sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$fecal_nadenoma,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$fecal_former,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$oral_nadenoma,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    genus.models.obj$oral_former,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))
# ANALYSIS: TAXA MODELS ORAL ASV ----------------------------------------

#oral #adenoma

set.seed(731)
oral_adenoma_models <- NULL

for (i in 1:ncol(oral_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(oral_adenoma_asv.df[, i] ~ Adenoma ,
             data = subset(polyp2_obj$meta, type == "Oral"))

    oral_adenoma_models[[colnames(oral_adenoma_asv.df)[i]]] <- fit1
  },
  error = function(e)
    e)
}

asv.models.obj$oral_nadenoma <- oral_adenoma_models


#oral former
set.seed(731)
oral_adenoma_models <- NULL

for (i in 1:ncol(oral_adenoma_asv.df)) {
  fit1 <- NULL
  tryCatch({
    fit1 <-
      glm(oral_adenoma_asv.df[, i] ~ factor(polyp) ,
             data = subset(polyp2_obj$meta, type == "Oral"))

    oral_adenoma_models[[colnames(oral_adenoma_asv.df)[i]]] <- fit1
  },
  error = function(e)
    e)
}

asv.models.obj$oral_former <- oral_adenoma_models
oral_adenoma_models <- NULL

#basic stats on fdr control
sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$fecal_nadenoma,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$fecal_former,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$oral_nadenoma,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

sum(na.omit(qvalue::qvalue(
  sapply(
    asv.models.obj$oral_former,
    FUN = function(x) {
      coefficients(summary(x))[2, 4]
    }
  )
)$qvalue < .1))

# ANALYSIS: GENUS MODEL AGGREGATION -------------------------------------

model_summary_genus.df <-
  as.data.frame(matrix(
    nrow =
      sum(sapply(genus.models.obj, FUN = function(x){length(x)})),
    ncol = 7 ))

colnames(model_summary_genus.df) <- c("taxa",
                                      "tissue",
                                      "test",
                                      "est",
                                      "stat_t_or_z",
                                      "pval",
                                      "qval")

counter <- 1
for(i in names(genus.models.obj)[1:3]){
  t1 <- counter
  for( j in 1:length(genus.models.obj[[i]])){
    fit <- summary(genus.models.obj[[i]][[j]])

    model_summary_genus.df[counter, 1 ] <- names(genus.models.obj[[i]])[j]
    model_summary_genus.df[counter, 2 ] <- "Mucosal"
    model_summary_genus.df[counter, 3 ] <- i
    model_summary_genus.df[counter, 4 ] <- fit$coefficients$cond[2, 1]
    model_summary_genus.df[counter, 5 ] <- fit$coefficients$cond[2, 3]
    model_summary_genus.df[counter, 6 ] <- fit$coefficients$cond[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_genus.df[t1:t2,7] <-
    p.adjust(model_summary_genus.df[t1:t2,6], method = "fdr")
    #qvalue(model_summary_genus.df[t1:t2,6])$qvalue

}

#mucosal location
#this is a little messy because we are really interested in the overall
#impact of location on the taxa abundance. to determine this effect we
#will use an anova and not report the estimates or tvals bc these will
#no be generated with chisq

t1 <- counter
for( j in 1:length(genus.models.obj[[4]])){
    fit <- car::Anova(genus.models.obj[[4]][[j]])

    model_summary_genus.df[counter, 1 ] <- names(genus.models.obj[[4]])[j]
    model_summary_genus.df[counter, 2 ] <- "Mucosal"
    model_summary_genus.df[counter, 3 ] <- "Location"
    model_summary_genus.df[counter, 4 ] <- NA
    model_summary_genus.df[counter, 5 ] <- NA
    model_summary_genus.df[counter, 6 ] <- fit$`Pr(>Chisq)`
    #move to populate the next row
    counter <- counter + 1
}

t2 <- counter - 1
#add qvalue
model_summary_genus.df[t1:t2,7] <-
  p.adjust(model_summary_genus.df[t1:t2,6], method = "fdr")
  #qvalue(model_summary_genus.df[t1:t2,6])$qvalue

t1 <- NULL
t2 <- NULL

#now on to fecal and oral which are different because they use glm.nb
#fecal
for(i in names(genus.models.obj)[5:6]){
  t1 <- counter
  for( j in 1:length(genus.models.obj[[i]])){
    fit <- summary(genus.models.obj[[i]][[j]])

    model_summary_genus.df[counter, 1 ] <- names(genus.models.obj[[i]])[j]
    model_summary_genus.df[counter, 2 ] <- "Fecal"
    model_summary_genus.df[counter, 3 ] <- i
    model_summary_genus.df[counter, 4 ] <- coefficients(fit)[2, 1]
    model_summary_genus.df[counter, 5 ] <- coefficients(fit)[2, 3]
    model_summary_genus.df[counter, 6 ] <- coefficients(fit)[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_genus.df[t1:t2,7] <-
    p.adjust(model_summary_genus.df[t1:t2,6], method = "fdr")
    #qvalue(model_summary_genus.df[t1:t2,6])$qvalue

}


for(i in names(genus.models.obj)[7:8]){
  t1 <- counter
  for( j in 1:length(genus.models.obj[[i]])){
    fit <- summary(genus.models.obj[[i]][[j]])

    model_summary_genus.df[counter, 1 ] <- names(genus.models.obj[[i]])[j]
    model_summary_genus.df[counter, 2 ] <- "Oral"
    model_summary_genus.df[counter, 3 ] <- i
    model_summary_genus.df[counter, 4 ] <- coefficients(fit)[2, 1]
    model_summary_genus.df[counter, 5 ] <- coefficients(fit)[2, 3]
    model_summary_genus.df[counter, 6 ] <- coefficients(fit)[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_genus.df[t1:t2,7] <-
    p.adjust(model_summary_genus.df[t1:t2,6], method = "fdr")
   # qvalue(model_summary_genus.df[t1:t2,6])$qvalue
}

#reorder
model_summary_genus.df <-
  model_summary_genus.df[order(model_summary_genus.df$test,
                               model_summary_genus.df$qval),]


# ANALYSIS: TAXA VISUALIZATIONS GENUS -------------------------------------

nadenoma_test_all.df <- subset(model_summary_genus.df,
                               test %in% c("nadenoma",
                                           "fecal_nadenoma",
                                           "oral_nadenoma"
                                           )
                               )

keeps.taxa <-
  unlist(unique(subset(nadenoma_test_all.df, qval < 0.1, select = taxa)))

nadenoma_test_all.df <-
  nadenoma_test_all.df[which(nadenoma_test_all.df$taxa %in% keeps.taxa),]


nadenoma_test_all.plot <- ggplot(nadenoma_test_all.df,
                                 aes(x = tissue,
                                     y = taxa,
                                     fill = est))
nadenoma_test_all.plot +
  geom_tile()+
  #coord_flip()+
  scale_fill_distiller(palette = "PRGn", limits= c(-20,20), oob=squish )+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "right"
  )+
  xlab("")+
  ylab("")

# ANALYSIS: ASV MODEL AGGREGATION -------------------------------------

model_summary_asv.df <-
  as.data.frame(matrix(
    nrow =
      sum(sapply(asv.models.obj, FUN = function(x){length(x)})),
    ncol = 7 ))

colnames(model_summary_asv.df) <- c("taxa",
                                    "tissue",
                                    "test",
                                    "est",
                                    "stat_t_or_z",
                                    "pval",
                                    "qval")

counter <- 1
for(i in names(asv.models.obj)[1:3]){
  t1 <- counter
  for( j in 1:length(asv.models.obj[[i]])){
    fit <- summary(asv.models.obj[[i]][[j]])

    model_summary_asv.df[counter, 1 ] <- names(asv.models.obj[[i]])[j]
    model_summary_asv.df[counter, 2 ] <- "Mucosal"
    model_summary_asv.df[counter, 3 ] <- i
    model_summary_asv.df[counter, 4 ] <- fit$coefficients$cond[2, 1]
    model_summary_asv.df[counter, 5 ] <- fit$coefficients$cond[2, 3]
    model_summary_asv.df[counter, 6 ] <- fit$coefficients$cond[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_asv.df[t1:t2,7] <-
    p.adjust(model_summary_asv.df[t1:t2,6], method = "fdr")
  #qvalue(model_summary_asv.df[t1:t2,6])$qvalue

}

#mucosal location
#this is a little messy because we are really interested in the overall
#impact of location on the taxa abundance. to determine this effect we
#will use an anova and not report the estimates or tvals bc these will
#no be generated with chisq

t1 <- counter
for( j in 1:length(asv.models.obj[[4]])){
  fit <- car::Anova(asv.models.obj[[4]][[j]])

  model_summary_asv.df[counter, 1 ] <- names(asv.models.obj[[4]])[j]
  model_summary_asv.df[counter, 2 ] <- "Mucosal"
  model_summary_asv.df[counter, 3 ] <- "Location"
  model_summary_asv.df[counter, 4 ] <- NA
  model_summary_asv.df[counter, 5 ] <- NA
  model_summary_asv.df[counter, 6 ] <- fit$`Pr(>Chisq)`
  #move to populate the next row
  counter <- counter + 1
}

t2 <- counter - 1
#add qvalue
model_summary_asv.df[t1:t2,7] <-
  p.adjust(model_summary_asv.df[t1:t2,6], method = "fdr")
#qvalue(model_summary_asv.df[t1:t2,6])$qvalue

t1 <- NULL
t2 <- NULL

#now on to fecal and oral which are different because they use glm.nb
#fecal
for(i in names(asv.models.obj)[5:6]){
  t1 <- counter
  for( j in 1:length(asv.models.obj[[i]])){
    fit <- summary(asv.models.obj[[i]][[j]])

    model_summary_asv.df[counter, 1 ] <- names(asv.models.obj[[i]])[j]
    model_summary_asv.df[counter, 2 ] <- "Fecal"
    model_summary_asv.df[counter, 3 ] <- i
    model_summary_asv.df[counter, 4 ] <- coefficients(fit)[2, 1]
    model_summary_asv.df[counter, 5 ] <- coefficients(fit)[2, 3]
    model_summary_asv.df[counter, 6 ] <- coefficients(fit)[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_asv.df[t1:t2,7] <-
    p.adjust(model_summary_asv.df[t1:t2,6], method = "fdr")
  #qvalue(model_summary_asv.df[t1:t2,6])$qvalue

}


for(i in names(asv.models.obj)[7:8]){
  t1 <- counter
  for( j in 1:length(asv.models.obj[[i]])){
    fit <- summary(asv.models.obj[[i]][[j]])

    model_summary_asv.df[counter, 1 ] <- names(asv.models.obj[[i]])[j]
    model_summary_asv.df[counter, 2 ] <- "Oral"
    model_summary_asv.df[counter, 3 ] <- i
    model_summary_asv.df[counter, 4 ] <- coefficients(fit)[2, 1]
    model_summary_asv.df[counter, 5 ] <- coefficients(fit)[2, 3]
    model_summary_asv.df[counter, 6 ] <- coefficients(fit)[2, 4]
    #move to populate the next row
    counter <- counter + 1
  }
  t2 <- counter - 1
  #add qvalue
  model_summary_asv.df[t1:t2,7] <-
    p.adjust(model_summary_asv.df[t1:t2,6], method = "fdr")
  # qvalue(model_summary_asv.df[t1:t2,6])$qvalue
}

#reorder
model_summary_asv.df <-
  model_summary_asv.df[order(model_summary_asv.df$test,
                             model_summary_asv.df$qval),]


# ANALYSIS: TAXA VISUALIZATIONS ASV -------------------------------------

asv_nadenoma_test_all.df <- subset(model_summary_asv.df,
                                   test %in% c("nadenoma",
                                               "fecal_nadenoma",
                                               "oral_nadenoma"
                                   )
)

keeps.taxa <-
  unlist(unique(subset(asv_nadenoma_test_all.df, qval < 0.1, select = taxa)))

asv_nadenoma_test_all.df <-
  asv_nadenoma_test_all.df[which(asv_nadenoma_test_all.df $taxa %in% keeps.taxa),]

t.vec <- NULL
for(i in 1:length(asv_nadenoma_test_all.df$taxa)){
  t <- NULL
  xdf <- polyp2_tax.df[(asv_nadenoma_test_all.df$taxa)[i],6:8]
  if(is.na(xdf["Species"])){
    if(is.na(xdf["Genus"])){
      t <-  paste0(xdf["Family"], "(F)")
    }else{
      t <- paste0(xdf["Genus"], "(G)")
    }
  }else{
    t <- paste0(xdf["Genus"]," ", xdf["Species"], "(S)")
  }
  t <- paste0(t, "-", (asv_nadenoma_test_all.df$taxa)[i])
  t.vec <- c(t.vec, t)
}

asv_nadenoma_test_all.df$tax.name <- t.vec
asv_nadenoma_test_all.df$sig <-
  sapply(asv_nadenoma_test_all.df$qval,
         FUN = function(x){if(x < 0.1 ){print("*")}else{print("")}}
         )

pdf("figs/asv_models_sig.pdf")
asv_nadenoma_test_all.plot <- ggplot(asv_nadenoma_test_all.df,
                                     aes(x = tissue,
                                         y = tax.name,
                                         fill = est))
asv_nadenoma_test_all.plot +
  geom_tile()+
  geom_text(aes(label = sig),nudge_y = -.1,
            size =6,
            color = "grey65")+
  #coord_flip()+
  scale_fill_distiller(palette = "PRGn",
                       limits= c(-10,20),
                       oob=squish )+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border     = element_rect(fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
  )+
  xlab("")+
  ylab("")+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))
dev.off()

# ANALYSIS: RANDOM FORESTS ------------------------------------------------

#Mucosal

set.seed(731)
mucosal_former_asv.rf <- randomForest(x = mucosal_adenoma_asv.df,
                                      y = factor(unlist(
                                        subset(polyp2_obj$meta,
                                               type == "Mucosal",
                                               select = polyp)
                                      )),
                                      na.action = na.omit,
                                      importance = T,
                                      ntree =5000
                                      )

mucosal_former_genus.rf <-
  randomForest(x = mucosal_adenoma_genus.df,
               y = factor(unlist(
                 subset(polyp2_obj$meta,
                        type == "Mucosal",
                        select = polyp)
               )),
               na.action = na.omit,
               importance = T,
               ntree =5000
  )

#Fecal
fecal_former_asv.rf <- randomForest(x = fecal_adenoma_asv.df,
                                    y = factor(unlist(
                                      subset(polyp2_obj$meta,
                                             type == "Fecal",
                                             select = polyp)
                                    )),
                                    na.action = na.omit,
                                    importance = T,
                                    ntree =5000
)

fecal_former_genus.rf <- randomForest(x = fecal_adenoma_genus.df,
                                      y = factor(unlist(
                                        subset(polyp2_obj$meta,
                                               type == "Fecal",
                                               select = polyp)
                                      )),
                                      na.action = na.omit,
                                      importance = T,
                                      ntree =5000
)

#Oral
oral_former_asv.rf <- randomForest(x = oral_adenoma_asv.df,
                                   y = factor(unlist(
                                     subset(polyp2_obj$meta,
                                            type == "Oral",
                                            select = polyp)
                                   )),
                                   na.action = na.omit,
                                   importance = T,
                                   ntree =5000
)

oral_former_genus.rf <- randomForest(x = oral_adenoma_genus.df,
                                     y = factor(unlist(
                                       subset(polyp2_obj$meta,
                                              type == "Oral", select = polyp)
                                     )),
                                     na.action = na.omit,
                                     importance = T,
                                     ntree =5000
)


# ANALYSIS: RANDOM FORESTS VISUALS ASV ROC CURVES ----------------------------

#mucosal
mucosal_asv.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Mucosal", select = polyp)
)), mucosal_former_asv.rf$votes[,2])

#fecal
fecal_asv.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Fecal", select = polyp)
)), fecal_former_asv.rf$votes[,2])

#oral
oral_asv.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Oral", select = polyp)
)), oral_former_asv.rf$votes[, 2])

#get auc

auc(mucosal_asv.roc)
auc(fecal_asv.roc)
auc(oral_asv.roc)

#make data frame for plotting

rf_roc.df <- NULL
rf_roc.df$variable    <- c(
  rep("Mucosal",
      times = length(mucosal_asv.roc$sensitivities)),
  rep("Fecal",
      times = length(fecal_asv.roc$sensitivities)),
  rep("Oral",
      times = length(oral_asv.roc$sensitivities))
)

#add sensitivity
rf_roc.df$sensitivity <-
  c(
    mucosal_asv.roc$sensitivities,
    fecal_asv.roc$sensitivities,
    oral_asv.roc$sensitivities
  )

#add specificity
rf_roc.df$specificity <-
  c(
    mucosal_asv.roc$specificities,
    fecal_asv.roc$specificities,
    oral_asv.roc$specificities
  )
rf_roc.df <- as.data.frame(rf_roc.df)


rf_roc.plot <- ggplot(rf_roc.df,
                      aes(x = specificity,
                          y = sensitivity,
                          color = variable))

pdf("figs/asv_rf_roc.pdf")

rf_roc.plot +
  geom_abline(slope =1, intercept = 1, size = 1, alpha = .4)+
  geom_path(size = 1.5)+
  scale_x_reverse(expand = c(0.01,0.02))+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
  ) +
  scale_y_continuous(expand = c(0.0, 0.01))+
  scale_color_brewer(palette = "Dark2")+
  ylab("Sensitivity")+
  xlab("Specificity")
dev.off()

# ANALYSIS: RANDOM FORESTS VISUALS GENUS ROC CURVES ----------------------------

#mucosal
mucosal_genus.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Mucosal", select = polyp)
)), mucosal_former_genus.rf$votes[,2])

#fecal
fecal_genus.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Fecal", select = polyp)
)), fecal_former_genus.rf$votes[,2])

#oral
oral_genus.roc <- roc(factor(unlist(
  subset(polyp2_obj$meta,
         type == "Oral", select = polyp)
)), oral_former_genus.rf$votes[, 2])

#get auc

auc(mucosal_genus.roc)
auc(fecal_genus.roc)
auc(oral_genus.roc)

#make data frame for plotting

genus_rf_roc.df <- NULL
genus_rf_roc.df$variable    <- c(
  rep("Mucosal",
      times = length(mucosal_genus.roc$sensitivities)),
  rep("Fecal",
      times = length(fecal_genus.roc$sensitivities)),
  rep("Oral",
      times = length(oral_genus.roc$sensitivities))
)

#add sensitivity
genus_rf_roc.df$sensitivity <-
  c(
    mucosal_genus.roc$sensitivities,
    fecal_genus.roc$sensitivities,
    oral_genus.roc$sensitivities
  )

#add specificity
genus_rf_roc.df$specificity <-
  c(
    mucosal_genus.roc$specificities,
    fecal_genus.roc$specificities,
    oral_genus.roc$specificities
  )
genus_rf_roc.df <- as.data.frame(rf_roc.df)


genus_rf_roc.plot <- ggplot(genus_rf_roc.df,
                      aes(x = specificity,
                          y = sensitivity,
                          color = variable))

pdf("figs/genus_rf_roc.pdf")
genus_rf_roc.plot +
  geom_abline(slope =1, intercept = 1, size = 1, alpha = .4)+
  geom_path(size = 1.5)+
  scale_x_reverse(expand = c(0.01,0.02))+
  theme(text = element_text(size=20, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "top"
  ) +
  scale_y_continuous(expand = c(0.0, 0.01))+
  scale_color_brewer(palette = "Dark2")+
  ylab("Sensitivity")+
  xlab("Specificity")
dev.off()

# ANALYSIS: RANDOM FORESTS VISUALS MUCOSAL ASV IMPORTANCE PLOT ----------------

#Now lets look at what might be important to these models
mucosal_asv_imp.df <- mucosal_former_asv.rf$importance
mucosal_asv_imp.df <- as.data.frame(mucosal_asv_imp.df)

mucosal_asv_imp.df <-
  mucosal_asv_imp.df[order(mucosal_asv_imp.df$MeanDecreaseAccuracy,
                           decreasing = T), ]

mucosal_asv_imp.df <- mucosal_asv_imp.df[1:20, 3, drop = F]

mucosal_asv_imp.df$names <- rownames(mucosal_asv_imp.df)
mucosal_asv_imp.df$sd <-
  mucosal_former_asv.rf$importanceSD[mucosal_asv_imp.df$names, 3]



mucosal_asv_imp.df$mda_scale <-
  mucosal_asv_imp.df$MeanDecreaseAccuracy /
  mucosal_asv_imp.df$sd

mucosal_asv_imp.df <-
  mucosal_asv_imp.df[order(mucosal_asv_imp.df$mda_scale, decreasing = T), ]

mucosal_asv_imp.df$names <- factor(mucosal_asv_imp.df$names,
                                   levels = rev(mucosal_asv_imp.df$names))

pdf("figs/mucosal_asv_varimp.pdf")

mucosal_asv_imp.plot <- ggplot(mucosal_asv_imp.df, aes(x = names,
                                                       y = mda_scale))

mucosal_asv_imp.plot +
  geom_point(size = 4,
             alpha = .9,
             color = "#D95F02") +
  coord_flip() +
  xlab("") +
  ylab("Mean Decrease in Accuracy") +
  ylim(c(0, 45)) +
  theme(
    text = element_text(size = 18),
    panel.border     = element_rect(fill = NA),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    aspect.ratio = 1
  )
dev.off()


# ANALYSIS: RANDOM FORESTS VISUALS MUCOSAL GENUS IMPORTANCE PLOT --------------

#Now lets look at what might be important to these models
mucosal_genus_imp.df <- mucosal_former_genus.rf$importance
mucosal_genus_imp.df <- as.data.frame(mucosal_genus_imp.df)

mucosal_genus_imp.df <-
  mucosal_genus_imp.df[order(mucosal_genus_imp.df$MeanDecreaseAccuracy,
                             decreasing = T), ]

mucosal_genus_imp.df <- mucosal_genus_imp.df[1:20, 3, drop = F]

mucosal_genus_imp.df$names <- rownames(mucosal_genus_imp.df)
mucosal_genus_imp.df$sd <-
  mucosal_former_genus.rf$importanceSD[mucosal_genus_imp.df$names, 3]



mucosal_genus_imp.df$mda_scale <-
  mucosal_genus_imp.df$MeanDecreaseAccuracy /
  mucosal_genus_imp.df$sd

mucosal_genus_imp.df <-
  mucosal_genus_imp.df[order(mucosal_genus_imp.df$mda_scale, decreasing = T), ]

mucosal_genus_imp.df$names <- factor(mucosal_genus_imp.df$names,
                                     levels = rev(mucosal_genus_imp.df$names))

pdf("figs/mucosal_genus_varimp.pdf")

mucosal_genus_imp.plot <-
  ggplot(mucosal_genus_imp.df, aes(x = names,
                                   y = mda_scale))

mucosal_genus_imp.plot +
  geom_point(size = 4,
             alpha = .9,
             color = "#D95F02") +
  coord_flip() +
  xlab("") +
  ylab("Mean Decrease in Accuracy") +
  ylim(c(0, 60)) +
  theme(
    text = element_text(size = 18),
    panel.border     = element_rect(fill = NA),
    panel.grid.major = element_line(color = "grey97"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    aspect.ratio = 1
  )
dev.off()


# ANALYSIS: RANDOM FOREST WITH LEAVE ONE OUT ASV  -----------------------------

# Although over fitting is less of a problem with random forest I will still
# attempt some cross validation here. I will try using a 90% test 10% validate
# split for all samples
set.seed(731)

error.vec <- NULL
for (i in 1:1000) {
  model <-
    sample(row.names(mucosal_adenoma_asv.df),
           size = round(nrow(mucosal_adenoma_asv.df) * .9))
  cs <-
    row.names(mucosal_adenoma_asv.df)[-which(row.names(mucosal_adenoma_asv.df) %in% model)]


  #make the train set
  mucosal_asv_model.rf <-
    mucosal_adenoma_asv.df[which(rownames(mucosal_adenoma_asv.df) %in% model), ]

  mucosal_asv_meta_model.rf <-
    polyp2_obj$meta[which(rownames(polyp2_obj$meta) %in% model), ]

  #make the test set
  mucosal_asv_cv.rf <-
    mucosal_adenoma_asv.df[which(rownames(mucosal_adenoma_asv.df) %in% cs), ]

  mucosal_asv_meta_cv.rf <-
    polyp2_obj$meta[which(rownames(polyp2_obj$meta) %in% cs), ]

  x <- randomForest(x = mucosal_asv_model.rf,
                    y = factor(mucosal_asv_meta_model.rf$polyp))

  y <- predict(x, newdata = mucosal_asv_cv.rf)

  z <-
    1 - (sum(y == factor(mucosal_asv_meta_cv.rf$polyp)) / length(y))
  error.vec <- c(error.vec, z)
}

mean(error.vec)

ggplot(data = as.data.frame(error.vec), aes(x = error.vec))+
  geom_histogram(binwidth = .05, color = "white")

# ANALYSIS: RANDOM FOREST WITH LEAVE ONE OUT GENUS ---------------------------

# Although over fitting is less of a problem with random forest I will still
# attempt some cross validation here. I will try using a 90% test 10% validate
# split for all samples
set.seed(731)

gerror.vec <- NULL
for (i in 1:1000) {
  model <-
    sample(row.names(mucosal_adenoma_genus.df),
           size = round(nrow(mucosal_adenoma_genus.df) * .9))
  cs <-
    row.names(mucosal_adenoma_genus.df)[-which(row.names(mucosal_adenoma_genus.df) %in% model)]


  #make the train set
  mucosal_genus_model.rf <-
    mucosal_adenoma_genus.df[which(rownames(mucosal_adenoma_genus.df) %in% model), ]

  mucosal_genus_meta_model.rf <-
    polyp2_obj$meta[which(rownames(polyp2_obj$meta) %in% model), ]

  #make the test set
  mucosal_genus_cv.rf <-
    mucosal_adenoma_genus.df[which(rownames(mucosal_adenoma_genus.df) %in% cs), ]

  mucosal_genus_meta_cv.rf <-
    polyp2_obj$meta[which(rownames(polyp2_obj$meta) %in% cs), ]

  x <- randomForest(x = mucosal_genus_model.rf,
                    y = factor(mucosal_genus_meta_model.rf$polyp))

  y <- predict(x, newdata = mucosal_genus_cv.rf)

  z <-
    1 - (sum(y == factor(mucosal_genus_meta_cv.rf$polyp)) / length(y))
  gerror.vec <- c(gerror.vec, z)
}

mean(gerror.vec)

ggplot(data = as.data.frame(gerror.vec), aes(x = gerror.vec))+
  geom_histogram(binwidth = .05, color = "white")



# SANDBOX -----------------------------------------------------------------
