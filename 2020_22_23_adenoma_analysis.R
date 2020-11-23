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
options("stringsAsFactors"=F)

# FUNCTIONS ---------------------------------------------------------------

###
#        Function diversity_analysis             #
###

#this function performs alpha and beta diversity analysis on the
#user provided data frame

#depends on vegan

diversity_analysis <- function(obj){
  obj$shannon <- vegan::diversity(obj$data, index = "shannon", MARGIN = 1)
  obj$simpson <- vegan::diversity(obj$data, index = "simpson", MARGIN = 1)
  obj$invsimpson <- vegan::diversity(obj$data, index = "invsimpson", MARGIN = 1)
  return(obj)
}

###
#              Function ordinate                 #
###

#This function creates ordination objects for plotting later

ordinate <- function(obj){
  obj_mmds <- metaMDS(obj$data,
                      k =5,
                      distance = "bray"
  )


  obj_prcomp <- prcomp(obj$data,
                       scale =T,
                       center = T
  )

  obj_mmds <- as.data.frame(obj_mmds$points)
  obj_mmds <- obj_mmds[rownames(obj$meta),]
  obj$mds <- obj_mmds

  obj_prcomp <- as.data.frame(obj_prcomp$x[,1:5])
  obj_prcomp <- obj_prcomp[rownames(obj$meta),]
  obj$prcomp <- obj_prcomp

  return(obj)
}

###
#             Function run_beta                  #
###

#calculates beta-div between all samples

run_beta <- function(obj){
  bdiv <- as.matrix(vegdist(obj$data))
  obj$bdiv <- bdiv
  return(obj)
}

###
#        Function phylotype_analysis             #
###

phylotype_analysis <- function(obj, tax){
  #obj: microbiome object with at least 1 slot (data)
  #tax: a tax object (named list taxa as names values in the list are seq ids)
  obj.out <- NULL
  for(h in 1:length(tax)){
    df <- NULL
    for( i in 1:length(tax[[h]])){
      v1       <- obj$data[,unlist(tax[[h]][[i]])]
      v2       <- names(tax[[h]])[i]
      if(is.null(dim(v1))){
        df[[v2]] <- v1
      }else{
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

group_means <- function(df, mapping, t = F){
  #df : data frame with rows as sample IDs and columns as objects (e.g. OTUs)
  #mapping : a mapping df with rownames df as rownames and group id in col2
  #t : boolean, defaults to FALSE. Use if df has sample ids as colnames
  if(t == T){
    df <- t(df)
  }
  groups <- base::unique(x=mapping[,2])
  my_df <- data.frame(matrix(nrow = length(groups), ncol = ncol(df)))
  for(i in 1:ncol(df)){
    tgvec <- NULL
    for(j in 1:length(groups)){
      s <- base::rownames(mapping[base::which(mapping[,2] %in% groups[j]),
                                  ,drop=F])
      m <- base::mean(df[s,i])
      tgvec <- c(tgvec, m)
    }
    my_df[,i] <- tgvec
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

polyp2_asv.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/rarefied_10000_adenoma_asv.txt",
                            sep = "\t",
                            row.names = 1,
                            header = T)
polyp2_asv_rel.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/rel_abd_adenoma_asv.txt",
                                sep = "\t",
                                row.names = 1,
                                header = T)
polyp2_tax.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_adenoma_tax.txt",
                            sep = "\t",
                            row.names = 1,
                            header = T)
polyp2_metadata.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_combined_metadata.txt",
                                 sep = "\t",
                                 row.names = 1,
                                 header = T)

polyp2_sequence_stats.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/prepped_data/filtered_adenoma_sequence_stats.txt",
                                 sep = "\t",
                                 row.names = 1,
                                 header = T)

# ANALYSIS: READ COUNT SUMMARY STATS --------------------------------------

#just some quick summary stats

sum(polyp2_sequence_stats.df$reads.in)
#[1] 45160707
mean(polyp2_sequence_stats.df$reads.in)
#[1] 69692.45
median(polyp2_sequence_stats.df$reads.in)
#[1] 68883

sum(polyp2_sequence_stats.df$reads.out)
#[1] 33206620
mean(polyp2_sequence_stats.df$reads.out)
#[1] 51244.78
median(polyp2_sequence_stats.df$reads.out)
#[1] 50533.5

# ANALYSIS: FILTERED READ COUNT SUMMARY STATS ---------------------------------

# remove samples filtered for rarefaction
filtered_polyp2_sequence_stats.df <- polyp2_sequence_stats.df[which(rownames(
  polyp2_sequence_stats.df) %in% rownames(polyp2_metadata.df)),]

filtered_polyp2_sequence_stats.df$type <- polyp2_metadata.df$type

mfiltered_polyp2_sequence_stats.df <- melt(filtered_polyp2_sequence_stats.df)

reads_by_type.hist <- ggplot(data = mfiltered_polyp2_sequence_stats.df[
  which(mfiltered_polyp2_sequence_stats.df$variable == "reads.out"),],
       aes(
           x = value,
           fill = type)
          )

pdf("/Users/cgaulke/Documents/research/ohsu_polyp_combined/analysis/figs/qcreads_by_tissue.pdf")
reads_by_type.hist +
  geom_histogram(alpha = .5, position = "identity", color = "black")+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )+
  scale_x_continuous(expand = c(.02,0))+
  scale_y_continuous(expand = c(.01,0))+
  ylab("Count")+
  xlab("Read Depth")+
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

lc.vec    <-  polyp2_metadata.df[polyp2_metadata.df$location == "LC","id"]
rc.vec    <-  polyp2_metadata.df[polyp2_metadata.df$location == "RC","id"]
rec.vec   <-  polyp2_metadata.df[polyp2_metadata.df$location == "REC","id"]
fecal.vec <-  polyp2_metadata.df[polyp2_metadata.df$location == "Fecal","id"]
oral.vec  <-  polyp2_metadata.df[polyp2_metadata.df$location == "Oral","id"]

#only the three gut tissues
complete_gut.cases <- unique(rec.vec[which(rec.vec %in% lc.vec & rec.vec %in% rc.vec)])

#three gut tissues plus fecal and oral
complete.cases <- unique(rec.vec[which(rec.vec %in% lc.vec &
                                         rec.vec %in% rc.vec &
                                          rec.vec %in% fecal.vec &
                                           rec.vec %in% oral.vec)]
                         )

# DATA: MAKE OBJECT -------------------------------------------------------

polyp2_obj <- NULL
polyp2_obj$data <- polyp2_asv.df
polyp2_obj$meta <- polyp2_metadata.df
polyp2_obj$rel_abd <- polyp2_asv_rel.df

polyp2_obj$data <- polyp2_obj$data[rownames(polyp2_obj$meta),]

# DATA: AGGREGATE PHYLOTYPES -----------------------------------------------------

polyp2_taxadf <- polyp2_tax.df[which(rownames(polyp2_tax.df) %in% colnames(polyp2_obj$data)),
                               ,drop=F]

kingdom.df <- replicate(length(unique(polyp2_taxadf[,2])), c())
names(kingdom.df) <- unique(polyp2_taxadf[,2])
phylum.df  <- replicate(length(unique(polyp2_taxadf[,3])), c())
names(phylum.df) <- unique(polyp2_taxadf[,3])
class.df   <- replicate(length(unique(polyp2_taxadf[,4])), c())
names(class.df) <- unique(polyp2_taxadf[,4])
order.df   <- replicate(length(unique(polyp2_taxadf[,5])), c())
names(order.df) <- unique(polyp2_taxadf[,5])
family.df  <- replicate(length(unique(polyp2_taxadf[,6])), c())
names(family.df) <- unique(polyp2_taxadf[,6])
genus.df   <- replicate(length(unique(polyp2_taxadf[,7])), c())
names(genus.df) <- unique(polyp2_taxadf[,7])

for(i in 1:nrow(polyp2_taxadf)){

  kingdom.df[[polyp2_taxadf[i, 2]]] <- c(kingdom.df[[polyp2_taxadf[i, 2]]], polyp2_taxadf[i,1])
  phylum.df[[polyp2_taxadf[i, 3]]]  <- c(phylum.df[[polyp2_taxadf[i, 3]]], polyp2_taxadf[i,1])
  class.df[[polyp2_taxadf[i, 4]]]   <- c(class.df[[polyp2_taxadf[i, 4]]], polyp2_taxadf[i,1])
  order.df[[polyp2_taxadf[i, 5]]]   <- c(order.df[[polyp2_taxadf[i, 5]]], polyp2_taxadf[i,1])
  family.df[[polyp2_taxadf[i, 6]]]  <- c(family.df[[polyp2_taxadf[i, 6]]], polyp2_taxadf[i,1])
  genus.df[[polyp2_taxadf[i, 7]]]   <- c(genus.df[[polyp2_taxadf[i, 7]]], polyp2_taxadf[i,1])

}

tax.obj <- NULL
tax.obj$kingdom <- kingdom.df
tax.obj$phylum  <- phylum.df
tax.obj$class   <- class.df
tax.obj$order   <- order.df
tax.obj$family  <- family.df
tax.obj$genus   <- genus.df

#aggregate phylotype counts (not really a df, actually an obj)
polyp2_obj$phylotype <- phylotype_analysis(polyp2_obj,tax = tax.obj)

# ANALYSIS: DIVERSITY --------------------------------------------------

set.seed(731)
polyp2_obj <- diversity_analysis(polyp2_obj)
polyp2_obj <- run_beta(polyp2_obj)
polyp2_obj <- ordinate(polyp2_obj)

# ANALYSIS: SHANNON -----------------------------------------------

polyp2_obj$shannon.df <- cbind(polyp2_obj$shannon,
                               polyp2_obj$meta[names(polyp2_obj$shannon),"type"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"location"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"polyp"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"Polyp_N"])

colnames(polyp2_obj$shannon.df) <- c("shannon", "tissue", "location", "polyp", "npolyp")
polyp2_obj$shannon.df <- as.data.frame(polyp2_obj$shannon.df)
polyp2_obj$shannon.df$shannon <- as.numeric(polyp2_obj$shannon.df$shannon)

polyp2_shannon.boxplot <- ggplot(na.omit(polyp2_obj$shannon.df),
                                 aes(x = tissue,
                                     y = shannon,
                                     fill = polyp)
)

#pdf("/Users/cgaulke/Documents/research/ohsu_polyp2/analysis/figs/shannon_tissue_polyp.pdf")
polyp2_shannon.boxplot +
  geom_boxplot()+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )+
  ylab("Shannon")+
  xlab("")+
  scale_fill_brewer(palette = "Dark2")
#dev.off()

