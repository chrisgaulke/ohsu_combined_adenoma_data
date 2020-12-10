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
#library(lme4)
#library(lmerTest)
library(MASS)
library(glmmTMB)

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

  obj$mds  <- obj_mmds
  obj_mmds <- as.data.frame(obj_mmds$points)
  obj_mmds <- obj_mmds[rownames(obj$meta),]
  obj$mds.df <- obj_mmds


  obj$prcomp <- obj_prcomp
  obj_prcomp    <- as.data.frame(obj_prcomp$x[,1:5])
  obj_prcomp    <- obj_prcomp[rownames(obj$meta),]
  obj$prcomp.df <- obj_prcomp


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

pdf("figs/qcreads_by_tissue.pdf")
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

# ANALYSIS: COMPUTE DIVERSITY ---------------------------------------------

set.seed(731)
polyp2_obj <- diversity_analysis(polyp2_obj)
polyp2_obj <- run_beta(polyp2_obj)
polyp2_obj <- ordinate(polyp2_obj)

# ANALYSIS: SHANNON ------------------------------------------------------

polyp2_obj$shannon.df <- cbind(polyp2_obj$shannon,
                               polyp2_obj$meta[names(polyp2_obj$shannon),"type"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"location"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"polyp"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"Adenoma"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"polyp.tissue"],
                               polyp2_obj$meta[names(polyp2_obj$shannon),"id"])

colnames(polyp2_obj$shannon.df) <- c("shannon", "tissue", "location", "former", "npolyp", "polyptissue","id")
polyp2_obj$shannon.df <- as.data.frame(polyp2_obj$shannon.df)
polyp2_obj$shannon.df$shannon <- as.numeric(polyp2_obj$shannon.df$shannon)
polyp2_obj$shannon.df$npolyp <- as.numeric(polyp2_obj$shannon.df$npolyp)

polyp2_shannon.boxplot <- ggplot(na.omit(polyp2_obj$shannon.df),
                                 aes(x = tissue,
                                     y = shannon,
                                     fill = former)
)

pdf("figs/shannon_tissue_polyp.pdf")
polyp2_shannon.boxplot +
  geom_boxplot()+
  geom_point(position = position_dodge(width= .75), color = "black", shape = 21, alpha = .5 )+
  theme(text = element_text(size=24, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1.5,
        legend.position = "top"
  )+
  ylab("Shannon")+
  xlab("")+
  scale_fill_brewer("", palette = "Dark2", labels = c("Non-Former", "Former"))+
  scale_color_brewer("", palette = "Dark2")
dev.off()

polyp2_shannon.boxplot <- ggplot(na.omit(polyp2_obj$shannon.df),
                                 aes(x = location,
                                     y = shannon,
                                     fill = former)
)

pdf("figs/shannon_location_boxplot.pdf")
polyp2_shannon.boxplot +
  geom_boxplot()+
  geom_point(position = position_dodge(width= .75), color = "black", shape = 21, alpha = .5 )+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )+
  ylab("Shannon")+
  xlab("")+
  scale_fill_brewer("Former", palette = "Dark2")+
  scale_color_brewer("Former", palette = "Dark2")
dev.off()

#this is the whole enchilada here. It might just be me but I find this really
#difficult to determine if this actually does what I want here.

lmefit1 <- lmer(shannon ~ npolyp *
                  factor(tissue) + (1 + tissue|id) ,
                na.omit(polyp2_obj$shannon.df))
anova(lmefit1)
summary(lmefit1) #sig

#mucosal former effects

lmefit2 <- lmer(shannon ~
                  npolyp + (1 |id),
                na.omit(subset(polyp2_obj$shannon.df,subset = tissue=="Mucosal")))
summary(lmefit2) # no sig

#mucosal adenomas with location effects

lmefit3 <- lmer(shannon ~
                  npolyp + location + (1 |id),
                na.omit(subset(polyp2_obj$shannon.df,subset = tissue=="Mucosal")))
anova(lmefit3)
summary(lmefit3) # no sig


#fecal npolyp
lmfit1 <- lm(shannon ~
               npolyp ,
             na.omit(subset(polyp2_obj$shannon.df,subset = tissue=="Fecal")))
summary(lmfit1) #sig


#oral npolyp
#(Note that these models only include main effects bc there aren't multi samples)

lmfit2 <- lm(shannon ~
                npolyp ,
             na.omit(subset(polyp2_obj$shannon.df,subset = tissue=="Oral")))
summary(lmfit2) #no sig


adenoma_shannon.point <- ggplot(na.omit(polyp2_obj$shannon.df),
                                aes(x = npolyp,
                                    y = shannon,
                                    color = tissue)
                                )
pdf("figs/shannon_regression.pdf",
   width = 14, height = 7)

adenoma_shannon.point +
  geom_point(size = 3, alpha = .7 ) +
  geom_smooth(method = "glm", aes(fill = tissue), alpha = .1)+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1
  )+
  ylab("Shannon")+
  xlab("Adenomas")+
  facet_wrap(.~tissue)+
  scale_color_brewer("Tissue",palette = "Dark2")+
  scale_fill_brewer("Tissue",palette = "Dark2")

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
                                   fill = former
                               )
)

pdf("figs/richness_tissue_polyp.pdf")
polyp2_richness.plot +
  geom_boxplot()+
  geom_point(position = position_dodge(width= .75), color = "black", shape = 21, alpha = .5 )+
  theme(text = element_text(size=24, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        aspect.ratio = 1.5,
        legend.position = "top"
  )+
  ylab("Richness")+
  xlab("")+
  scale_fill_brewer("", palette = "Dark2", labels = c("Non-Former", "Former"))+
  scale_color_brewer("", palette = "Dark2")
dev.off()


polyp2_richness.boxplot <- ggplot(na.omit(polyp2_obj$richness),
                                 aes(x = location,
                                     y = richness,
                                     fill = former)
)

pdf("figs/richness_location_boxplot.pdf")
polyp2_richness.boxplot +
  geom_boxplot()+
  geom_point(position = position_dodge(width= .75), color = "black", shape = 21, alpha = .5 )+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black")
  )+
  ylab("Richness")+
  xlab("")+
  scale_fill_brewer("Former", palette = "Dark2")+
  scale_color_brewer("Former", palette = "Dark2")
dev.off()


adenoma_richness.point <- ggplot(na.omit(polyp2_obj$richness),
                                aes(x = npolyp,
                                    y = richness,
                                    color = tissue)
)

# So I think that what we can say here is that the there are tissue level differences
# in richness (no surprise here). The relationship with number of polyps is
# interesting but I am not sure exactly how to interpret the tissue level
# interaction terms.

lmerich_fit1 <- lmer(richness ~ npolyp * factor(tissue) +
                (1 +tissue|id) ,na.omit(polyp2_obj$richness))
anova(lmerich_fit1)
summary(lmerich_fit1)#sig

#mucosal former effects
lmerich_fit2 <- lmer(richness ~
                       npolyp + (1|id),
                     na.omit(subset(polyp2_obj$richness,subset = tissue=="Mucosal")))
summary(lmerich_fit2)#no sig

#fecal former
#(Note that these models only include main effects bc there aren't multi samples)
lmrich_fit1 <- lm(richness ~
                    npolyp ,
                  na.omit(subset(polyp2_obj$richness,subset = tissue=="Fecal")))
summary(lmrich_fit1)#sig

#oral former
#(Note that these models only include main effects bc there aren't multi samples)

lmrich_fit2 <- lm(richness ~
                    npolyp ,
                  na.omit(subset(polyp2_obj$richness,subset = tissue=="Oral")))
summary(lmrich_fit2) #no sig


#mucosal adenomas with location effects

lmerichfit3 <- lmer(richness ~
                  npolyp + location + (1 |id),
                na.omit(subset(polyp2_obj$richness,subset = tissue=="Mucosal")))
anova(lmerichfit3)
summary(lmerichfit3) # marginal sig

# These results mirror those from the original (separate) analyses which indicate
# that only fecal alpha diversity associates with number of polyps


pdf("figs/richness_regression.pdf", width = 14, height = 7)

adenoma_richness.point +
  geom_point(size = 3, alpha = .7 ) +
  geom_smooth(method = "glm", aes(fill = tissue), alpha = .1)+
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1
  )+
  ylab("Richness")+
  xlab("Adenomas")+
  facet_wrap(.~tissue)+
  scale_color_brewer("Tissue",palette = "Dark2")+
  scale_fill_brewer("Tissue",palette = "Dark2")

dev.off()


# ANALYSIS: BETA DIVERSITY NMDS ---------------------------------------------

#Start with a simple all tissue ordination

#add metadata
polyp2_obj$mds.df$tissue       <- polyp2_obj$meta$type
polyp2_obj$mds.df$former       <- polyp2_obj$meta$polyp
polyp2_obj$mds.df$npolyp       <- polyp2_obj$meta$Adenoma
polyp2_obj$mds.df$polyp.tissue <- polyp2_obj$meta$polyp.tissue
polyp2_obj$mds.df$location     <- polyp2_obj$meta$location
polyp2_obj$mds.df$nbin         <- cut(polyp2_obj$mds.df$npolyp, # for plotting
                                      breaks = c(-1, 0.1,2.1,5,10, 20 ),
                                      labels = c("0","1-2", "3-5","6-10","10+"))

tissue_mds.ord <- ggplot(na.omit(polyp2_obj$mds.df),
                            aes(x = MDS1,
                                y = MDS2,
                                color = tissue)
)

pdf("figs/tissue_nmds.pdf")
tissue_mds.ord +
  geom_point( size = 3, alpha = .4 )+
  theme(text = element_text(size=18, colour = "black"),
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
  scale_color_brewer("",palette = "Dark2")
dev.off()


mucosal_mds.ord <- ggplot(na.omit(subset(polyp2_obj$mds.df,subset = tissue == "Mucosal")),
                         aes(x = MDS1,
                             y = MDS2,
                             color = nbin
                             )
)

pdf("figs/mucosal_nmds.pdf")
mucosal_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size =3 ) +
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border     = element_rect(fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1
  )+
  scale_alpha_manual("", breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
                     values =  c(0.1, 0.7, 0.7, 0.7, .7))+
  scale_color_brewer("Adenomas", palette = "Set1")+
  guides(alpha = FALSE)
dev.off()

#fecal

fecal_mds.ord <- ggplot(na.omit(subset(polyp2_obj$mds.df,subset = tissue == "Fecal")),
                          aes(x = MDS1,
                              y = MDS2,
                              color = nbin
                          )
)

pdf("figs/fecal_nmds.pdf")
fecal_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size =3 ) +
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border     = element_rect(fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1
  )+
  scale_alpha_manual("", breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
                     values =  c(0.1, 0.7, 0.7, 0.7, .7))+
  scale_color_brewer("Adenomas", palette = "Set1")+
  guides(alpha = FALSE)
dev.off()



#oral

oral_mds.ord <- ggplot(na.omit(subset(polyp2_obj$mds.df,subset = tissue == "Oral")),
                        aes(x = MDS1,
                            y = MDS2,
                            color = nbin
                        )
)

pdf("figs/oral_nmds.pdf")
oral_mds.ord +
  geom_point(aes(alpha = as.character(nbin)),
             size =3 ) +
  theme(text = element_text(size=18, colour = "black"),
        panel.grid.major = element_line(color = "grey97"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border     = element_rect(fill = NA),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        aspect.ratio = 1
  )+
  scale_alpha_manual("", breaks = c("0", "1-2", "3-5" , "6-10", "10+"),
                     values =  c(0.1, 0.7, 0.7, 0.7, .7))+
  scale_color_brewer("Adenomas", palette = "Set1")+
  guides(alpha = FALSE)
dev.off()


# ANALYSIS: BETA DIVERSITY ADONIS (ASV) -------------------------------------

#start with one big model that includes all the data
set.seed(731)
tissue_npolyp.adonis <- adonis(polyp2_obj$data ~  polyp2_obj$meta$Adenoma * polyp2_obj$meta$type,
                                                    permutations = 5000,
                                                    method = "bray")
tissue_npolyp.adonis

#this model indicates that there are significant impacts of tissue and number of
#adenomas on the microbiome as well as an interaction between these two effects


#lets look at mucosa separately
set.seed(731)
mucosal_npolyp.adonis <-
  adonis(
    polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), ]),] ~
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
  adonis(polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), ]),] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

fecal_npolyp.adonis #not significant

#now oral
set.seed(731)

oral_npolyp.adonis <-
  adonis(polyp2_obj$data[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), ]),] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

oral_npolyp.adonis #not significant

# ANALYSIS: BETA DIVERSITY ADONIS (Genus) -------------------------------------

#start with one big model that includes all the data
set.seed(731)
gtissue_npolyp.adonis <- adonis(polyp2_obj$phylotype$genus ~  polyp2_obj$meta$Adenoma * polyp2_obj$meta$type,
                               permutations = 5000,
                               method = "bray")
gtissue_npolyp.adonis

#this model indicates that there are significant impacts of tissue and number of
#adenomas on the microbiome as well as an interaction between these two effects


#lets look at mucosa separately
set.seed(731)
gmucosal_npolyp.adonis <-
  adonis(
    polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"), ]),] ~
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
  adonis(polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), ]),] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

gfecal_npolyp.adonis #not significant

#now oral
set.seed(731)

goral_npolyp.adonis <-
  adonis(polyp2_obj$phylotype$genus[rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), ]),] ~
           polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"), "Adenoma"] ,
         permutations = 5000,
         method = "bray")

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
       rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),]
        )
  ),
  ]

#Filter
fecal_adenoma_asv.df <-
  fecal_adenoma_asv.df[,apply(fecal_adenoma_asv.df, 2,
                function(x){sum(as.logical(x)) >
                    round(nrow(fecal_adenoma_asv.df) * .2)})
  ]

#Fecal Genus
fecal_adenoma_genus.df <-
  polyp2_obj$phylotype$genus[which(rownames(polyp2_obj$phylotype$genus) %in%
           rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Fecal"),]
                                              )
                                   ),
                             ]

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
        rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),]
                )
  ),
  ]

#Filter
mucosal_adenoma_asv.df <-
  mucosal_adenoma_asv.df[,apply(mucosal_adenoma_asv.df, 2,
                      function(x){sum(as.logical(x)) >
                          round(nrow(mucosal_adenoma_asv.df) * .2)})
  ]

#MucosalGenus
mucosal_adenoma_genus.df <-
  polyp2_obj$phylotype$genus[which(rownames(polyp2_obj$phylotype$genus) %in%
      rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Mucosal"),]
               )
  ),
  ]

#Filter
mucosal_adenoma_genus.df <-
  mucosal_adenoma_genus.df[,apply(mucosal_adenoma_genus.df, 2,
              function(x){sum(as.logical(x)) >
                  round(nrow(mucosal_adenoma_genus.df) * .2)})
  ]


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
              rownames(polyp2_obj$meta[which(polyp2_obj$meta$type == "Oral"),]
                      )
  ),
  ]

#Filter

oral_adenoma_genus.df <-
  oral_adenoma_genus.df[,apply(oral_adenoma_genus.df, 2,
                      function(x){
                        sum(as.logical(x)) >
                          round(nrow(oral_adenoma_genus.df) * .2)})
]

#make master objects
asv.models.obj   <- NULL
genus.models.obj <- NULL

# ANALYSIS: TAXA MODELS MUCOSAL GENUS -------------------------------------

#nadenoma models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_genus.df[, i] ~ Adenoma + (1 |id),
             data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1);
  mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <- fit1}, error=function(e) e

  )
}

#assign to proper slot so we can recycle the object
genus.models.obj$nadenoma <- mucosal_adenoma_models

#former models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_genus.df[, i] ~ factor(polyp) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1)
    ;
  mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <- fit1}, error=function(e) e

  )
}

#assign to proper slot so we can recycle the object
genus.models.obj$former <- mucosal_adenoma_models

#polyp.tissue models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_genus.df[, i] ~ factor(polyp.tissue) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1)
  ;
  mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <- fit1}, error=function(e) e

  )
}

#assign to proper slot so we can recycle the object
genus.models.obj$polyp.tissue <- mucosal_adenoma_models

#location models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_genus.df[, i] ~ factor(location) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1)
  ;
  mucosal_adenoma_models[[colnames(mucosal_adenoma_genus.df)[i]]] <- fit1}, error=function(e) e

  )
}

#assign to proper slot so we can recycle the object
genus.models.obj$location <- mucosal_adenoma_models

#basic stats on fdr control genus
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$nadenoma,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$former,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$polyp.tissue,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$location,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))


# ANALYSIS: TAXA MODELS MUCOSAL ASV ---------------------------------------

#nadenoma models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_asv.df[, i] ~ Adenoma + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1);
  mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <- fit1},
    error=function(e) e

  )
}

asv.models.obj$nadenoma <- mucosal_adenoma_models

#former models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_asv.df[, i] ~ factor(polyp) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1);
  mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$former <- mucosal_adenoma_models

#polyp.tissue models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_asv.df[, i] ~ factor(polyp.tissue) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1);
  mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$polyp.tissue <- mucosal_adenoma_models

#location models
set.seed(731)
mucosal_adenoma_models <- NULL

for(i in 1:ncol(mucosal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glmmTMB(mucosal_adenoma_asv.df[, i] ~ factor(location) + (1 |id),
            data = subset(polyp2_obj$meta, type == "Mucosal"), family=nbinom1);
  mucosal_adenoma_models[[colnames(mucosal_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$location <- mucosal_adenoma_models

#basic stats on fdr control ASV
sum(na.omit(qvalue::qvalue(sapply(asv.models.obj$nadenoma,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(asv.models.obj$former,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(asv.models.obj$polyp.tissue,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(asv.models.obj$location,
                                  FUN = function(x){ coefficients(summary(x))$cond[2,4]}))$qvalue < .1 ))

# ANALYSIS: TAXA MODELS FECAL GENUS -------------------------------------

# Fecal and oral models will be a little different. This is because there is
# only 1 samples / ID and thus using ID as a random effect is probably not a
# good idea. I will just use a standard neg binomial glm for these models.

#fecal #adenoma

set.seed(731)
fecal_adenoma_models <- NULL

for(i in 1:ncol(fecal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(fecal_adenoma_genus.df[, i] ~ Adenoma ,
            data = subset(polyp2_obj$meta, type == "Fecal")
           );
  fecal_adenoma_models[[colnames(fecal_adenoma_genus.df)[i]]] <- fit1},
  error=function(e) e

  )
}

genus.models.obj$fecal_nadenoma <- fecal_adenoma_models

#fecal former
set.seed(731)
fecal_adenoma_models <- NULL

for(i in 1:ncol(fecal_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(fecal_adenoma_genus.df[, i] ~ factor(polyp) ,
           data = subset(polyp2_obj$meta, type == "Fecal")
    );
  fecal_adenoma_models[[colnames(fecal_adenoma_genus.df)[i]]] <- fit1},
  error=function(e) e

  )
}

genus.models.obj$fecal_former <- fecal_adenoma_models

# ANALYSIS: TAXA MODELS FECAL ASV ---------------------------------------

#fecal asv #adenoma

set.seed(731)
fecal_adenoma_models <- NULL

for(i in 1:ncol(fecal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(fecal_adenoma_asv.df[, i] ~ Adenoma ,
           data = subset(polyp2_obj$meta, type == "Fecal")
    );
  fecal_adenoma_models[[colnames(fecal_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$fecal_nadenoma <- fecal_adenoma_models


#fecal former
set.seed(731)
fecal_adenoma_models <- NULL

for(i in 1:ncol(fecal_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(fecal_adenoma_asv.df[, i] ~ factor(polyp) ,
           data = subset(polyp2_obj$meta, type == "Fecal")
    );
  fecal_adenoma_models[[colnames(fecal_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$fecal_former <- fecal_adenoma_models
fecal_adenoma_models <- NULL

# ANALYSIS: TAXA MODELS ORAL GENUS --------------------------------------

#oral #adenoma

set.seed(731)
oral_adenoma_models <- NULL

for(i in 1:ncol(oral_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(oral_adenoma_genus.df[, i] ~ Adenoma ,
           data = subset(polyp2_obj$meta, type == "Oral")
    );
  oral_adenoma_models[[colnames(oral_adenoma_genus.df)[i]]] <- fit1},
  error=function(e) e

  )
}

genus.models.obj$oral_nadenoma <- oral_adenoma_models

#oral former
set.seed(731)
oral_adenoma_models <- NULL

for(i in 1:ncol(oral_adenoma_genus.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(oral_adenoma_genus.df[, i] ~ factor(polyp) ,
           data = subset(polyp2_obj$meta, type == "Oral")
    );
  oral_adenoma_models[[colnames(oral_adenoma_genus.df)[i]]] <- fit1},
  error=function(e) e

  )
}

genus.models.obj$oral_former <- oral_adenoma_models

#basic stats on fdr control genus
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$fecal_nadenoma,
                                  FUN = function(x){ coefficients(summary(x))[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$fecal_former,
                                  FUN = function(x){ coefficients(summary(x))[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$oral_nadenoma,
                                  FUN = function(x){ coefficients(summary(x))[2,4]}))$qvalue < .1 ))
sum(na.omit(qvalue::qvalue(sapply(genus.models.obj$oral_former,
                                  FUN = function(x){ coefficients(summary(x))[2,4]}))$qvalue < .1 ))

# ANALYSIS: TAXA MODELS ORAL ASV ----------------------------------------

#oral #adenoma

set.seed(731)
oral_adenoma_models <- NULL

for(i in 1:ncol(oral_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(oral_adenoma_asv.df[, i] ~ Adenoma ,
           data = subset(polyp2_obj$meta, type == "Oral")
    );
  oral_adenoma_models[[colnames(oral_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$oral_nadenoma <- oral_adenoma_models


#oral former
set.seed(731)
oral_adenoma_models <- NULL

for(i in 1:ncol(oral_adenoma_asv.df)){
  fit1 <- NULL
  tryCatch({fit1 <-
    glm.nb(oral_adenoma_asv.df[, i] ~ factor(polyp) ,
           data = subset(polyp2_obj$meta, type == "Oral")
    );
  oral_adenoma_models[[colnames(oral_adenoma_asv.df)[i]]] <- fit1},
  error=function(e) e

  )
}

asv.models.obj$oral_former <- oral_adenoma_models
oral_adenoma_models<- NULL

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

# ANALYSIS: RANDOM FORESTS ------------------------------------------------

#Mucosal
mucosal_former_asv.rf <- randomForest(x = mucosal_adenoma_asv.df,
                                  y = factor(unlist(subset(polyp2_obj$meta,
                                                           type == "Mucosal",
                                                           select = polyp))))

mucosal_former_genus.rf <- randomForest(x = mucosal_adenoma_genus.df,
                                  y = factor(unlist(subset(polyp2_obj$meta,
                                                           type == "Mucosal",
                                                           select = polyp))))

#Fecal
fecal_former_asv.rf <- randomForest(x = fecal_adenoma_asv.df,
                                    y = factor(unlist(subset(polyp2_obj$meta,
                                                             type == "Fecal",
                                                             select = polyp))))

fecal_former_genus.rf <- randomForest(x = fecal_adenoma_genus.df,
                                    y = factor(unlist(subset(polyp2_obj$meta,
                                                             type == "Fecal",
                                                             select = polyp))))

#Oral
oral_former_asv.rf <- randomForest(x = oral_adenoma_asv.df,
                                    y = factor(unlist(subset(polyp2_obj$meta,
                                                             type == "Oral",
                                                             select = polyp))))

oral_former_genus.rf <- randomForest(x = oral_adenoma_genus.df,
                                    y = factor(unlist(subset(polyp2_obj$meta,
                                     type == "Oral", select = polyp))))


