# Title and author information --------------------------------------------
#!/usr/bin/R

#################################################
#                                               #
#    2020_11_19_polyp_combined_file_maker.R     #
#                                               #
#################################################

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

#Purpose: The purpose of this script is to prep the raw data file in this project
#for down stream analysis.

# Functions ---------------------------------------------------------------

###
#           Function adjust_colnames             #
###

#adjust colnames of dada2 sequence table to be colnames, not nt sequences

adjust_colnames <- function(df){
  x <- paste0(rep("seq", times = ncol(df)), c(1:ncol(df)))
  colnames(df) <- x
  return(df)
}

###
#             Function make_dict                 #
###

#This function makes a lookup table which maps IDs
#generated as in adjust_colnames with the sequence variant names

make_dict <- function(df){
  x <- cbind(id = paste0(rep("seq", times = ncol(df)), c(1:ncol(df))),
             sequence = colnames(df)
  )
  return(x)
}

###
#              Function filter_df                #
###

filter <- function(df) {
  df <- df[which(rowSums(df) > 0), which(colSums(df) > 0)]
  return(df)
}

###
#             Function normalize                 #
###

#normalize counts either by relative abundance of rarefying
#depends on vegan

normalize <- function(df, method="rel", depth=depth){
  #default method = relative abundance
  if(method == "rare"){
    if( is.null(depth)){
      depth <- min(rowSums(df))
    }
    ndf <- rrarefy(df, depth)
    ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }else{
    ndf <- sweep(df,
                 1,
                 rowSums(df),
                 `/`)
    #  ndf <- ndf[,which(colSums(ndf) > 0), drop =F]
  }
  return(ndf)
}


# SET ENVIRONMENT ---------------------------------------------------------

library(vegan)
options(stringsAsFactors = F)

# IMPORT: DATA ------------------------------------------------------------
#import sequence stats

adenoma_sequence_stats.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/polyp_combined_out/logs/filter_stats.txt",
                                 sep = "\t",
                                 row.names = 1,
                                 header = T)


#import data: ASV table
adenoma_asv.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/polyp_combined_out/sequence_table.txt",
                            sep = "\t",
                            row.names = 1,
                            header = T)

#import data: tax table

adenoma_tax.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/data/polyp_combined_out/tax_final.txt",
                           sep = "\t",
                           header = T,
                           row.names = NULL)

# Because these data are combined it becomes a little tricky once we get to the
# metadata so I will try to document everything I do here in good detail.

# we need to start by going back to the original metadata document that was
# provided by OHSU.


#import metadata

all_metadata.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/metadata/Microbiome_SharedData02132020_edit.txt",
                                 sep = "\t",
                                 header = T,
                                 row.names = 1)
# DATA: Filter Data -------------------------------------------------------

rownames(all_metadata.df) <- sapply(rownames(all_metadata.df),FUN = function(x) {(strsplit(x, split = "-")[[1]][2])})



#gather information for filtering data
adenoma_asv.names <- rownames(adenoma_asv.df)
adenoma.shortnames <- sapply(rownames(adenoma_asv.df),
                        FUN = function(x) { paste(strsplit(x, "-")[[1]][-(1:4)], collapse = "-")})

names.df <- cbind(adenoma_asv.names,adenoma.shortnames)

# since there are some pesky interlopers we will remove those completely now. We
# will also remove the biopsies that were sequenced in the first run and then
# again in run 2 or 3 and also the kit and reagent blanks because in the prior
# analyses these failed to reach thresholds at various steps indicating that
# contamination was not likely an issue.

#rownames(names.df[grep("B-", names.df[,2]),])
#rownames(names.df[grep(pattern = "blank", x = names.df[,2], ignore.case = T),])
#rownames(names.df[grep(pattern = "kit", x = names.df[,2], ignore.case = T),])


filter.names <- c( rownames(names.df[grep("B-",
                                          names.df[,2]),]) , #biopsy
                   rownames(names.df[grep(pattern = "blank",
                                          x = names.df[,2],
                                          ignore.case = T),]) , #blanks
                   rownames(names.df[grep(pattern = "kit",
                                          x = names.df[,2],
                                          ignore.case = T),]) , #kit blanks
                   rownames(names.df[which(names.df[,2] %in% c("3086",
                  "3087",
                  "3088",
                  "3089",
                  "3090",
                  "3091",
                  "3092",
                  "3093a",
                  "3093b",
                  "3094",
                  "9095",
                  "9096",
                  "9097",
                  "3100",
                  "3101",
                  "3102",
                  "3105a",
                  "3105b",
                  "P1V4",
                  "P2V5",
                  "P3V6",
                  "P44-200nm",
                  "P55-200nm",
                  "P66-20nm",
                  "P74-20nm",
                  "P85-20nm",
                  "P96-200nm")
),] ) # interlopers
)


adenoma_asv.df <- adenoma_asv.df[-which(rownames(adenoma_asv.df) %in% filter.names),]
names.df <- names.df[-which(rownames(names.df) %in% filter.names),]

#now we need to switch the F-sample# and S-sample# to be consistent with the
#rest of the data which is sample#-location

names.df[grep("S-", names.df[,2]),2] <- paste0(gsub(pattern = "S-", replacement = "", x = names.df[grep("S-", names.df[,2]),2]), "-S")
names.df[grep("F-", names.df[,2]),2] <- paste0(gsub(pattern = "F-", replacement = "", x = names.df[grep("F-", names.df[,2]),2]), "-F")

names.df <- as.data.frame(names.df)

names.df$id <- sapply(names.df[,2], FUN = function(x) {(strsplit(x, "-"))[[1]][1]})

names.df$type <- "Mucosal"

names.df$type[grep("-F", names.df$adenoma.shortnames)] <- "Fecal"
names.df$type[grep("-S$", names.df$adenoma.shortnames)] <- "Oral"


# if there is a better way to accomplish what I have done here I don't want to
# know

combined_metadata <- matrix(ncol = ncol(all_metadata.df), nrow = nrow(names.df))

for(i in 1:nrow(names.df)){
  combined_metadata[i,] <- unlist(all_metadata.df[names.df[i,3],])

}

#now add the row and column names
rownames(combined_metadata) <- names.df[,1]
colnames(combined_metadata) <- colnames(all_metadata.df)
combined_metadata <- as.data.frame(combined_metadata)

#add a patient ID
combined_metadata$id <- names.df$id

#add short sample name

combined_metadata$shortname <- names.df$adenoma.shortnames
combined_metadata$type <- names.df$type


#Now we have to add in some data for the mucosal samples because the have some
#unique characteristics

adenoma_metadata.df <- read.table("/Users/cgaulke/Documents/research/ohsu_polyp_combined/metadata/OHSU_polyp2_basic_metadata_filtered.txt",
                                 sep = "\t",
                                 header = T,
                                 row.names = NULL)



adenoma_metadata.df$id <- gsub(pattern = "_", replacement = "-",adenoma_metadata.df$Sample)

afilter <- c("3086",
             "3087",
             "3088",
             "3089",
             "3090",
             "3091",
             "3092",
             "3093a",
             "3093b",
             "3094",
             "9095",
             "9096",
             "9097",
             "3100",
             "3101",
             "3102",
             "3105a",
             "3105b",
             "P1V4",
             "P2V5",
             "P3V6",
             "P44-200nm",
             "P55-200nm",
             "P66-20nm",
             "P74-20nm",
             "P85-20nm",
             "P96-200nm",
             "Blank1",
             "Blank7",
             "KitB1",
             "KitB2")

adenoma_metadata.df <- adenoma_metadata.df[-which(adenoma_metadata.df$id %in% afilter),]

combined_metadata$location <- "Oral"
combined_metadata$location[which(combined_metadata$type == "Fecal")] <- "Fecal"
combined_metadata$polyp.tissue <- 0
combined_metadata$file_name <- rownames(combined_metadata)

#some sample have multiple polyp samples from the same part of the colon. Some
# are not labeled so we will label here.

combined_metadata$shortname[492] <- "106-P-LC-2"
adenoma_metadata.df$id[399] <- "106-P-LC-2"
rownames(combined_metadata) <- combined_metadata$shortname
rownames(adenoma_metadata.df) <- adenoma_metadata.df$id

#add the final bit of metadata
for(i  in 1:length(adenoma_metadata.df$id)){
  combined_metadata[adenoma_metadata.df$id[i],"polyp.tissue"] <-
    adenoma_metadata.df[adenoma_metadata.df$id[i], "Polyp.tissue"]
  combined_metadata[adenoma_metadata.df$id[i],"location"] <-
    adenoma_metadata.df[adenoma_metadata.df$id[i], "location"]
}

rownames(combined_metadata) <- combined_metadata$file_name

#now the metadata is ready to play nice with the rest of the data

# the last thing we need to do is filter the sequence statistics to remove unused
# sample data

rownames(adenoma_sequence_stats.df) <- sapply(rownames(
  adenoma_sequence_stats.df),
    FUN = function(x) {(strsplit(x, split = "_")[[1]][1])})

adenoma_sequence_stats.df <- adenoma_sequence_stats.df[which(
  rownames(adenoma_sequence_stats.df) %in% rownames(combined_metadata)),]

# DATA WRANGLE: FIX FILE HEADERS ----------------------------------------------

# One of the "features" of dada2 is the use of very long sequence variant names
# (i.e., the nt sequence itself). I fix this by assigning a sequential name to
# each asv (e.g., seq1, seq2, seq3, etc.), however I keep the information about
# the sequence in a dictionary that can be exported and used later.

# lets start with the asv table. We know that the colnames of the asv table
# should match the polyp2_tax.df$row.names, but it never hurts to be thorough.
all(adenoma_tax.df$row.names == colnames(adenoma_asv.df))

#fix dada2 colnames
adenoma_asv.dict <- make_dict(adenoma_asv.df) # we keep a record just in case
adenoma_asv.df <- adjust_colnames(adenoma_asv.df)

#double check to make sure we haven't mangled anything
colnames(adenoma_asv.df)

#Good. Now onto the taxa data frame

rownames(adenoma_tax.df) <-
  paste0("seq", seq(from = 1, to = nrow(adenoma_tax.df)))
adenoma_tax.df$row.names <- rownames(adenoma_tax.df)
colnames(adenoma_tax.df)[1] <- "id"

# Next I will collect all of the sequences that are annotated to Eukaryota.
# Given that the V4 primers used in this experiment are specific for the
# classification of bacteria it is unclear if the region they would amplify
# from Euks would be taxonomically informative.Looking at the data this seems
# to be the case. I will remove them before rarefaction so it does not impact
#downstream analysis as much.

euks <- rownames(adenoma_tax.df[which(adenoma_tax.df$Kingdom =="Eukaryota"),])


# DATA TRANSFORMATION: FILTER AND NORMALIZE -------------------------------

# First we remove these Euk ASVs from the table for the reasons outlined above.

adenoma_asv.df  <- adenoma_asv.df[,-which(colnames(adenoma_asv.df) %in% euks),drop=F]

#get rid of 0 sum columns (i.e., those with no counts)
adenoma_asv.df <- filter(adenoma_asv.df)
# To select an appropriate rarefaction depth I will generate some rarefaction
# plots first

#n species in each sample
nspec <- specnumber(adenoma_asv.df)

#Check spread
hist(nspec)

# get read depth

nrsums <- rowSums(adenoma_asv.df)

#Check spread
hist(nrsums)

#min number observations in the df
raresamp <- min(rowSums(adenoma_asv.df))

#now we can plot the curves
rarecurve(adenoma_asv.df, step = 1000,
          sample = raresamp,
          label = F,
          col = "lightblue")

abline(v = 500, col = "red")
abline(v = 1000, col = "red")
abline(v = 2500, col = "red")
abline(v = 5000, col = "yellow")
abline(v = 7500, col = "yellow")
abline(v = 10000, col = "black")

# so it looks like a read depth of 10,000 sequences is about appropriate for our
# data. Lower depths such as 7,500 and 5,000 might also be acceptable, however,
# the intercept of the vline is more frequently in the increasing phase of the
# curve with these depths. A depth of 10,000 reads avoids this but also will
# filter more sequences. We will err on the side of caution here.

#Gather the samples with < 10000 reads and remove them. This is about 4% of
#samples (25)

filter.names <- names(rowSums(adenoma_asv.df)[rowSums(adenoma_asv.df) < 10000])

# I am removing patient # 56 because we don't have metadata for this individual
filter.names <- c(filter.names,"lane1-s208-index-CCTAACGGTCCA-056-RC",
                  "lane1-s183-index-CGAATACTGACA-S-056" )

adenoma_asv.df <- adenoma_asv.df[-which(rownames(adenoma_asv.df) %in% filter.names), ]

# now sync asv table with metadata

combined_metadata <- combined_metadata[which(rownames(combined_metadata) %in% rownames(adenoma_asv.df)), ]

#finally filter and rarefy the asv table
adenoma_asv.filt <- filter(adenoma_asv.df)

#now lets check the read depth, # species, and rarefaction curves again.

nspec <- specnumber(adenoma_asv.df)

pdf("figs/combined_nspecies_hist.pdf")
hist(nspec)
dev.off()
# get read depth

nrsums <- rowSums(adenoma_asv.df)

pdf("figs/read_depth_hist.pdf")
hist(nrsums)
dev.off()

#now we can plot the curves
pdf("figs/rarecurves.pdf")
rarecurve(adenoma_asv.df, step = 1000,
          sample = 0,
          label = F,
          col = "lightblue")

abline(v = 500, col = "red")
abline(v = 1000, col = "red")
abline(v = 2500, col = "red")
abline(v = 5000, col = "yellow")
abline(v = 7500, col = "yellow")
abline(v = 10000, col = "black")
dev.off()

#now move on to normalization and filter

set.seed(731)
#normalize data
rare_adenoma_asv <- normalize(adenoma_asv.filt,
                                 method = "rare",
                                 depth = 10000)

rel_adenoma_asv <- normalize(adenoma_asv.filt,
                                method = "rel",
                                depth = 10000)

#filter again
rare_adenoma_asv <- filter(rare_adenoma_asv)
rel_adenoma_asv <- filter(rel_adenoma_asv)

#all done here


# DATA WRANGLE: RECODE METADATA -------------------------------------------

# There is inconsistent labeling of location used in this metadata. We need to
# clean this up a bit

unique(combined_metadata$location)
# all polyp tissue will be denoted with a "P_"
#location will be assigned regardless of polyp status

lc <- c("P_Ds", "P_RS", "LC", "P_LC-2x", "P_LC", "P_SF1", "P_SF2", "LC", "P_LC-2x", "P_LC")
rc <- c("P_HF2", "P_As", "P_Ce3", "RC", "P_Ce4", "P_Ce", "P_Rc", "P_RC1", "P_As1", "P_RC", "P_Ce1", "P_RC2", "P_As2", "P_HF1", "P_Ce2", "P_HF","P_RC3")
tv <-  c("P_Tv", "P_TV1","TVnormal","P_ProxTV", "P_TV2", "PTV", "P_TV")
sgd <- c("P_Sg", "P_Sg1", "Sig", "P_Sg2", "P-Sg")
rec <- c("P_Re1", "P_Re2", "Re", "P_Re", "Re-duplicate")


combined_metadata[which(combined_metadata$location %in% lc), "location"] <- "LC"
combined_metadata[which(combined_metadata$location %in% rc), "location"] <- "RC"
combined_metadata[which(combined_metadata$location %in% tv), "location"] <- "TV"
combined_metadata[which(combined_metadata$location %in% sgd), "location"] <- "SGD"
combined_metadata[which(combined_metadata$location %in% rec), "location"] <- "REC"
combined_metadata$polyp = as.numeric(((combined_metadata$Polyp_N) > 0))

#final filter to remove samples that meet exclusion criteria

#antibiotics use == 1

combined_metadata <- combined_metadata[-which(combined_metadata$AntiBio == 1 ),]

#polyps but no adenomas
combined_metadata <- combined_metadata[-which(combined_metadata$Polyp_N > 0 & combined_metadata$Adenoma == 0),]


#now the final filtering of the other files
rare_adenoma_asv <- rare_adenoma_asv[which(rownames(rare_adenoma_asv) %in% rownames(combined_metadata)),]
rare_adenoma_asv <- filter(rare_adenoma_asv)

rel_adenoma_asv <- rel_adenoma_asv[which(rownames(rel_adenoma_asv) %in% rownames(combined_metadata)),]
rel_adenoma_asv <- filter(rel_adenoma_asv)

adenoma_sequence_stats.df <- adenoma_sequence_stats.df[which(rownames(adenoma_sequence_stats.df) %in% rownames(combined_metadata)),]

# EXPORT DATA  ------------------------------------------------------------

# The code below is run only once to prevent inadvertent modifying these files
# later on.
#
#
# write.table(x = rare_adenoma_asv,
#             file = "data/prepped_data/rarefied_10000_adenoma_asv.txt",
#             quote = F,
#             row.names = T,
#             col.names = T,
#             sep = "\t")
#
#
# write.table(x = rel_adenoma_asv,
#             file = "data/prepped_data/rel_abd_adenoma_asv.txt",
#             quote = F,
#             row.names = T,
#             col.names = T,
#             sep = "\t")
#
# write.table(x = combined_metadata,
#            file = "data/prepped_data/filtered_combined_metadata.txt",
#            quote = F,
#            row.names = T,
#            col.names = T,
#            sep = "\t")
#
# write.table(x = adenoma_tax.df,
#             file = "data/prepped_data/filtered_adenoma_tax.txt",
#             quote = F,
#             row.names = T,
#             col.names = T,
#             sep = "\t")
#
# write.table(x = adenoma_sequence_stats.df,
#             file = "data/prepped_data/filtered_adenoma_sequence_stats.txt",
#             quote = F,
#             row.names = T,
#             col.names = T,
#             sep = "\t")
