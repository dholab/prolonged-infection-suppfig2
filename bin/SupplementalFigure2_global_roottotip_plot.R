#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING INTRAHOST SINGLE-NUCLEOTIDE VARIANT (iSNV) FREQUENCES
# UPDATED: 02-Aug-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will do the following:

# 1) Import VCFs called from GenBank subsample.
# 2) Import aligned FASTA sequences from GenBank subsample, plus their metadata.
# 3) Polish these three files.
# 4) Count single-nucleotide variations from Wuhan-1 (GenBank: MN908947.3)
#    present at each time point, as well as in each GenBank sample.
# 5) Plot these mutations through time.

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### PREPARING THE ENVIRONMENT ####
library(tidyverse)
data_filepath = args[1]  ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)

patient_day0 <- as.Date("2020-12-21")-113

if (length(list.files(getwd(), "*.vcf"))>0){
  
  print("Processing metadata.")
  
  vcf_list <- list.files(getwd(), "*.vcf")
  vcf_list <- sort(vcf_list)
  genbank_meta <- read.delim("b12_enriched_global_subsampled_metadata.tsv")
  genbank_meta <- genbank_meta[order(genbank_meta$strain),] ; rownames(genbank_meta) <- NULL
  genbank_meta$date <- as.Date(genbank_meta$date)
  distances <- data.frame("sample" = rep(NA, times = length(vcf_list)),
                          "distance" = rep(NA, times = length(vcf_list)),
                          "date" = rep(as.Date("2000-01-01"), times = length(vcf_list)),
                          "day" = rep(NA, times = length(vcf_list)),
                          "lineage" = rep(NA, times = length(vcf_list)))
  
  for (i in 1:length(vcf_list)){
    
    print(paste("Importing and processing", vcf_list[i], sep = " "))
    
    vcf <- read.delim(paste(getwd(), vcf_list[i], sep = "/"), skip = 55)
    vcf <- vcf[str_length(vcf$REF)==1,] # filtering to insertions and substitutions only
    vcf <- vcf[str_length(vcf$ALT)==1,] # filtering to substitutions only
    vcf <- vcf[vcf$QUAL>1,] # filtering to quality variants likely to be real
    sample <- colnames(vcf)[ncol(vcf)]
    sample <- str_replace_all(sample, fixed("."), fixed("-"))
    sample <- str_replace_all(sample, fixed("_"), fixed("/"))
    distances[i, "sample"] <- sample
    distances[i, "distance"] <- nrow(vcf)
    distances[i, "date"] <- as.Date(genbank_meta$date[i])
    distances[i, "day"] <- as.numeric(as.Date(genbank_meta$date[i]) - patient_day0)
    distances[i, "lineage"] <- genbank_meta$Pango.lineage[i]
    
    remove(vcf)
    
    print(paste(round((i/length(vcf_list))*100), "% of VCFs processed.", sep = ""))
    
  }
  
  write.csv(distances, 
            paste(getwd(),
                  "/", as.character(Sys.Date()), ".csv", sep = ""),
            quote = F, row.names = F)
  paste(getwd(), vcf_list, sep = "") %>%
    file.remove() %>%
    table() %>%
    as.numeric() %>%
    paste("files were cleared from disk space.", sep = " ") %>%
    print()
  
} else if (length(list.files(getwd(), "genbank_subsample_root2tip_distances")) > 0) {
  
  distance_tables <- list.files(getwd(), "genbank_subsample_root2tip_distances")
  distances <- read.csv(
    paste(getwd(), distance_tables[length(distance_tables)], sep = "/")
  )
  distances$date <- as.Date(distances$date)
  genbank_meta <- read.delim("b12_enriched_global_subsampled_metadata.tsv")
  genbank_meta <- genbank_meta[order(genbank_meta$strain),] ; rownames(genbank_meta) <- NULL
  genbank_meta$date <- as.Date(genbank_meta$date)
  
} else {
  
  print(paste("There are no VCFs or pre-processed distance tables in the current working directory, ",
        getwd(), ", so this script cannot continue plotting.", sep = ""))
  
}

if (nrow(distances)!=nrow(genbank_meta)){
  genbank_meta <- genbank_meta[match(distances$sample, genbank_meta$strain),]
}


### IDENTIFYING PANGO LINEAGES IN GenBank SUBSAMPLE ####
lineages <- as.data.frame(table(genbank_meta$Pango.lineage))
colnames(lineages) <- c("lineage", "count")
lineages <- lineages[order(lineages$count, decreasing = T),] ; rownames(lineages) <- NULL
# sum(lineages[1:5,"count"])/total_count # try to get at least 50% of the samples colored by lineage

palette <- read.delim("SupplementalFigure2_palette.txt", header = FALSE, sep = "\t")
palette$color <- c(rep(1, 5),
                   rep(2, 5),
                   rep(3, 5),
                   rep(4, 5))

# palette <- rev(c("#090c08","#757083","#8a95a5","#b9c6ae","#d6e5e3","#9fd8cb","#b4cded")) # palette from https://coolors.co/090c08-757083-8a95a5-b9c6ae-d6e5e3-9fd8cb-b4cded

lineages$color <- ""
lineages$color[1:4] <- palette[1:4,1] # brewer.pal(6, "Reds")[2:6]
lineages$color[5:nrow(lineages)] <- palette[5,1] # rep(brewer.pal(6, "Reds")[1], times = length(lineages$color[6:nrow(lineages)]))
lineages$label <- ""
lineages$label[1:4] <- as.character(lineages$lineage[1:4])
lineages$label[5:nrow(lineages)] <- rep("other", times = length(lineages$label[5:nrow(lineages)]))



### MOVING LINEAGE DATA INTO FASTA DF ####
distances$lineage <- genbank_meta$Pango.lineage
distances$color <- ""
distances$label <- ""
for (i in 1:nrow(distances)){
  
  if(is.na(distances[i, "lineage"])){
    distances[i,"color"] <- "#08266D"
    distances[i,"label"] <- "other"
    
  } else {
    distances[i,"color"] <- lineages[lineages$lineage==distances$lineage[i],"color"]
    distances[i,"label"] <- lineages[lineages$lineage==distances$lineage[i],"label"]
  }
}



### COUNTING MUTATIONS #####
variant_counts <- read.csv("patient_variant_counts.csv")



### PLOTTING ####
pdf(file = "figsupp2_global_roottotip_plot.pdf", 
    width = 9, height = 5)
plot(1:500, 1:500, 
     xlim = c(min(distances$day), max(distances$day)), ylim = c(0,90),
     xlab = "Days Post-Diagnosis", ylab = "Genetic Distance from Wuhan-1", 
     frame.plot = F, cex.axis = 0.8, cex.lab = 0.85, las = 1, pch = 20,
     cex = 3, col = "#4B7395", type = "n")
grid()

points(distances$day, distances$distance,
       pch = 20, cex = 1, col = distances$color)

text(198, 80, labels = "Bamlanivimab\nTreatment", cex = 0.8, bty = "l")
segments(198, y0 = -4, y1 = 75, col = "red", lty = 2)
segments(198, y0 = 85, y1 = 95, col = "red", lty = 2)

points(variant_counts$day, variant_counts$no_of_mutations,
       pch = 20, cex = 2, col = palette[11,1])

legend("topleft", legend = c(lineages$label[1:5], "patient"),
       col = c(lineages$color[1:5], palette[11,1]), bty="n",
       pch = 16, ncol = 2, xpd = T, xjust = 0.5, cex = 0.8)

dev.off()
