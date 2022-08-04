#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# collecting inputs/parameters
full_list <- read.delim(args[1])
sample_size <- as.numeric(args[2])
min_date <- as.Date(args[3])
max_date <- as.Date(args[4])

# filtering
full_list <- full_list[full_list$Source.database=="GenBank",]
full_list$Isolate.Collection.date <- as.Date(full_list$Isolate.Collection.date)
full_list <- full_list[!is.na(full_list$Isolate.Collection.date),]
full_list <- full_list[full_list$Isolate.Collection.date >= min_date &
                         full_list$Isolate.Collection.date <= max_date,]

# randomly sampling subset of accessions to pull from GenBank
set.seed(args[5])
keepers <- sample(rownames(full_list), 
                 sample_size, 
                 replace = FALSE)

# formatting and exporting include list
include_list <- full_list[keepers, c("Accession", "Isolate.Collection.date", "Virus.Pangolin.Classification")]
colnames(include_list) <- c("accession", "date", "pango")

write.csv(include_list, file = "include_list.csv", 
          quote = F, row.names = F)
