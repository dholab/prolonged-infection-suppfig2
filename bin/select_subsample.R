#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

full_list <- read.delim(args[1])
sample_size <- as.numeric(args[2])

set.seed(14)
keepers <- sample(rownames(full_list), 
                 sample_size, 
                 replace = FALSE)

include_list <- full_list[keepers, "accession"]

write.csv(include_list, file = "include_list.csv", quote = F, row.names = F)