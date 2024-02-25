#!/bin/bash

# Load required R libraries
echo "library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)" > script.R

# Set working directory to the location of the files
echo 'setwd("./")' >> script.R

# List all the files in the directory
echo 'files <- list.files(pattern = "\\.txt$")' >> script.R

# Create an empty list to store the tables
echo 'table_list <- list()' >> script.R

# Loop through each file, read in the data and add it to the list
echo 'for (file in files) {
  # read in the data
  data <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  # check if there is a value in the second column
  if (any(data[, 2] != "")) {
    # add the file name as a column and add to the list
    table_list[[file]] <- data %>%
      mutate(file = file) %>%
      select(file, V1, V2)
  }
}' >> script.R

# Combine all tables into one data frame
echo 'combined <- bind_rows(table_list)' >> script.R

# Pivot to long format
echo 'long_format <- combined %>%
  # filter(V2 <= 78)  %>% uncomment and change value if you want to remove large values in the long range or like
  pivot_longer(-c(file, V2), names_to = "filename", values_to = "value")' >> script.R

# Generate plot
echo 'pdf("facet_wrap.pdf")
ggplot(long_format, aes(x = V2, y = value)) +
  geom_point() +
  facet_wrap(~file, ncol = 2, scales = "free_y") +
  labs(x = "V2", y = "value") +
  xlab("Read length in bp") +
  ylab("Number of reads")
dev.off()' >> script.R

# Run R script
Rscript script.R
