library(dplyr)

# Search for gene symbols that contain 'H3', 'H4', 'H2A', 'H2B', or start with HIST
histone_pattern <- "^(HIST|H3|H4|H2A|H2B)"

histone_peaks <- df %>% filter(grepl(histone_pattern, SYMBOL))

# Show unique histone gene symbols found
unique(histone_peaks$SYMBOL)

#Plot it

library(ggplot2)

histone_counts <- histone_peaks %>%
  count(SYMBOL, sort = TRUE)

ggplot(histone_counts, aes(x = reorder(SYMBOL, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Histone Genes with Unique WT Peaks",
       x = "Histone Gene Symbol",
       y = "Number of Peaks") +
  theme_minimal()
