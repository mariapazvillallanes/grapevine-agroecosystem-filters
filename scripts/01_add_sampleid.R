library(tidyverse)

library(readr)
library(dplyr)

df <- read_delim("Hongos.csv", delim = ";", show_col_types = FALSE)
colnames(df)

df <- df %>%
  mutate(
    SampleID = paste(Management, Medium, Organ, PlantID, sep = "_")
  )

df %>%
  count(SampleID) %>%
  arrange(desc(n))

write.csv2(
  df,
  "Hongos_withSampleID.csv",
  row.names = FALSE
)

df %>%
  count(SampleID, PlantID, Management, Medium, Organ) %>%
  filter(n > 1)
