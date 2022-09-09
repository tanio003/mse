# get_ecotypeDatabase.R
# Script to create ecotype-averaged database from orgDatabase
library(dplyr)
library(readxl)
library(writexl)
library(naniar)
library(tidyverse)
library(plyr)
library(dplyr)
library(lubridate)
library(R.matlab)
library(rlist)
library(stringr)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(pracma)
library(zoo)
library(fuzzyjoin)
library(ggpubr)

orgDatabase <- read.csv("~/Desktop/Project/mse/data/db/orgDatabase.csv") %>%
  as_tibble() %>% dplyr::filter(Include == 1)


OGTDat <- read.csv("~/Desktop/Project/mse/data/db/OGTDat.csv")
orgDatabase <- orgDatabase %>%
  left_join(OGTDat, by = c("StrainName" = "Strain"))
  
# orgDatabase$Ecotype <- as.factor(orgDatabase$Ecotype)

ecotypeDatabase <- orgDatabase %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::summarize(Size_Mb = mean(Size_Mb),
                   Cell_radius = mean(Cell_radius),
                   GCpct = mean(GCpct),
                   Pmax = mean(Pmax),
                   OGT = mean(OGT),
                   n = n())
ecotypeDatabase

## Specific curation for HLIII_IV
# Assume Cell size, radius, and OGT for HLIII_IV are same as that of HLII as the
# experimental data do not exits 
orgDatabase_HLIII_IV <- read.csv("~/Desktop/Project/mse/data/db/orgDatabase.csv") %>%
  as_tibble() %>% 
  dplyr::filter(Ecotype == "HLIII_HLIV") %>%
  dplyr::select(Ecotype, Size_Mb, Cell_radius, GCpct, Pmax) %>%
  dplyr::mutate(OGT = NA) %>%
  dplyr::group_by(Ecotype) %>%
  dplyr::summarize(Size_Mb = mean(Size_Mb),
                   Cell_radius = mean(Cell_radius),
                   GCpct = mean(GCpct),
                   Pmax = mean(Pmax),
                   OGT = mean(OGT),
                   n = n())
# orgDatabase_HLIII_IV$Ecotype <- as.factor(orgDatabase_HLIII_IV$Ecotype)
orgDatabase_HLIII_IV$OGT = ecotypeDatabase$OGT[ecotypeDatabase$Ecotype == "HLII"]
orgDatabase_HLIII_IV$Cell_radius = ecotypeDatabase$Cell_radius[ecotypeDatabase$Ecotype == "HLII"]
orgDatabase_HLIII_IV$Pmax = ecotypeDatabase$Pmax[ecotypeDatabase$Ecotype == "HLII"]

ecotypeDatabase <- rbind(ecotypeDatabase,orgDatabase_HLIII_IV)
ecotypeDatabase <- ecotypeDatabase %>% 
  dplyr::arrange(Ecotype)

ecotypeDatabase$Ecotype <- gsub('HLIII_HLIV','HLIII_IV', ecotypeDatabase$Ecotype)
ecotypeDatabase$Ecotype <- gsub('LLII_LLIII','LLII_III', ecotypeDatabase$Ecotype)

write.csv(ecotypeDatabase,
          "~/Desktop/Project/mse/data/db/ecotypeDatabase.csv",
          quote = FALSE,
          row.names = FALSE,
          na = "NaN"
          )
