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
  
orgDatabase$Ecotype <- as.factor(orgDatabase$Ecotype)

ecotypeDatabase <- orgDatabase %>% 
  dplyr::group_by(Ecotype) %>% 
  dplyr::summarize(Size_Mb = mean(Size_Mb),
                   Cell_radius = mean(Cell_radius),
                   GCpct = mean(GCpct),
                   Pmax = mean(Pmax),
                   OGT = mean(OGT),
                   n = n())
ecotypeDatabase

write.csv(ecotypeDatabase,
          "~/Desktop/Project/mse/data/db/ecotypeDatabase.csv",
          quote = FALSE,
          row.names = FALSE,
          na = "NaN"
          )
