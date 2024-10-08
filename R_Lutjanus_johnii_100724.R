#### Chapter 2: Growth dynamics of Tropical Marine Snappers from the Indo-Pacific 
## Author: Clement Ng 

#### Lutjanus johnii ####

## Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

## Libraries
## Packages used for wrangling 
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(reshape2)

library(rerddap)
require(ncdf4)
library(splitstackshape)
library(utils)

# Packages used to check model assumptions
library(visdat)
library(naniar)
library(DHARMa)
library(GGally)
library(performance)

# IAPE and CV
library(FSA)

# Linear mixed models and model selection
library(arm)
library(effects)
library(lme4)
library(bbmle)
library(MuMIn)
library(AICcmodavg)
library(piecewiseSEM)

# Sliding window
library(climwin)

# Packages used for plotting
library(ggplot2)
library(ggpubr)
library(viridis)
library(lattice)
#library(rgdal)
library(sf)
library(terra)

cbPalette <- c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "#F0E442", "#CC79A7", "#999999") 

c. <- function (x) scale(x, scale = FALSE) 

#load("lutjanus_johnii.RData")

#################################################################################################################################
#################################################################################################################################

#### Data coding ####
## Lutjanus johnii from Southeast Asia 
SEA_lutjanus_johnii_txt_files <- list.files(path = ".", recursive = TRUE,
                                            pattern = "\\.txt$", 
                                            full.names = TRUE)
SEA_lutjanus_johnii_increment_data <- rbindlist(sapply(SEA_lutjanus_johnii_txt_files, fread, simplify = FALSE),
                                                use.names = TRUE, fill = TRUE, idcol = "TxtFileName")

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  subset(select = -c(V1, ypos, xpos)) %>%
  dplyr::rename("CalYear" = "IncNum v1.3c") %>%
  rename("Width" = "Thickness (mm)")

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  tidyr::separate(col="TxtFileName", sep=c(2, 7), into=c("delim1", "FishID", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim1, delim2))

SEA_lutjanus_johnii_increment_data$Width[SEA_lutjanus_johnii_increment_data$Width == 0.000] <- NA
SEA_lutjanus_johnii_increment_data <- na.omit(SEA_lutjanus_johnii_increment_data) # remove points on the edge

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  group_by(FishID) %>%
  arrange(CalYear, .by_group = TRUE) %>% 
  mutate(Age = 1:n()+1) %>% # Add one point because the first ring was not included in the growth measurements
  ungroup()

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Age == max(Age) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>% # Extract the edge as a new column in the data frame
  ungroup()

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  group_by(FishID) %>% 
  mutate(EdgeRatio = case_when(Age == (max(Age) - 1) ~ Width)) %>% 
  fill(EdgeRatio, .direction = "downup") %>% # Establish the ratio between the edge and the thickness of the previous ring to identify if rings have been missed
  mutate(EdgeRatio = Edge / EdgeRatio) %>%
  mutate(EdgeType = case_when(EdgeRatio < 0.2 ~ "Narrow", 
                              EdgeRatio >= 0.2 & EdgeRatio <= 0.7 ~ "Intermediate",
                              EdgeRatio > 0.7 ~ "Wide")) %>% 
  ungroup()

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge_Delete = case_when(Age == (max(Age)) ~ "Edge")) 

SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data[- grep("Edge", SEA_lutjanus_johnii_increment_data$Edge_Delete),]
SEA_lutjanus_johnii_increment_data <- subset(SEA_lutjanus_johnii_increment_data, select = -c(Edge_Delete))
SEA_lutjanus_johnii_increment_data <- SEA_lutjanus_johnii_increment_data %>% rename("SampleID" = "FishID") 

#### Inputting fish data into the increment data frame ####
SEA_lutjanus_johnii_fish_data <- read.csv("lutjanus_johnii_fish_data.csv")
colnames(SEA_lutjanus_johnii_fish_data)
SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  rename(Year = Capture.Year, Month = Capture.Month, Day = Capture.Day, Method = Harvest.Method, 
         TL = Total.Length..mm., FL = Fork.Length..mm., SL = Standard.Length..mm., Weight = Wet.weight..g.,
         LOtoWt = L.Otolith.Weight..g., ROtoWt = R.Otolith.Weight..g.) 

SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>% 
  mutate(LOtoWt = case_when(State.of.Left.Otolith == "Chipped" | State.of.Left.Otolith == "Broken" ~ NA, TRUE ~ LOtoWt)) %>% # remove otolith weights if they are chipped
  mutate(ROtoWt = case_when(State.of.Right.Otolith == "Chipped" | State.of.Right.Otolith == "Broken" ~ NA, TRUE ~ ROtoWt)) %>% 
  subset(select = c(1, 3:5, 12, 17:23, 29:30)) %>%
  rename(Stage = `Sex.Stage`)
SEA_lutjanus_johnii_fish_data$Stage <- as.character(SEA_lutjanus_johnii_fish_data$Stage)

unique(SEA_lutjanus_johnii_fish_data$Location)

SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  mutate(Region = case_when(Location == "Eastern Makassar" ~ "Sulawesi", 
                            
                            Location == "Melacca Straits" ~ "Malacca Straits", 
                            Location == "Strait of Malacca" ~ "Malacca Straits", 
                            Location == "Malacca" ~ "Malacca Straits", 
                            Location == "Medan" ~ "Malacca Straits", 
                            
                            Location == "Anambas" ~ "Riau Archipelago", 
                            Location == "Batam" ~ "Riau Archipelago", 
                            Location == "Kijang" ~ "Riau Archipelago", 
                            Location == "South-Eastern Sumatra" ~ "Riau Archipelago", 
                            
                            Location == "Unknown" ~ "",
                            
                            TRUE ~ Location))
unique(SEA_lutjanus_johnii_fish_data$Region)

## Reclassify the names into grouping categories
SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  mutate(Location = case_when(Location == "Melacca Straits" ~ "Malacca Straits", 
                              Location == "Strait of Malacca" ~ "Malacca Straits", 
                              Location == "Malacca" ~ "Malacca Straits", 
                              Location == "Medan" ~ "Malacca Straits", 
                              TRUE ~ Location))

SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  mutate(Location = case_when(Location == "Anambas" ~ "Riau Archipelago", 
                              Location == "Batam" ~ "Riau Archipelago", 
                              Location == "Kijang" ~ "Riau Archipelago", 
                              Location == "South-Eastern Sumatra" ~ "Riau Archipelago", 
                              TRUE ~ Location))

SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  mutate(Location = case_when(Location == "Unknown" ~ "", 
                              TRUE ~ Location))

SEA_lutjanus_johnii_fish_data <- SEA_lutjanus_johnii_fish_data %>%
  mutate(Location = case_when(Location == "" ~ "", 
                              TRUE ~ Location))
unique(SEA_lutjanus_johnii_fish_data$Location)

#### Merge to form a dataset ####
SEA_lutjanus_johnii_full_data <- merge(SEA_lutjanus_johnii_increment_data, SEA_lutjanus_johnii_fish_data, by = "SampleID")
SEA_lutjanus_johnii_full_data <- SEA_lutjanus_johnii_full_data %>% rename(Site = Location)
SEA_lutjanus_johnii_full_data$birth_month = 9 # Lutjanus johnii annuli were completed in early Sept based on Cappo et al., (2000)
SEA_lutjanus_johnii_full_data <- as.data.frame(SEA_lutjanus_johnii_full_data) 

SEA_lutjanus_johnii_full_data$TL <- as.integer(SEA_lutjanus_johnii_full_data$TL)
SEA_lutjanus_johnii_full_data$FL <- as.integer(SEA_lutjanus_johnii_full_data$FL)
SEA_lutjanus_johnii_full_data$SL <- as.integer(SEA_lutjanus_johnii_full_data$SL)

SEA_lutjanus_johnii_full_data$Year <- as.character(SEA_lutjanus_johnii_full_data$Year)
SEA_lutjanus_johnii_full_data$Month <- as.character(SEA_lutjanus_johnii_full_data$Month)
SEA_lutjanus_johnii_full_data$Day <- as.character(SEA_lutjanus_johnii_full_data$Day)
SEA_lutjanus_johnii_full_data$birth_month <- as.character(SEA_lutjanus_johnii_full_data$birth_month)

head(SEA_lutjanus_johnii_full_data, 3) ; str(SEA_lutjanus_johnii_full_data) # looks good
SEA_lutjanus_johnii_full_data %>% distinct(SampleID) # Here, we see that we have 115 unique individuals that were aged

#################################################################################################################################
#################################################################################################################################

### Input data into the model 
## Lutjanus johnii from Australia
WA_lutjanus_johnii_txt_files <- list.files(path = "../../2022 WA DPIRD/2022_DPIRD_Images/Lutjanus johnii/ImageJ", recursive = TRUE,
                                           pattern = "\\.txt$", 
                                           full.names = TRUE)
WA_lutjanus_johnii_increment_data <- rbindlist(sapply(WA_lutjanus_johnii_txt_files, fread, simplify = FALSE),
                                               use.names = TRUE, fill = TRUE, idcol = "TxtFileName") 
head(WA_lutjanus_johnii_increment_data, 3) 
## There is a mixture of files that are saved as Thickness (mm) and Thickness (µm)

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data %>% filter(!is.na(`Thickness (mm)`)) # files saved as mm
WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data %>% filter(is.na(`Thickness (mm)`)) # files saved as µm

## Data frame 1
## Processing increment width files saved as mm
WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  subset(select = -c(`Thickness (µm)`, ypos, xpos, V1)) %>%
  rename("CalYear" = "IncNum v1.3c") %>%
  rename("Width" = "Thickness (mm)")

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  tidyr::separate(col="TxtFileName", sep=c(61), into=c("delim1", "FishID"), remove=TRUE) %>% 
  subset(select = -c(delim1))

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  tidyr::separate(col="FishID", sep= "_", into=c("FishID", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim2))

WA_lutjanus_johnii_increment_data_1$Width[WA_lutjanus_johnii_increment_data_1$Width == 0.000] <- NA
WA_lutjanus_johnii_increment_data_1 <- na.omit(WA_lutjanus_johnii_increment_data_1) # remove points on the edge

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  group_by(FishID) %>%
  arrange(CalYear, .by_group = TRUE) %>% 
  mutate(Age = 1:n()+1) %>% # Add one point because the first ring was not included in the growth measurements
  ungroup()

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Age == max(Age) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>%
  ungroup()

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  group_by(FishID) %>% 
  mutate(EdgeRatio = case_when(Age == (max(Age) - 1) ~ Width)) %>% 
  fill(EdgeRatio, .direction = "downup") %>%
  mutate(EdgeRatio = Edge / EdgeRatio) %>% # Establish the ratio between the edge and the thickness of the previous ring to identify if rings have been missed
  mutate(EdgeType = case_when(EdgeRatio < 0.2 ~ "Narrow",
                              EdgeRatio >= 0.2 & EdgeRatio <= 0.7 ~ "Intermediate",
                              EdgeRatio > 0.7 ~ "Wide")) %>%
  ungroup()

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% 
  group_by(FishID) %>% 
  mutate(Edge_Delete = case_when(Age == (max(Age)) ~ "Edge"))

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1[- grep("Edge", WA_lutjanus_johnii_increment_data_1$Edge_Delete),]
WA_lutjanus_johnii_increment_data_1 <- subset(WA_lutjanus_johnii_increment_data_1, select = -c(Edge_Delete))
WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% rename("SampleID" = "FishID") 
WA_lutjanus_johnii_increment_data_1 <- as.data.frame(WA_lutjanus_johnii_increment_data_1)

## Data frame 2
## Processing increment width files saved as µm
WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  subset(select = -c(`Thickness (mm)`, ypos, xpos, V1)) %>%
  rename("CalYear" = "IncNum v1.3c") %>%
  rename("Width" = "Thickness (µm)")

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  tidyr::separate(col="TxtFileName", sep=c(61), into=c("delim1", "FishID"), remove=TRUE) %>% 
  subset(select = -c(delim1))

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  tidyr::separate(col="FishID", sep= "_", into=c("Species", "FishID", "PhotoNumber", "FileType"), remove=TRUE) %>% 
  subset(select = -c(Species, PhotoNumber, FileType))

WA_lutjanus_johnii_increment_data_2$Width[WA_lutjanus_johnii_increment_data_2$Width == 0.000] <- NA
WA_lutjanus_johnii_increment_data_2 <- na.omit(WA_lutjanus_johnii_increment_data_2) # remove points on the edge

WA_lutjanus_johnii_increment_data_2$Width <- WA_lutjanus_johnii_increment_data_2$Width / 1000 # Converting µm into mm

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  group_by(FishID) %>%
  arrange(CalYear, .by_group = TRUE) %>% 
  mutate(Age = 1:n()+1) %>% # Add one point because the first ring was not included in the growth measurements
  ungroup()

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Age == max(Age) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>% 
  ungroup()

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  group_by(FishID) %>% 
  mutate(EdgeRatio = case_when(Age == (max(Age) - 1) ~ Width)) %>% 
  fill(EdgeRatio, .direction = "downup") %>%
  mutate(EdgeRatio = Edge / EdgeRatio) %>% # Establish the ratio between the edge and the thickness of the previous ring to identify if rings have been missed
  mutate(EdgeType = case_when(EdgeRatio < 0.2 ~ "Narrow",
                              EdgeRatio >= 0.2 & EdgeRatio <= 0.7 ~ "Intermediate",
                              EdgeRatio > 0.7 ~ "Wide")) %>%
  ungroup()

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% 
  group_by(FishID) %>% 
  mutate(Edge_Delete = case_when(Age == (max(Age)) ~ "Edge"))

WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2[- grep("Edge", WA_lutjanus_johnii_increment_data_2$Edge_Delete),]
WA_lutjanus_johnii_increment_data_2 <- subset(WA_lutjanus_johnii_increment_data_2, select = -c(Edge_Delete))
WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% rename("SampleID" = "FishID") 
WA_lutjanus_johnii_increment_data_2 <- as.data.frame(WA_lutjanus_johnii_increment_data_2)

#### Inputting fish data into the increment data frame ####
## There are two data sheets containing the data -- sheet 1 contains sectioned otolith metadata and sheet 4 includes the unsectioned otolith metadata
## Sectioned otolith metadata
WA_lutjanus_johnii_fish_data_1 <- read_xlsx("../../2022 WA DPIRD/DPIRD_Snapper_Biological_Data/DPIRD_Lutjanus_johnii.xlsx", 1) 
WA_lutjanus_johnii_fish_data_1 <- WA_lutjanus_johnii_fish_data_1 %>%
  rename(SampleID = `Slide Number`, Capture_Year = Year, Capture_Month = Month, Weight = TotalWt) %>%
  subset(select = -c(1, 3)) 

head(WA_lutjanus_johnii_fish_data_1, 3) ; str(WA_lutjanus_johnii_fish_data_1)
WA_lutjanus_johnii_fish_data_1$SampleID <- as.character(WA_lutjanus_johnii_fish_data_1$SampleID)
WA_lutjanus_johnii_fish_data_1$Capture_Year <- as.character(WA_lutjanus_johnii_fish_data_1$Capture_Year)
WA_lutjanus_johnii_fish_data_1$Capture_Month <- as.character(WA_lutjanus_johnii_fish_data_1$Capture_Month)
WA_lutjanus_johnii_fish_data_1$Stage <- as.character(WA_lutjanus_johnii_fish_data_1$Stage)
WA_lutjanus_johnii_fish_data_1 <-  as.data.frame(WA_lutjanus_johnii_fish_data_1)

colnames(WA_lutjanus_johnii_fish_data_1)
setcolorder(WA_lutjanus_johnii_fish_data_1, c("SampleID", "Region", "Location", 
                                              "Capture_Year", "Capture_Month", 
                                              "TL",  "FL", "SL", "Weight", "Sex","GonadWt"))

## Unsectioned otolith metadata
WA_lutjanus_johnii_fish_data_2 <- read_xlsx("../../2022 WA DPIRD/2022_Wakefield_Fish Measurements/2022_Lutjanus_johnii_Chronology.xlsx", 4)
WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>%
  subset(select = -c(1, 3, 4, 6, 8, 9, 11, 12, 20:29)) %>%
  rename(Stage = `Macro Stage`)
WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>%
  rename(SampleID = `Fish code`, Region = Region_2, Location = Site_3, TL = `TL (mm)`, FL = `FL (mm)`, SL = `SL (mm)`, Weight = `Whole Weight (g)`)

## Format columns
WA_lutjanus_johnii_fish_data_2 <- as.data.frame(WA_lutjanus_johnii_fish_data_2)
WA_lutjanus_johnii_fish_data_2$TL <- as.numeric(WA_lutjanus_johnii_fish_data_2$TL)
WA_lutjanus_johnii_fish_data_2$FL <- as.numeric(WA_lutjanus_johnii_fish_data_2$FL)
WA_lutjanus_johnii_fish_data_2$SL <- as.numeric(WA_lutjanus_johnii_fish_data_2$SL)
WA_lutjanus_johnii_fish_data_2$Weight <- as.numeric(WA_lutjanus_johnii_fish_data_2$Weight)

WA_lutjanus_johnii_fish_data_2$Capture_Year <- format(as.Date(WA_lutjanus_johnii_fish_data_2$Date), "%Y")
WA_lutjanus_johnii_fish_data_2$Capture_Month <- format(as.Date(WA_lutjanus_johnii_fish_data_2$Date), "%m")
WA_lutjanus_johnii_fish_data_2$Capture_Day <- format(as.Date(WA_lutjanus_johnii_fish_data_2$Date), "%d")

WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>% subset(select = -c(Date))
WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>% mutate(SampleID = gsub(" ", "", SampleID))

colnames(WA_lutjanus_johnii_fish_data_2)
setcolorder(WA_lutjanus_johnii_fish_data_2, c("SampleID","Region","Location","Method", 
                                              "Capture_Year","Capture_Month","Capture_Day", 
                                              "TL","FL","SL","Weight","Sex","Stage"))

## Amend format of the increment data to facilitate the merging of data frames
WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% mutate(SampleID = gsub("-", "", SampleID))
WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% mutate(SampleID = gsub("-", "", SampleID))
WA_lutjanus_johnii_fish_data_1 <- WA_lutjanus_johnii_fish_data_1 %>% mutate(SampleID = gsub("-", "", SampleID))
WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>% mutate(SampleID = gsub("-", "", SampleID))

WA_lutjanus_johnii_increment_data_1 <- WA_lutjanus_johnii_increment_data_1 %>% mutate(SampleID = gsub("LJ", "", SampleID))
WA_lutjanus_johnii_increment_data_2 <- WA_lutjanus_johnii_increment_data_2 %>% mutate(SampleID = gsub("LJ", "", SampleID))
WA_lutjanus_johnii_fish_data_1 <- WA_lutjanus_johnii_fish_data_1 %>% mutate(SampleID = gsub("LJ", "", SampleID))
WA_lutjanus_johnii_fish_data_2 <- WA_lutjanus_johnii_fish_data_2 %>% mutate(SampleID = gsub("LJ", "", SampleID))

#### Merge to form a dataset ####
WA_lutjanus_johnii_full_data_1 <- merge(WA_lutjanus_johnii_increment_data_1, WA_lutjanus_johnii_fish_data_2, by = "SampleID")
WA_lutjanus_johnii_full_data_2 <- merge(WA_lutjanus_johnii_increment_data_2, WA_lutjanus_johnii_fish_data_1, by = "SampleID")
WA_lutjanus_johnii_full_data_3 <- merge(WA_lutjanus_johnii_increment_data_2, WA_lutjanus_johnii_fish_data_2, by = "SampleID")

WA_lutjanus_johnii_full_data_2$TL <- as.numeric(WA_lutjanus_johnii_full_data_2$TL)
WA_lutjanus_johnii_full_data_2$Weight <- as.numeric(WA_lutjanus_johnii_full_data_2$Weight)

WA_lutjanus_johnii_full_data <- dplyr::bind_rows(WA_lutjanus_johnii_full_data_1, WA_lutjanus_johnii_full_data_2, WA_lutjanus_johnii_full_data_3)
WA_lutjanus_johnii_full_data <- WA_lutjanus_johnii_full_data %>% rename(Site = Location)

WA_lutjanus_johnii_full_data$birth_month = 9
WA_lutjanus_johnii_full_data$birth_month <- as.character(WA_lutjanus_johnii_full_data$birth_month)
WA_lutjanus_johnii_full_data$TL <- as.integer(WA_lutjanus_johnii_full_data$TL)
WA_lutjanus_johnii_full_data$FL <- as.integer(WA_lutjanus_johnii_full_data$FL)
WA_lutjanus_johnii_full_data$SL <- as.integer(WA_lutjanus_johnii_full_data$SL)
WA_lutjanus_johnii_full_data$Weight <- as.integer(WA_lutjanus_johnii_full_data$Weight)

WA_lutjanus_johnii_full_data <- WA_lutjanus_johnii_full_data %>% 
  rename(Year = Capture_Year, Month = Capture_Month, Day = Capture_Day) 
WA_lutjanus_johnii_full_data <- as.data.frame(WA_lutjanus_johnii_full_data)

head(WA_lutjanus_johnii_full_data, 3) ; str(WA_lutjanus_johnii_full_data) # looks good
WA_lutjanus_johnii_full_data %>% distinct(SampleID) # Here, we see that we have 451 unique individuals that were aged

WA_lutjanus_johnii_full_data <- WA_lutjanus_johnii_full_data %>% arrange(SampleID, CalYear)

## Additional check for missing metadata 
checkMissing <- WA_lutjanus_johnii_full_data %>%
  filter(!is.na(Width)) %>% filter(is.na(Region)) # Filter Regions without any information but has increment Width data
checkMissing

checkMissing %>%
  filter(!is.na(Width)) %>% filter(is.na(Region)) %>% # Filter Regions without any information but has increment Width data
  distinct(SampleID) 
#data.table::fwrite(checkMissing, file = "missingMetadata.csv") 

#################################################################################################################################
#################################################################################################################################

#### Combine the Southeast Asia and Western Australia data sets ####
lutjanus_johnii_full_data <- merge(SEA_lutjanus_johnii_full_data, WA_lutjanus_johnii_full_data, all=TRUE)
lutjanus_johnii_full_data <- as.data.frame(lutjanus_johnii_full_data)
head(lutjanus_johnii_full_data, 3) ; str(lutjanus_johnii_full_data)

lutjanus_johnii_full_data %>% distinct(SampleID) 

lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>% 
  group_by(SampleID) %>% 
  mutate(AAC = max(Age)) # fully formed rings

lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>%
  mutate(Cohort = CalYear-Age)

lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>% rename(FishYear = CalYear)

lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>%
  mutate(MaturityStage = case_when(Age < 5 ~ "Immature",
                                   Age >= 5 ~ "Mature"))
## Personal communication -- Corey mentioned in our meeting that the age at maturity of Lutjanus johnii was 4 years old. 

unique(lutjanus_johnii_full_data$Region)
lutjanus_johnii_full_data$Region[lutjanus_johnii_full_data$Region == "North Kimberley"] <- "Kimberley"
lutjanus_johnii_full_data$Region[lutjanus_johnii_full_data$Region == ""] <- NA
lutjanus_johnii_full_data %>% group_by(Region) %>% summarise(n=n_distinct(SampleID))

## Discussions with DPIRD suggested that I should combine the Regions due to the lack of genetic differences between populations for Lutjanus sebae, L. malabaricus, and P. multidens.
## However, Taillebois et al., (2021) found small-scale spatial genetic differences in L. johnii across the northern coastline of Australia. 
lutjanus_johnii_full_data$Site <- lutjanus_johnii_full_data$Region
unique(lutjanus_johnii_full_data$Site)

## Remove additional Sites from the sampling due to inadequate sample sizes for constructing a growth chronology
lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>% 
  filter(Site == "Kimberley" | Site == "Malacca Straits" | Site == "Riau Archipelago")

lutjanus_johnii_full_data <- lutjanus_johnii_full_data %>%
  mutate(Region = case_when(Site == "Kimberley" ~ "NW Australia",
                            Site == "Malacca Straits" | Site == "Riau Archipelago" ~ "Southeast Asia"))
unique(lutjanus_johnii_full_data$Region)

unique(lutjanus_johnii_full_data$Sex)
lutjanus_johnii_full_data$Sex[lutjanus_johnii_full_data$Sex == ""] <- NA
lutjanus_johnii_full_data$Sex[lutjanus_johnii_full_data$Sex == "m"] <- "M"
lutjanus_johnii_full_data$Sex[lutjanus_johnii_full_data$Sex == "f"] <- "F"
lutjanus_johnii_full_data$Sex[lutjanus_johnii_full_data$Sex == "U"] <- NA

### Create grouping variables
lutjanus_johnii_full_data <- as.data.frame(lutjanus_johnii_full_data)
lutjanus_johnii_full_data$SampleID <- as.factor(lutjanus_johnii_full_data$SampleID)
lutjanus_johnii_full_data$Region <- as.factor(lutjanus_johnii_full_data$Region)
lutjanus_johnii_full_data$Site <- as.factor(lutjanus_johnii_full_data$Site)
lutjanus_johnii_full_data$Cohort <- as.factor(lutjanus_johnii_full_data$Cohort)

lutjanus_johnii_full_data$region_year <- as.factor(paste(lutjanus_johnii_full_data$Region, lutjanus_johnii_full_data$FishYear, sep = "_")) 
lutjanus_johnii_full_data$site_year <- as.factor(paste(lutjanus_johnii_full_data$Site, lutjanus_johnii_full_data$FishYear, sep = "_")) 
lutjanus_johnii_full_data$cohort_year <- as.factor(paste(lutjanus_johnii_full_data$Cohort, lutjanus_johnii_full_data$FishYear, sep = "_")) 
lutjanus_johnii_full_data$region_cohort <- as.factor(paste(lutjanus_johnii_full_data$Region, lutjanus_johnii_full_data$Cohort, sep = "_")) 
lutjanus_johnii_full_data$site_cohort <- as.factor(paste(lutjanus_johnii_full_data$Site, lutjanus_johnii_full_data$Cohort, sep = "_")) 

lutjanus_johnii_full_data$sAge <- c.(log(lutjanus_johnii_full_data$Age))
lutjanus_johnii_full_data$sAAC <- c.(log(lutjanus_johnii_full_data$AAC))

colnames(lutjanus_johnii_full_data)
setcolorder(lutjanus_johnii_full_data, c("SampleID","FishYear", 
                                         "Year","Month","Day","birth_month", 
                                         "Method","Region","region_year","Site","site_year", 
                                         "Cohort","cohort_year","site_cohort","region_cohort",
                                         "TL","FL","SL","Weight","Sex","Stage","GonadWt",
                                         "LOtoWt","ROtoWt",
                                         "Width","Age","sAge","AAC","sAAC",
                                         "Edge","EdgeRatio","EdgeType"))

#################################################################################################################################
#################################################################################################################################

## Intrinsic variables
## Note that this data set is based on increments, so we need to split it based on distinct individuals based on the function, summarise(n_distinct(SampleID))
## Age at capture table and distribution
lutjanus_johnii_full_data %>%
  group_by(AAC) %>% summarise(n_distinct(SampleID)) %>% # summarise based on long format
  rename(n = `n_distinct(SampleID)`) %>%
  #spread(AAC, n) %>% # change to long format, if required
  print(n=Inf)

lutjanus_johnii_full_data %>% 
  group_by(AAC) %>% summarise(n=n_distinct(SampleID)) %>%
  ggplot(., aes(x=AAC, y=n))+ # histogram to highlight age distribution
  geom_bar(stat='identity',color='grey5',fill='grey45')+
  labs(y= 'Count')+ theme_bw()

## Region
lutjanus_johnii_full_data %>% 
  group_by(Region) %>% summarise(n=n_distinct(SampleID)) %>%
  print(n=Inf)

lutjanus_johnii_full_data %>% 
  group_by(AAC, Region) %>% summarise(n=n_distinct(SampleID)) %>%
  ggplot(., aes(x=AAC, y=n))+ # histogram to highlight AAC distribution based on Region
  geom_bar(stat='identity',color='grey5',fill='grey45')+
  scale_y_continuous(limits=c(0, 120))+
  scale_x_continuous(limits=c(0, 30))+
  facet_wrap(.~Region, ncol = 2)+
  labs(y= 'Count', x= 'Age-at-Capture')+ theme_bw()

## Total length
## Replace with Standard and Fork lengths respectively. 
lutjanus_johnii_full_data %>% 
  group_by(TL) %>% summarise(n=n_distinct(SampleID)) %>%
  print(n=Inf)

lutjanus_johnii_full_data %>% 
  group_by(TL, Region) %>% summarise(n=n_distinct(SampleID)) %>% drop_na(TL) %>%
  ggplot(., aes(y=as.numeric(TL), x=Region, group=Region))+ # histogram to highlight age distribution
  geom_boxplot(color='grey5',fill='grey45')+ geom_jitter()+
  #facet_wrap(.~Region)+
  labs(y= 'Total Length (mm)')+ theme_bw()

#################################################################################################################################
#################################################################################################################################

#### Check missing values ####
colSums(is.na(lutjanus_johnii_full_data)) # note that this is repeated for all increments, not unique values
# There are a large number of missing values for the left and right otolith weights, as they were not recorded in Australia. 
# There are also missing methods of capture, FL and SL measurements. This indicates that we should ideally use TL measurements. 

## Create an additional table
checkTable <- lutjanus_johnii_full_data %>% 
  filter(Age=='2') 
colSums(is.na(checkTable))
# Check for missing values based on two year old fish
# This works because 2 is the minimum age of fish that are being sampled in my study
checkTable <- as.data.frame(checkTable)

## Visualise missing values
colnames(checkTable)
naniar::gg_miss_var(checkTable, show_pct = T) 
visdat::vis_miss(checkTable) 
# Based on the missing variables function, there are almost no values for gonad weight. 
# There is also a large number of missing values for Right and Left otolith weights and the margin status of the otolith. 
# A proportion of the data entries are missing Day, Method and the Standard Length. 

## Assess the correlation between fish length and otolith weight
checkTable <- checkTable %>% drop_na(TL, ROtoWt)
cor.test(checkTable$TL, checkTable$ROtoWt, method = "pearson") 

## Assess the relationship between left and right otoliths
checkTable <- checkTable %>% drop_na(LOtoWt, ROtoWt)
t.test(checkTable$LOtoWt, checkTable$ROtoWt) 

### Check for outliers in the response variable ####
## Examine the relationship between increment width and fish age for each individual
xyplot(Width ~ Age, group=SampleID, lutjanus_johnii_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus johnii)", ylab = "Width (mm)") # Overall pattern

xyplot(Width ~ Age | Region, group=SampleID, lutjanus_johnii_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus johnii)", ylab = "Width (mm)")

xyplot(Width ~ Age | Site, group=SampleID, lutjanus_johnii_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus johnii)", ylab = "Width (mm)") # Regional pattern
# There may potentially be three outliers in the data frame, present in Riau Archipelago, and two data points in the Malacca Straits. 
# There may also be an outlier in the Kimberley data frame, a 7 year old individual with a large increment width.
# Worth going back to the data frame to check the data and see if a ring was missed. 

xyplot(log(Width) ~ log(Age), group=SampleID, lutjanus_johnii_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus johnii)", ylab = "Width (mm)") # Log transformed

### Check for dependency 
lutjanus_johnii_full_data %>% 
  ggplot()+ 
  geom_point(aes(y=Width, x=Region, group=SampleID), position=position_dodge(0.8))+
  theme(text = element_text(size = 12), 
        legend.position = "none", 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 0.7))

## Create additional grouping variables 
lutjanus_johnii_full_data$region_stage_year <- as.factor(paste(lutjanus_johnii_full_data$Region, lutjanus_johnii_full_data$MaturityStage, lutjanus_johnii_full_data$FishYear, sep = "_")) 
lutjanus_johnii_full_data$region_stage_cohort <- as.factor(paste(lutjanus_johnii_full_data$Region, lutjanus_johnii_full_data$MaturityStage, lutjanus_johnii_full_data$Cohort, sep = "_")) 

#################################################################################################################################
#################################################################################################################################

lutjanus_johnii_data <- lutjanus_johnii_full_data

unique(lutjanus_johnii_data$Method)
lutjanus_johnii_data <- lutjanus_johnii_data %>%
  mutate(Method = case_when(Method == "" ~ NA, 
                            Method == "Traps" ~ "Trap",
                            TRUE ~ Method))

lutjanus_johnii_data <- lutjanus_johnii_data %>% subset(select = -c(GonadWt))

sum(is.na(lutjanus_johnii_data$Site))
lutjanus_johnii_data <- lutjanus_johnii_data %>% drop_na(Site) 

## For reading 1 
#AAC_Measurements <- lutjanus_johnii_data %>% distinct(SampleID, Region, AAC) 
#write.csv(AAC_Measurements, 'LJ_AAC_Measurements_4.csv')

### Split the data frame into the different geographic regions 
unique(lutjanus_johnii_data$Region)
lutjanus_johnii_data <- lutjanus_johnii_data %>%
  mutate(Region = case_when(Region == "NW Australia" ~ "North Western Australia", 
                            TRUE ~ Region))
lutjanus_johnii_data

#################################################################################################################################
#################################################################################################################################

## Following our discussion on 31 May 2024, we thought it would be best to split my models into region specific ones. 
## This is because we have unequal sample sizes for Kimberley versus Malacca Straits versus Riau Archipelago.
## The underlying driver of growth would then be explained by Kimberley, which has a much longer time series, compared to the Southeast Asian samples. 

Ljohnii_Australia <- lutjanus_johnii_data %>% filter(Region == "North Western Australia")
Ljohnii_Asia <- lutjanus_johnii_data %>% filter(Region == "Southeast Asia")

### Construct the basic intrinsic linear mixed models 
## Australian samples
M1a <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID), 
            Ljohnii_Australia, REML= T)

M1b <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID), 
            Ljohnii_Australia, REML= T)

M1c <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|FishYear), 
            Ljohnii_Australia, REML= T)

M1d <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Australia, REML= T)

M1e <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|Cohort), 
            Ljohnii_Australia, REML= T)

M1f <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Ljohnii_Australia, REML= T)

bbmle::AICctab(M1a,M1b,M1c,M1d,M1e,M1f, base=T,logLik=T,weights=T) 

### Extract the conditional and marginal R-squared values using the piecewiseSEM package
rsquared(M1a)
rsquared(M1b)
rsquared(M1c)
rsquared(M1d)
rsquared(M1e)
rsquared(M1f)

### Refitting with Maximum Likelihood Estimation (ML) 
## Fit the most complex model and perform automated model selection on it to identify the best model structure
M1 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Ljohnii_Australia, REML= F)

options(na.action=na.fail)
M1_dredge<-dredge(M1, trace=2) 

subset(M1_dredge)

## Optimal intrinsic mixed model structure comprises of
M1 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Ljohnii_Australia, REML= T)

plot(M1) # evenly distributed
qqnorm(residuals(M1)) ; qqline(residuals(M1)) # slight skewness in the lower values and the high values

summary(M1) 
anova(M1)
rsquared(M1)

M1_conf <- sim(M1, 5000)
M1_conf <- apply(M1_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M1_conf

p <- predictorEffect("sAge", (M1)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))
p <- predictorEffect("sAAC", (M1)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))

#### Year model for Australian samples ###
M1_year <- ranef(M1)$FishYear[,1]
M1_year_se <-sqrt (attr(ranef(M1, postVar=TRUE) [["FishYear"]], "postVar")[1,1,])

modelyear1 <- data.frame(y=M1_year)
modelyear1$upper <- (modelyear1$y + M1_year_se)
modelyear1$lower <- (modelyear1$y - M1_year_se)
modelyear1$year <- rownames(ranef(M1)$FishYear)
modelyear1$year <- as.integer(modelyear1$year)
sample_depth_1 <- Ljohnii_Australia %>% 
  group_by(FishYear) %>% 
  summarise(N=n())
colnames(sample_depth_1)[1] <- "year"
modelyear1 <- modelyear1 %>% 
  left_join(sample_depth_1) 
modelyear1

#### Plot relative annual otolith growth of Australian samples over time ###
plot1 <- modelyear1 %>% 
  ggplot(aes(y=y, x=year))+
  geom_point(size=1)+ 
  geom_line(linewidth=0.5)+
  geom_ribbon(aes(y=y, ymin=lower, ymax=upper), alpha=0.4)+
  scale_x_continuous(limits=c(1975,2023), breaks=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020), labels=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020))+
  labs(y="Relative annual otolith growth", x="Year")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot1

#### Predict how Lutjanus johnii growth from Australia varies with age ####
Ljohnii_Australia %>% 
  dplyr::slice(which.max(Age)) %>% 
  pull(Age)

## Age plot
Ljohnii_Australia_Age <- as.data.frame (Effect (c('sAge'), M1, xlevels = list(sAge=seq(-0.974, 3.91 ,by=0.01)) )) 
mean_Australia_age <- mean(Ljohnii_Australia$Age) ; mean_Australia_age
sd_Australia_age <- sd(Ljohnii_Australia$Age) ; sd_Australia_age
Ljohnii_Australia_Age$age <- (Ljohnii_Australia_Age$sAge * sd_Australia_age + mean_Australia_age)

Ljohnii_Australia_Age$transfit<-exp(Ljohnii_Australia_Age$fit)
Ljohnii_Australia_Age$transupper<-exp(Ljohnii_Australia_Age$upper)
Ljohnii_Australia_Age$translower<-exp(Ljohnii_Australia_Age$lower)

plot2 <- Ljohnii_Australia_Age %>% 
  filter(age <= 20.0) %>%
  ggplot(aes(y=transfit, x=age)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age", y = "Predicted annual otolith growth (mm)", fill = "Site", color = "Site")+ 
  scale_x_continuous(limits=c(0,20), breaks = c(0,5,10,15,20), labels = c('0','5','10','15','20'))+
  scale_y_continuous(limits=c(0.0,0.8), breaks = c(0.000, 0.200, 0.400, 0.600, 0.800), labels = c('0.0','0.2','0.4','0.6','0.8'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.position.inside = c(1,1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot2

## AAC plot
Ljohnii_Australia_AAC <- as.data.frame (Effect (c('sAAC'), M1, xlevels = list(sAAC=seq(-1.458, 2.17 ,by=0.01)) )) 
mean_Australia_AAC <- mean(Ljohnii_Australia$AAC) ; mean_Australia_AAC
sd_Australia_AAC <- sd(Ljohnii_Australia$AAC) ; sd_Australia_AAC
Ljohnii_Australia_AAC$AAC <- (Ljohnii_Australia_AAC$sAAC * sd_Australia_AAC + mean_Australia_AAC)

Ljohnii_Australia_AAC$transfit<-exp(Ljohnii_Australia_AAC$fit)
Ljohnii_Australia_AAC$transupper<-exp(Ljohnii_Australia_AAC$upper)
Ljohnii_Australia_AAC$translower<-exp(Ljohnii_Australia_AAC$lower)

plot3 <- Ljohnii_Australia_AAC %>% 
  filter(AAC <= 20.0) %>%
  ggplot(aes(y=transfit, x=AAC)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age-at-Capture", y = "Predicted annual otolith growth (mm)")+ 
  scale_x_continuous(limits=c(0,20), breaks = c(0,5,10,15,20), labels = c('0','5','10','15','20'))+
  scale_y_continuous(limits=c(0.0,0.6), breaks = c(0.000, 0.200, 0.400, 0.600), labels = c('0.0','0.2','0.4','0.6'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.position = "none",
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14)
  )
plot3

ggpubr::ggarrange(plot2, plot3, labels = c("A", "B"))

#################################################################################################################################
#################################################################################################################################

### Construct the basic intrinsic linear mixed models 
## Southeast Asian samples
M2a <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID), 
            Ljohnii_Asia, REML= T)

M2b <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID), 
            Ljohnii_Asia, REML= T)

M2c <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|FishYear), 
            Ljohnii_Asia, REML= T)

M2d <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Asia, REML= T)

M2e <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|Cohort), 
            Ljohnii_Asia, REML= T)

M2f <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Ljohnii_Asia, REML= T)

bbmle::AICctab(M2a,M2b,M2c,M2d,M2e,M2f, base=T,logLik=T,weights=T) 

### Extract the conditional and marginal R-squared values using the piecewiseSEM package
rsquared(M2a)
rsquared(M2b)
rsquared(M2c)
rsquared(M2d)
rsquared(M2e)
rsquared(M2f)

### Refitting with Maximum Likelihood Estimation (ML) 
## Fit the most complex model and perform automated model selection on it to identify the best model structure
M2 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Ljohnii_Asia, REML= F)

options(na.action=na.fail)
M2_dredge<-dredge(M2, trace=2) 

subset(M2_dredge)

## Optimal intrinsic mixed model structure comprises of
M2 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Ljohnii_Asia, REML= T)

plot(M2) # evenly distributed
qqnorm(residuals(M2)) ; qqline(residuals(M2)) # slight skewness in the lower values and the high values

summary(M2) 
anova(M2)
rsquared(M2)

M2_conf <- sim(M2, 5000)
M2_conf <- apply(M2_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M2_conf

p <- predictorEffect("sAge", (M2)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))
p <- predictorEffect("sAAC", (M2)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))

#### Year model for Southeast Asian samples ###
M2_year <- ranef(M2)$FishYear[,1]
M2_year_se <-sqrt (attr(ranef(M2, postVar=TRUE) [["FishYear"]], "postVar")[1,1,])

modelyear2 <- data.frame(y=M2_year)
modelyear2$upper <- (modelyear2$y + M2_year_se)
modelyear2$lower <- (modelyear2$y - M2_year_se)
modelyear2$year <- rownames(ranef(M2)$FishYear)
modelyear2$year <- as.integer(modelyear2$year)
sample_depth_2 <- Ljohnii_Asia %>% 
  group_by(FishYear) %>% 
  summarise(N=n())
colnames(sample_depth_2)[1] <- "year"
modelyear2 <- modelyear2 %>% 
  left_join(sample_depth_2) 
modelyear2

#### Plot relative annual otolith growth of Southeast Asian samples over time ###
plot4 <- modelyear2 %>% 
  ggplot(aes(y=y, x=year))+
  geom_point(size=1)+ 
  geom_line(linewidth=0.5)+
  geom_ribbon(aes(y=y, ymin=lower, ymax=upper), alpha=0.4)+
  scale_x_continuous(limits=c(1975,2023), breaks=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020), labels=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020))+
  scale_y_continuous(limits=c(-0.20,0.20))+
  labs(y="Relative annual otolith growth", x="Year")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot4

ggpubr::ggarrange(plot1, plot4, labels = c("A", "B"), nrow=2)

#### Predict how Lutjanus johnii growth from Asia varies with age ####
Ljohnii_Asia %>% 
  dplyr::slice(which.max(Age)) %>% 
  pull(Age)

## Age plot
Ljohnii_Asia_Age <- as.data.frame (Effect (c('sAge'), M2, xlevels = list(sAge=seq(-1.12, 4.17 ,by=0.01)) )) 
mean_Asia_age <- mean(Ljohnii_Asia$Age) ; mean_Asia_age
sd_Asia_age <- sd(Ljohnii_Asia$Age) ; sd_Asia_age
Ljohnii_Asia_Age$age <- (Ljohnii_Asia_Age$sAge * sd_Asia_age + mean_Asia_age)

Ljohnii_Asia_Age$transfit<-exp(Ljohnii_Asia_Age$fit)
Ljohnii_Asia_Age$transupper<-exp(Ljohnii_Asia_Age$upper)
Ljohnii_Asia_Age$translower<-exp(Ljohnii_Asia_Age$lower)

plot5 <- Ljohnii_Asia_Age %>% 
  filter(age <= 27.0) %>%
  ggplot(aes(y=transfit, x=age)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age", y = "Predicted annual otolith growth (mm)", fill = "Site", color = "Site")+ 
  scale_x_continuous(limits=c(0,27), breaks = c(0,5,10,15,20,25), labels = c('0','5','10','15','20','25'))+
  scale_y_continuous(limits=c(0.0,0.8), breaks = c(0.000, 0.200, 0.400, 0.600, 0.800), labels = c('0.0','0.2','0.4','0.6','0.8'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.position.inside = c(1,1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot5

## AAC plot
Ljohnii_Asia_AAC <- as.data.frame (Effect (c('sAAC'), M2, xlevels = list(sAAC=seq(-1.825, 2.49 ,by=0.01)) )) 
mean_Asia_AAC <- mean(Ljohnii_Asia$AAC) ; mean_Asia_AAC
sd_Asia_AAC <- sd(Ljohnii_Asia$AAC) ; sd_Asia_AAC
Ljohnii_Asia_AAC$AAC <- (Ljohnii_Asia_AAC$sAAC * sd_Asia_AAC + mean_Asia_AAC)

Ljohnii_Asia_AAC$transfit<-exp(Ljohnii_Asia_AAC$fit)
Ljohnii_Asia_AAC$transupper<-exp(Ljohnii_Asia_AAC$upper)
Ljohnii_Asia_AAC$translower<-exp(Ljohnii_Asia_AAC$lower)

plot6 <- Ljohnii_Asia_AAC %>% 
  filter(AAC <= 27.0) %>%
  ggplot(aes(y=transfit, x=AAC)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age-at-Capture", y = "Predicted annual otolith growth (mm)")+ 
  scale_x_continuous(limits=c(0,27), breaks = c(0,5,10,15,20,25), labels = c('0','5','10','15','20','25'))+
  scale_y_continuous(limits=c(0.0,0.4), breaks = c(0.000, 0.100, 0.200, 0.300, 0.400), labels = c('0.0','0.1','0.2','0.3','0.4'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.position = "none",
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14)
  )
plot6

ggpubr::ggarrange(plot5, plot6, labels = c("A", "B"))

plot7 <- ggplot() + 
  geom_line(aes(y = transfit, x = age, colour="Australia"), linewidth=0.5, data = Ljohnii_Australia_Age) + 
  geom_ribbon(aes(y = transfit, x = age, ymin = translower, ymax = transupper, fill="Australia"), linewidth=0.4, alpha=0.3, colour=0, data = Ljohnii_Australia_Age) +
  geom_line(aes(y = transfit, x = age, colour="Southeast Asia"), linewidth=0.5,data = Ljohnii_Asia_Age) + 
  geom_ribbon(aes(y = transfit, x = age, ymin = translower, ymax = transupper, fill="Southeast Asia"), linewidth=0.4, alpha=0.3, colour=0, data = Ljohnii_Asia_Age) + 
  scale_color_manual(name='Legend', breaks=c("Australia","Southeast Asia"), values=c('Australia'="#D55E00",'Southeast Asia'="#0072B2"))+
  scale_fill_manual(name='Legend', breaks=c("Australia","Southeast Asia"), values=c('Australia'="#D55E00",'Southeast Asia'="#0072B2"))+
  labs(x = "Age", y = "Predicted growth (mm)") + 
  scale_x_continuous(limits=c(0,30), breaks = c(0,5,10,15,20,25,30), labels = c('0','5','10','15','20','25','30'))+
  scale_y_continuous(limits=c(0.0,0.8), breaks = c(0.000, 0.200, 0.400, 0.600, 0.800), labels = c('0.0','0.2','0.4','0.6','0.8'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = "white"),
        legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot7

plot8 <- ggplot() + 
  geom_line(aes(y = transfit, x = AAC, colour="Australia"), linewidth=0.5, data = Ljohnii_Australia_AAC) + 
  geom_ribbon(aes(y = transfit, x = AAC, ymin = translower, ymax = transupper, fill="Australia"), linewidth=0.4, alpha=0.3, colour=0, data = Ljohnii_Australia_AAC) +
  geom_line(aes(y = transfit, x = AAC, colour="Southeast Asia"), linewidth=0.5,data = Ljohnii_Asia_AAC) + 
  geom_ribbon(aes(y = transfit, x = AAC, ymin = translower, ymax = transupper, fill="Southeast Asia"), linewidth=0.4, alpha=0.3, colour=0, data = Ljohnii_Asia_AAC) + 
  scale_color_manual(name='Legend', breaks=c("Australia","Southeast Asia"), values=c('Australia'="#D55E00",'Southeast Asia'="#0072B2"))+
  scale_fill_manual(name='Legend', breaks=c("Australia","Southeast Asia"), values=c('Australia'="#D55E00",'Southeast Asia'="#0072B2"))+
  labs(x = "AAC", y = "Predicted growth (mm)") + 
  scale_x_continuous(limits=c(0,30), breaks = c(0,5,10,15,20,25,30), labels = c('0','5','10','15','20','25','30'))+
  scale_y_continuous(limits=c(0.0,0.8), breaks = c(0.000, 0.200, 0.400, 0.600, 0.800), labels = c('0.0','0.2','0.4','0.6','0.8'))+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        legend.key = element_rect(colour = "white"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot8

ggpubr::ggarrange(plot7, plot8, labels = c("A", "B"))

#################################################################################################################################
#################################################################################################################################

#### Incorporating environmental variables ####
### Large-scale climate indices ####
## Pacific Decadal Oscillation: https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/ersst.v5.pdo.dat
## ReadME file for PDO: https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/index/Readme 

## Nino 4 ENSO Index: https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii 
## More information from https://www.cpc.ncep.noaa.gov/data/indices/ 
## Note that Nino 4 values were extracted and converted into the table format. 

## Dipole Mode Index: https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data 
## More information from https://psl.noaa.gov/gcos_wgsp/Timeseries/DMI/ 

## Download date: 08 Feb 2023

## Pacific Decadal Oscillation
pdo_data <- readLines('../../Climate/ersst.v5.pdo.dat.txt')
pdo_data <- as.data.frame(pdo_data, stringsAsFactors = F, header = F)
pdo_data <- cSplit(pdo_data, names(pdo_data)," ")
pdo_data <- pdo_data[-1,]
pdo_data <- as.data.frame(pdo_data)
names(pdo_data) <- NULL

names(pdo_data) <- pdo_data[1,]
pdo_data <- pdo_data[-1,]
pdo_data <- pdo_data %>% 
  mutate_if(is.character, as.numeric) 
pdo_data <- pdo_data[pdo_data$Year >= '1975',] # filter based on the minimum FishYear value of 1977, and we allow for a one year lag period
pdo_data <- pdo_data %>% as.data.frame(row.names = 1:nrow(.))

pdo_data <- pdo_data %>% 
  gather(Month, PDO, Jan:Dec) %>% 
  arrange(Year)
pdo_data$Month <- match(pdo_data$Month, month.abb) 
pdo_data$Day <- 1 # set arbitrary day
pdo_data <- pdo_data %>% mutate(ymd = make_date(Year, Month, Day))

## Remove cell values from the future
pdo_data$PDO[pdo_data$PDO == 99.99] <- NA
pdo_data <- na.omit(pdo_data)
setcolorder(pdo_data, c('ymd', 'Year', 'Month', 'Day', 'PDO'))

## El Nino Southern Oscillation
all_enso <- read.table("../../Climate/ersst5.nino.mth.91-20.ascii.txt", header = T)

## El Nino Southern Oscillation 4
enso_data <- all_enso %>% 
  dplyr::select(YR, MON, NINO4) %>% 
  rename(Year= YR, Month = MON) %>% 
  arrange(Year)
enso_data$Day <- 1
enso_data <- enso_data %>% mutate(ymd = make_date(Year, Month, Day))
enso_data <- enso_data[enso_data$Year >= '1975',] # filter based on the minimum FishYear value of 1977, and we allow for a one year lag period
enso_data <- enso_data %>% as.data.frame(row.names = 1:nrow(.))
enso_data <- na.omit(enso_data)
setcolorder(enso_data, c('ymd', 'Year', 'Month', 'Day', 'NINO4'))

## Dipole Mode Index
dmi_data <- readLines('../../Climate/dmi.had.long.data.txt')
dmi_data <- dmi_data[-c(1, 155:160)]
dmi_data <- as.data.frame(dmi_data,stringsAsFactors = F, header = T)
dmi_data <- cSplit(dmi_data, names(dmi_data)," ")
dmi_data <- dmi_data %>% mutate_if(is.character, as.numeric) 
dmi_data <- as.data.frame(dmi_data)

colnames(dmi_data) <- c("Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
dmi_data <- dmi_data %>% 
  gather(Month, DMI, Jan:Dec) %>% 
  arrange(Year)
dmi_data$Month <- match(dmi_data$Month, month.abb) 
dmi_data$Day <- '1'
dmi_data <- dmi_data %>% mutate(ymd = make_date(Year, Month, Day))
dmi_data <- dmi_data[dmi_data$Year >= '1975',] # filter based on the minimum FishYear value of 1977, and we allow for a one year lag period
dmi_data <- dmi_data %>% as.data.frame(row.names = 1:nrow(.))

dmi_data[dmi_data == -9999.000] <- NA
dmi_data <- na.omit(dmi_data)
setcolorder(dmi_data, c('ymd', 'Year', 'Month', 'Day', 'DMI'))

#################################################################################################################################
#################################################################################################################################

#### Download additional environmental variables ####
### Daily Optimum Interpolation Sea Surface Temperature Version 2.1 
## Dataset ID: ncdcOisst21Agg_LonPM180
## Dataset Title: SST, Daily Optimum Interpolation (OI), AVHRR Only, Version 2.1, Final, Global,0.25°, 1981-present, Lon+/-180
## Institution: NOAA NCEI
rerddap::info("ncdcOisst21Agg_LonPM180") 

## Defining spatial boundaries 
OISST_sub_dl <- function(time_df){
  OISST_dat <- griddap(datasetx = "ncdcOisst21Agg_LonPM180", 
                       url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                       time = c(time_df$start, time_df$end), 
                       zlev = c(0, 0),
                       latitude = c(-21.5, 7.5), 
                       longitude = c(90.5, 130.5), 
                       fields = "sst")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(date = time, temp = sst) %>% 
    dplyr::select(longitude, latitude, date, temp) %>% 
    na.omit()
} 

## Defining temporal boundaries 
dl_years1 <- data.frame(date_index = 1:5,
                        start = as.Date(c("1981-09-01", "1990-01-01", 
                                          "2000-01-01", "2010-01-01", 
                                          "2020-01-01")),
                        end = as.Date(c("1989-12-31", "1999-12-31", 
                                        "2009-12-31", "2019-12-31", 
                                        "2022-12-31"))
)

## The following lines of code extract data from NOAA ERDDAP Optimum Interpolation SST. 
oisst_dat <- dl_years1 %>% 
  dplyr::group_by(date_index) %>% 
  dplyr::group_modify(~OISST_sub_dl(.x)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(longitude, latitude, date, temp)

### Monthly Precipitation 
## Dataset ID: chirps20GlobalMonthlyP05_Lon0360
## Dataset Title: CHIRPS Version 2.0, Precipitation, Global, 0.05°, Monthly, 1981-present, Lon0360
## Institution: UCSB Climate Hazards Group
rerddap::info("chirps20GlobalMonthlyP05_Lon0360") 

## Defining spatial boundaries 
precip_sub_dl <- function(time_df){
  precip_dat <- griddap(datasetx = "chirps20GlobalMonthlyP05_Lon0360", 
                        url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                        time = c(time_df$start, time_df$end), 
                        latitude = c(-21.5, 7.5), 
                        longitude = c(90.5, 130.5), 
                        fields = "precip")$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(date = time) %>% 
    na.omit()
} 

## Defining temporal boundaries 
dl_years2 <- data.frame(date_index = 1:5,
                        start = as.Date(c("1981-01-01", "1990-01-01", 
                                          "2000-01-01", "2010-01-01", 
                                          "2020-01-01")),
                        end = as.Date(c("1989-12-31", "1999-12-31", 
                                        "2009-12-31", "2019-12-31", 
                                        "2022-12-01"))
)

## The following lines of code extract monthly precipitation data from UCSB Climate Hazards Group
precip_dat <- dl_years2 %>% 
  dplyr::group_by(date_index) %>% 
  dplyr::group_modify(~precip_sub_dl(.x)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(longitude, latitude, date, precip)

### ESA CCI Ocean Colour Product
## Dataset ID: pmlEsaCCI50OceanColorMonthly_Lon0360
## Dataset Title: ESA CCI Ocean Colour Product (CCI ALL-v5.0-MONTHLY), 0.04166666°, 1997-present, Lon0360 
## Institution: Plymouth Marine Laboratory 
rerddap::info("pmlEsaCCI50OceanColorMonthly_Lon0360") 

## Defining spatial boundaries 
chla_sub_dl <- function(time_df){
  chla_dat <- griddap(datasetx = "pmlEsaCCI50OceanColorMonthly_Lon0360", 
                      url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                      time = c(time_df$start, time_df$end), 
                      latitude = c(-21.5, 7.5), 
                      longitude = c(90.5, 130.5),
                      fields = c("chlor_a"))$data %>% 
    mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>% 
    dplyr::rename(date = time) %>% 
    na.omit()
} 

## Defining temporal boundaries 
dl_years3 <- data.frame(date_index = 1:6,
                        start = as.Date(c("1997-09-04", "2000-01-01", "2005-01-01", 
                                          "2010-01-01", "2015-01-01", "2020-01-01")),
                        end = as.Date(c("1999-12-31", "2004-12-31", "2009-12-31", 
                                        "2014-12-31", "2019-12-31", "2021-12-01"))
)

## The following lines of code extract monthly sea colour data from Plymouth Marine Laboratory 
chla_dat <- dl_years3 %>% 
  dplyr::group_by(date_index) %>% 
  dplyr::group_modify(~chla_sub_dl(.x)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(longitude, latitude, date, chlor_a) 

#################################################################################################################################
#################################################################################################################################

### Daily Optimum Interpolation Sea Surface Temperature Version 2.1 
## Data exploration
setDT(oisst_dat)
head(oisst_dat, 3) ; str(oisst_dat)
colnames(oisst_dat)
oisst_dat$longitude <- as.numeric(oisst_dat$longitude)
oisst_dat$latitude <- as.numeric(oisst_dat$latitude)

oisst_dat <- oisst_dat[ !(longitude < 97.5 | latitude > 6.0), ] # reduce the size of the bounding box

oisst_dat <- oisst_dat %>% # To amend the latitude and longitude boxes for the region maps
  dplyr::mutate(Region = case_when(latitude <= 5.0 & latitude >= 0.5 & longitude >= 105.0 & longitude <= 109.7 ~ "Riau Archipelago",
                                   latitude <= 1.5 & latitude >= 0.35 & longitude >= 104.0 & longitude <= 105.0 ~ "Riau Archipelago",
                                   
                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Malacca Straits",
                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Malacca Straits", 
                                   
                                   latitude <= -3.5 & latitude >= -6.0 & longitude >= 117.0 & longitude <= 119.5 ~ "Eastern Makassar",
                                   
                                   latitude <= -0.7 & latitude >= -5.0 & longitude >= 103.5 & longitude <= 108.5 ~ "South Eastern Sumatra",
                                   
                                   latitude <= -17.1 & latitude >= -21.5 & longitude >= 114.5 & longitude <= 116.5 ~ "Pilbara",
                                   latitude <= -15.0 & latitude >= -21.5 & longitude >= 116.5 & longitude <= 120.0 ~ "Pilbara",
                                   
                                   latitude <= -14.0 & latitude >= -20.0 & longitude >= 120.0 & longitude <= 122.0 ~ "Kimberley",
                                   latitude <= -13.0 & latitude >= -20.0 & longitude >= 122.0 & longitude <= 123.0 ~ "Kimberley",
                                   latitude <= -12.25 & latitude >= -17.0 & longitude >= 123.0 & longitude <= 126.0 ~ "Kimberley",
                                   latitude <= -11.25 & latitude >= -15.0 & longitude >= 126.0 & longitude <= 129.0 ~ "Kimberley",
                                   
                                   FALSE ~ "others" # Creates NA files 
  ))
oisst_data <- na.omit(oisst_dat) # remove the additional NA files
oisst_data <- as.data.frame(oisst_data)
head(oisst_data, 3) ; str(oisst_data)

### Monthly Precipitation 
## Data exploration
setDT(precip_dat)
head(precip_dat, 3) ; str(precip_dat)
colnames(precip_dat)
precip_dat$date <- as.Date(precip_dat$date)
precip_dat$longitude <- as.numeric(precip_dat$longitude)
precip_dat$latitude <- as.numeric(precip_dat$latitude)

precip_dat <- precip_dat %>% # To amend the latitude and longitude boxes for the region maps
  dplyr::mutate(Region = case_when(latitude <= 5.0 & latitude >= 0.5 & longitude >= 105.0 & longitude <= 109.7 ~ "Riau Archipelago",
                                   latitude <= 1.5 & latitude >= 0.35 & longitude >= 104.0 & longitude <= 105.0 ~ "Riau Archipelago",
                                   
                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Malacca Straits",
                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Malacca Straits", 
                                   
                                   latitude <= -3.5 & latitude >= -6.0 & longitude >= 117.0 & longitude <= 119.5 ~ "Eastern Makassar",
                                   
                                   latitude <= -0.7 & latitude >= -5.0 & longitude >= 103.5 & longitude <= 108.5 ~ "South Eastern Sumatra",
                                   
                                   latitude <= -17.1 & latitude >= -21.5 & longitude >= 114.5 & longitude <= 116.5 ~ "Pilbara",
                                   latitude <= -15.0 & latitude >= -21.5 & longitude >= 116.5 & longitude <= 120.0 ~ "Pilbara",
                                   
                                   latitude <= -14.0 & latitude >= -20.0 & longitude >= 120.0 & longitude <= 122.0 ~ "Kimberley",
                                   latitude <= -13.0 & latitude >= -20.0 & longitude >= 122.0 & longitude <= 123.0 ~ "Kimberley",
                                   latitude <= -12.25 & latitude >= -17.0 & longitude >= 123.0 & longitude <= 126.0 ~ "Kimberley",
                                   latitude <= -11.25 & latitude >= -15.0 & longitude >= 126.0 & longitude <= 129.0 ~ "Kimberley",
                                   
                                   FALSE ~ "others" # Creates NA files 
  ))
sum(is.na(precip_dat))
precip_data <- na.omit(precip_dat) # remove the additional NA files
precip_data <- as.data.frame(precip_data)
head(precip_data, 3) ; str(precip_data)

#################################################################################################################################
#################################################################################################################################

## Decided to remove the Chlorophyll A data set from further analyses due to the short temporal resolution of the data. 
## The inclusion of the data points causes the data frame to be truncated from 3000 values to <2000 values. 

### ESA CCI Ocean Colour Product
## Data exploration
#setDT(chla_dat)
#head(chla_dat, 3) ; str(chla_dat)
#colnames(chla_dat)
#chla_dat$date <- as.Date(chla_dat$date)
#chla_dat$longitude <- as.numeric(chla_dat$longitude)
#chla_dat$latitude <- as.numeric(chla_dat$latitude)

#chla_dat <- chla_dat %>% # To amend the latitude and longitude boxes for the region maps
#  dplyr::mutate(Region = case_when(latitude <= 5.0 & latitude >= 0.5 & longitude >= 105.0 & longitude <= 109.7 ~ "Riau Archipelago",
#                                   latitude <= 1.5 & latitude >= 0.35 & longitude >= 104.0 & longitude <= 105.0 ~ "Riau Archipelago",
#                                   
#                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Malacca Straits",
#                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Malacca Straits", 
#                                   
#                                   latitude <= -3.5 & latitude >= -6.0 & longitude >= 117.0 & longitude <= 119.5 ~ "Eastern Makassar",
#                                   
#                                   latitude <= -0.7 & latitude >= -5.0 & longitude >= 103.5 & longitude <= 108.5 ~ "South Eastern Sumatra",
#                                   
#                                   latitude <= -17.1 & latitude >= -21.5 & longitude >= 114.5 & longitude <= 116.5 ~ "Pilbara",
#                                   latitude <= -15.0 & latitude >= -21.5 & longitude >= 116.5 & longitude <= 120.0 ~ "Pilbara",
#                                   
#                                   latitude <= -14.0 & latitude >= -20.0 & longitude >= 120.0 & longitude <= 122.0 ~ "Kimberley",
#                                   latitude <= -13.0 & latitude >= -20.0 & longitude >= 122.0 & longitude <= 123.0 ~ "Kimberley",
#                                   latitude <= -12.25 & latitude >= -17.0 & longitude >= 123.0 & longitude <= 126.0 ~ "Kimberley",
#                                   latitude <= -11.25 & latitude >= -15.0 & longitude >= 126.0 & longitude <= 129.0 ~ "Kimberley",
#                                   
#                                  FALSE ~ "others" # Creates NA files 
#  ))
#sum(is.na(chla_dat))
#chla_data <- na.omit(chla_dat) # remove the additional NA files, note that we are working with the raw files 
#chla_data <- as.data.frame(chla_data)
#head(chla_data, 3) ; str(chla_data)

#################################################################################################################################
#################################################################################################################################

## Save as .csv file
#data.table::fwrite(oisst_data, file = "../../Climate/oisst_data.csv") 
#data.table::fwrite(precip_data, file = "../../Climate/precip_data.csv") 
#data.table::fwrite(chla_data, file = "../../Climate/chla_data.csv") 

## Remove large environmental datasets 
rm(oisst_dat)
rm(precip_dat)
#rm(chla_dat)

##save.image("lutjanus_johnii.RData")
##save.image("lutjanus_johnii_backup.RData")

#################################################################################################################################
################################################################################################################################# 

years <- c(1950:2023)

### Daily Optimum Interpolation Sea Surface Temperature Version 2.1 
oisst_data <- fread("../../Climate/oisst_data.csv")

oisst_data <- as.data.table(oisst_data)
## Amend the Region variable to Site to align it to individual locations. 
oisst_data <- oisst_data %>%
  rename(Site = Region)

## Create year, month, and day variables
head(oisst_data, 3) ; str(oisst_data)
oisst_data$date <- as.Date(oisst_data$date)
oisst_data$Site[oisst_data$Site == "Eastern Makassar"] <- "Sulawesi"
oisst_data$Year <- format(as.Date(oisst_data$date), "%Y")
oisst_data$Month <- format(as.Date(oisst_data$date), "%m")
oisst_data$Day <- format(as.Date(oisst_data$date), "%d")

oisst_data$FishYear <- as.character(as.factor(cut(oisst_data$date, breaks=as.Date(paste(years,"-09-01",sep="")), labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))))
setDT(oisst_data)[, paste0("FishYear", 1:2) := tstrsplit(FishYear, "/")]
oisst_data$FishYear <- as.numeric(as.character(oisst_data$FishYear1)) 
oisst_data[, c("FishYear1","FishYear2"):=NULL]
oisst_data$LagYear <- oisst_data$FishYear +1 
setcolorder(oisst_data, c('date','FishYear','LagYear','Year','Month','Day','longitude','latitude','Site','temp'))

### Create additional grouping variables
oisst_data$site_year <- as.factor(paste(oisst_data$Site, oisst_data$Year, sep = "_")) 
oisst_data$site_fyear <- as.factor(paste(oisst_data$Site, oisst_data$FishYear, sep = "_")) # note that this is fish year
oisst_data$site_month_year <- as.factor(paste(oisst_data$Site, oisst_data$Month, oisst_data$Year, sep = "_")) 
setcolorder(oisst_data, c('date','FishYear','LagYear','Year','Month','Day','longitude','latitude','Site','site_year','site_fyear','site_month_year','temp'))

## Monthly site specific sea surface temperature
oisst_month_data <- oisst_data %>% 
  dplyr::select(c("date","Year","Month","Day","Site","temp"))

### Precipitation 
precip_data <- fread("../../Climate/precip_data.csv")

precip_data <- as.data.table(precip_data)
## Amend the Region variable to Site to align it to individual locations. 
precip_data <- precip_data %>%
  rename(Site = Region)

## Create year, month, and day variables
head(precip_data, 3) ; str(precip_data)
precip_data$date <- as.Date(precip_data$date)
precip_data$Site[precip_data$Site == "Eastern Makassar"] <- "Sulawesi"
precip_data$Year <- format(as.Date(precip_data$date), "%Y")
precip_data$Month <- format(as.Date(precip_data$date), "%m")
precip_data$Day <- "1"

## Monthly rainfall data frame
precip_month_data <- precip_data %>% 
  dplyr::select(c("date","Year","Month","Day","Site","precip"))

precip_data$FishYear <- as.character(as.factor(cut(precip_data$date, breaks=as.Date(paste(years,"-09-01",sep="")), labels=paste(years[-length(years)], years[-length(years)]+1,sep="/"))))
setDT(precip_data)[, paste0("FishYear", 1:2) := tstrsplit(FishYear, "/")]
precip_data$FishYear <- as.numeric(as.character(precip_data$FishYear1)) 
precip_data[, c("FishYear1","FishYear2"):=NULL]
precip_data$LagYear <- precip_data$FishYear +1 
setcolorder(precip_data, c('date','FishYear','LagYear','Year','Month','longitude','latitude','Site','precip'))

### Create additional grouping variables
precip_data$site_year <- as.factor(paste(precip_data$Site, precip_data$Year, sep = "_")) 
precip_data$site_fyear <- as.factor(paste(precip_data$Site, precip_data$FishYear, sep = "_")) # note that this is fish year
setcolorder(precip_data, c('date','FishYear','LagYear','Year','Month','longitude','latitude','Site','site_year','site_fyear','precip'))

## Monthly site specific precipitation
precip_month_data <- precip_data %>% 
  dplyr::select(c("date","Year","Month","Day","Site","precip"))

## Decided to remove the Chlorophyll A data set from further analyses due to the short temporal resolution of the data. 
### Chlorophyll A 
#chla_data <- fread("../../Climate/chla_data.csv")

## Create year, month, and day variables
#setDT(chla_data)
#head(chla_data, 3) ; str(chla_data)
#chla_data$date <- as.Date(chla_data$date)
#chla_data$Region[chla_data$Region == "Eastern Makassar"] <- "Sulawesi"
#chla_data$Year <- format(as.Date(chla_data$date), "%Y")
#chla_data$Month <- format(as.Date(chla_data$date), "%m")

#chla_data$FishYear <- as.character(as.factor(cut(chla_data$date, breaks=as.Date(paste(years,"-09-01",sep="")), labels=paste(years[-length(years)], years[-length(years)]+1,sep="/"))))
#setDT(chla_data)[, paste0("FishYear", 1:2) := tstrsplit(FishYear, "/")]
#chla_data$FishYear <- as.numeric(as.character(chla_data$FishYear1)) 
#chla_data[, c("FishYear1","FishYear2"):=NULL]
#chla_data$LagYear <- chla_data$FishYear +1 
#setcolorder(chla_data, c('date','FishYear','LagYear','Year','Month','longitude','latitude','Region','chlor_a'))

### Create additional grouping variables
#chla_data$region_year <- as.factor(paste(chla_data$Region, chla_data$Year, sep = "_")) 
#chla_data$region_fyear <- as.factor(paste(chla_data$Region, chla_data$FishYear, sep = "_")) # note that this is fish year
#setcolorder(chla_data, c('date','FishYear','LagYear','Year','Month','longitude','latitude','Region','region_year','region_fyear','chlor_a'))

#################################################################################################################################
################################################################################################################################# 

## Combine the regional climate indices into a single data frame 
env_data <- merge(pdo_data, enso_data[,c('ymd', 'NINO4')], by = c('ymd')) %>%
  merge(dmi_data[,c('ymd', 'DMI')], by = c('ymd')) 

## Format the local environmental variables into the same format
oisst_month_data <- oisst_month_data %>% 
  filter(Site == "Riau Archipelago" | Site == "Malacca Straits" | Site == "Kimberley") %>% 
  rename(ymd = date) 
oisst_month_data[, temp := mean(temp) , by = .(Year, Month, Site)][] 
oisst_month_data <- oisst_month_data %>%
  distinct(temp, Year, Month, .keep_all = TRUE) 

precip_month_data <- precip_month_data %>% 
  filter(Site == "Riau Archipelago" | Site == "Malacca Straits" | Site == "Kimberley") %>% 
  rename(ymd = date) 
precip_month_data[, precip := mean(precip) , by = .(Year, Month, Site)][] 
precip_month_data <- precip_month_data %>%
  distinct(precip, Year, Month, .keep_all = TRUE) 

env_data <- merge(env_data, oisst_month_data[,c('ymd','Site','temp')], by = c('ymd')) 
env_data <- merge(env_data, precip_month_data[,c('ymd','Site','precip')], by = c('ymd','Site')) 

env_data$ymd <- ymd(env_data$ymd)
env_data$FishYear <- as.character(as.factor(cut(env_data$ymd, breaks=as.Date(paste(years,"-09-01",sep="")), labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))))
setDT(env_data)[, paste0("FishYear", 1:2) := tstrsplit(FishYear, "/")]
env_data$FishYear <- as.numeric(as.character(env_data$FishYear1)) 
env_data[, c("FishYear1","FishYear2"):=NULL]

colnames(env_data)
setcolorder(env_data, c('ymd','Site','Year','Month','Day','FishYear','PDO','NINO4','DMI','temp','precip'))
ggpairs(env_data, 7:11) # Looks good, no discerable pattern present

## Variance Inflation Factor files provided by Highland Statistics Ltd.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}

### Variance inflation factor (VIF)
# The threshold for an ideal VIF varies between literature books. 
# Some books claiming that 10 is the threshold while others state that it is 2. 
corvif(env_data[, c(7:11)]) # The VIF values vary between environmental covariates, for both regional and large scale variables. 

#Mypairs
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 6) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) points(x, y, 
                                             pch = 16, cex = 0.8, 
                                             col = gray(0.1)))
  #print(P)
}

Mypairs(env_data[, c(7:11)]) # looks good

#################################################################################################################################
#################################################################################################################################

### Plot the data to observe trends in the environmental data #### 
## Annual optimum interpolation sea surface temperature
temp_plot <- env_data %>% 
  ggplot(aes(y=temp, x=ymd, colour=Site))+ 
  geom_line(linewidth=0.2)+ geom_point()+ 
  geom_smooth(method = 'loess')+
  facet_wrap(~fct_relevel(Site), ncol=3)+
  labs(y= 'Sea surface temperature (˚C)')+ 
  theme_bw()
temp_plot

## Annual precipitation
rain_plot <- env_data %>% 
  ggplot(aes(y=precip, x=ymd, colour=Site))+ 
  geom_line(linewidth=0.2)+ geom_point()+ 
  geom_smooth(method = 'loess')+
  facet_wrap(~fct_relevel(Site), ncol=3)+
  labs(y= 'Rainfall (mm)')+ 
  theme_bw()
rain_plot

## Pacific Decadal Oscillation index
pdo_plot <- env_data %>% 
  ggplot(aes(y=PDO, x=ymd))+ 
  geom_line()+ 
  labs(y= 'Pacific Decadal Oscillation')+ 
  theme_bw()
pdo_plot

## El Niño Southern Oscillation 4 index
enso_plot <- env_data %>% 
  ggplot(aes(y=NINO4, x=ymd))+ 
  geom_line()+ 
  labs(y= 'ENSO (Niño 4 Index)')+ 
  theme_bw()
enso_plot

## Dipole Mode index
dmi_plot <- env_data %>% 
  ggplot(aes(y=DMI, x=ymd))+ 
  geom_line()+ 
  labs(y= 'Dipole Mode Index')+ 
  theme_bw()
dmi_plot

## Scale all environmental variables
env_data$s.temp <- scale(env_data$temp)
env_data$s.precip <- scale(env_data$precip)

env_data$s.pdo <- scale(env_data$PDO)
env_data$s.enso <- scale(env_data$NINO4)
env_data$s.dmi <- scale(env_data$DMI)

## Investigate the change in mean DMI across months
dmi_month_plot <- env_data %>% 
  subset(select = c(ymd,Year,Month,DMI,s.dmi)) %>% 
  distinct(ymd, .keep_all = TRUE) 
dmi_month_plot <- dmi_month_plot %>% 
  group_by(Month) %>% 
  summarise(Mean=mean(DMI), SD=sd(DMI))
dmi_month_plot$lower <- dmi_month_plot$Mean - dmi_month_plot$SD
dmi_month_plot$upper <- dmi_month_plot$Mean + dmi_month_plot$SD

dmi_month_plot %>% 
  ggplot(aes(y=Mean, x=Month))+ 
  geom_line(linewidth=1)+
  geom_point(size=2.5)+
  scale_x_continuous(limits=c(1,12), breaks=c(1,2,3,4,5,6,7,8,9,10,11,12))+
  labs(y="Mean DMI")+ 
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
# Note that the upper and lower SD can also be included but it was not included due to the large error bars. 

ggpubr::ggarrange(temp_plot, rain_plot, labels = c("A", "B"), nrow=2)

ggpubr::ggarrange(pdo_plot, enso_plot, dmi_plot, labels = c("A", "B", "C"), nrow=3)

## Check data frame for missing data points
colSums(is.na(env_data)) # note that this is repeated for all increments, not unique values
visdat::vis_miss(env_data) 

## Check for linearity and correlations among local and regional environmental covariates 
ggpairs(env_data, 7:11) # Linearity between environmental covariates across all Regions
ggpairs(env_data, ggplot2::aes(colour=Site),7:11) # Linearity between covariates separated between Sites
Mypairs(env_data[, c(7:11)]) 
## There appears to be a relationship between sst and precipitation, where cool sites have lower rainfall compared to warmer sites which have higher rainfall quantities. 
## Given that they are different environmental variables, we can use them. 

## Check the Variance Inflation Factor
corvif(env_data[, c(7:11)]) # All less than 5

Australia_env_data <- env_data %>% filter(Site == "Kimberley")
Asia_env_data <- env_data %>% filter(Site == "Malacca Straits" | Site == "Riau Archipelago")

corvif(Australia_env_data[, c(7:11)]) 
corvif(Asia_env_data[, c(7:11)]) 

## Check the environmental data frame
env_data <- as.data.frame(env_data)
env_data

#################################################################################################################################
#################################################################################################################################

### Prepare the data for the exploratory sliding window analysis
Ljohnii_Australia <- Ljohnii_Australia %>%
  rename(capture_year = Year, capture_month = Month, capture_day = Day)

## Dummy variable for the alignment of environmental variables
Ljohnii_Australia$Month <- "9" 
Ljohnii_Australia$Month <- as.numeric(Ljohnii_Australia$Month)
Ljohnii_Australia <- Ljohnii_Australia %>% mutate(ymd = make_date(FishYear, Month)) 

## Preliminary analyses of the data indicated that we have biological data which starts in 1976. However, our environmental data only starts in September 1981. 
## Truncate the data to 1982 to ensure that we have environmental data coverage for all increment years. 
Ljohnii_Australia <- Ljohnii_Australia %>% filter(ymd >= '1983-09-01')
Ljohnii_Australia # Note the new data frame name
## Note that we truncated the data to 1983 to allow us to explore if growth was affected by environmental conditions in 1981 (two years).

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3a <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Australia, REML= F)

## Applying the sliding window to the environmental data ####
M3a_sliding_window <- slidingwin(baseline = M3a, 
                                 xvar = list(sst = env_data$s.temp,
                                             rain = env_data$s.precip,
                                             pdo = env_data$s.pdo,
                                             enso = env_data$s.enso,
                                             dmi = env_data$s.dmi),  
                                 type = "absolute",                                 # "absolute" or "relative", relative climate window=number of days before each biological record is measured; absolute climate window=number of days before a set point (refday) in time.
                                 range = c(23, 0),                                  # the furthest and closest number of time intervals (set by cinterval)
                                 stat = c("mean"),                                  # aggregate statistics used
                                 func = c("lin"),                                   # functions used to fit the climate variable. Can be linear ("lin"), quadratic ("quad"), cubic ("cub"), inverse ("inv") or log ("log")
                                 refday = c(1, 9),                                  # day and month respectively of the year from which the absolute window analysis will start
                                 cinterval = "month",                               # resolution at which climate window analysis will be conducted. 
                                 cdate = env_data$ymd,                              # climate date variable (dd/mm/yyyy)
                                 bdate = Ljohnii_Australia$ymd,
                                 cohort = Ljohnii_Australia$FishYear,       # cohort is required because our environmental data spans across two years (see Advanced Vignette for more information)
                                 spatial = list(Ljohnii_Australia$Site, env_data$Site))

## As there is little temporal variation as to when the increment was completed, I used the absolute sliding window to analyse the data. 
## This means that there is a fixed reference time for all climate windows. 
## The climate year is important for the alignment of biological data, but the month is irrelevant for the biological data, as we assumed that the opaque zone completion period is consistent every year.
## On the other hand, bdate represents the biological date variable for the alignment of cdate. 
## The day and month of bdate is irrelevant and is only required for the alignment of cdate. 

## Note: Delta AICc is generated in comparison with the base null model M3. 

### Sliding window data ####
## The 'combos' function provides an overview of the best fitted climate window, the delta AICc value, the start and end window date, and slope of the best window. 
M3a_sliding_window$combos 

### Applying the random window approach 
## The random window is repeated 1000 times to estimate the probability that a given result represents a false positive. 
## It randomises a given dataset (i.e., removes any climate signal) and conducts a sliding window analysis to extract a value of delta AIC. 
M3a_random_window <- randwin(baseline = M3a, 
                             repeats = 1000, # Randomisation of 1000 times 
                             xvar = list(sst = env_data$s.temp,
                                         rain = env_data$s.precip,
                                         pdo = env_data$s.pdo,
                                         enso = env_data$s.enso,
                                         dmi = env_data$s.dmi), 
                             type = "absolute",                       
                             range = c(23, 0),                      
                             stat = c("mean"),                          
                             func = c("lin"),                         
                             refday=c(1, 9),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Ljohnii_Australia$ymd,
                             cohort = Ljohnii_Australia$FishYear,
                             spatial = list(Ljohnii_Australia$Site, env_data$Site),
                             window= "sliding")

## There are warning messages when running the random window. 
## Fixed effect model matrix is rank deficient so dropping 1 column / coefficient. 
## This error message does not occur with the earlier mixed effects models constructed (i.e., M2), appearing only in M3a, likely suggesting that the climate data is not comprehensive enough for the sliding window approach. 

## As the reference date does not matter in this scenario, the absolute sliding window indicates that annual growth is dependent on the average annual temperature at depth conditions across the Gascoyne region over the Fish Year. 

## Calculate the p value to establish significance 
pvalue(datasetrand = M3a_random_window[[1]], dataset = M3a_sliding_window[[1]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3a_random_window[[2]], dataset = M3a_sliding_window[[2]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3a_random_window[[3]], dataset = M3a_sliding_window[[3]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3a_random_window[[4]], dataset = M3a_sliding_window[[4]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3a_random_window[[5]], dataset = M3a_sliding_window[[5]]$Dataset, metric = "AIC") 

### Exploratory plots ####
## Note based on the warning that Pc will be overly conservative when the sample size is greater than 47. 
## Site-specific sea surface temperature
plotall(datasetrand = M3a_random_window[[1]],
        dataset = M3a_sliding_window[[1]]$Dataset, 
        bestmodel = M3a_sliding_window[[1]]$BestModel,
        bestmodeldata = M3a_sliding_window[[1]]$BestModelData)

## Site-specific precipitation  
plotall(datasetrand = M3a_random_window[[2]],
        dataset = M3a_sliding_window[[2]]$Dataset, 
        bestmodel = M3a_sliding_window[[2]]$BestModel,
        bestmodeldata = M3a_sliding_window[[2]]$BestModelData)

## Regional Pacific Decadal Oscillation
plotall(datasetrand = M3a_random_window[[3]],
        dataset = M3a_sliding_window[[3]]$Dataset, 
        bestmodel = M3a_sliding_window[[3]]$BestModel,
        bestmodeldata = M3a_sliding_window[[3]]$BestModelData)

## Regional El Niño Southern Oscillation
plotall(datasetrand = M3a_random_window[[4]],
        dataset = M3a_sliding_window[[4]]$Dataset, 
        bestmodel = M3a_sliding_window[[4]]$BestModel,
        bestmodeldata = M3a_sliding_window[[4]]$BestModelData)

## Regional Dipole Mode Index
plotall(datasetrand = M3a_random_window[[5]],
        dataset = M3a_sliding_window[[5]]$Dataset, 
        bestmodel = M3a_sliding_window[[5]]$BestModel,
        bestmodeldata = M3a_sliding_window[[5]]$BestModelData)

### Sliding window for rainfall data #### 
M3a_plot <- as.data.frame(M3a_sliding_window[[5]]$Dataset)

plot9 <- M3a_plot %>% 
  ggplot(aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+
  geom_tile(aes(fill=deltaAICc))+
  geom_tile(data = M3a_plot[c(1),], fill = NA, color = "white", linewidth = 1.5) +
  scale_x_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec\n(t-1)','Nov\n(t-1)','Oct\n(t-1)','Sep\n(t-1)','Aug\n(t-1)','Jul\n(t-1)','Jun\n(t-1)','May\n(t-1)','Apr\n(t-1)','Mar\n(t-1)','Feb\n(t-1)','Jan\n(t-1)','Dec\n(t-2)','Nov\n(t-2)','Oct\n(t-2)'), breaks=c(0:23), expand=c(0,0))+
  scale_y_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec (t-1)','Nov (t-1)','Oct (t-1)','Sep (t-1)','Aug (t-1)','Jul (t-1)','Jun (t-1)','May (t-1)','Apr (t-1)','Mar (t-1)','Feb (t-1)','Jan (t-1)','Dec (t-2)','Nov (t-2)','Oct (t-2)'), breaks=c(0:23), expand=c(0,0))+
  scale_fill_viridis(option = 'magma')+
  labs(x="Window Close", y="Window Open")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot9

### Distribution of AICc values plotted in a histogram
M3a_hist <- as.data.frame(M3a_random_window[[5]])

plot10 <- M3a_hist %>%
  ggplot(aes(x= deltaAICc))+
  geom_histogram(fill='grey20', bins = 50)+
  geom_hline(yintercept = 0, colour = 'grey80')+
  geom_vline(xintercept = -11.93, colour = 'red', linetype = 'dashed')+
  labs(x="Delta AICc", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot10

#################################################################################################################################
#################################################################################################################################

### Construct the optimal growth model for Lutjanus johnii from Northwestern Australia ####
## Our results indicate that growth was best explained by the IOD based on the p values.

## Site-specific sea surface temperature
head(arrange(M3a_sliding_window[[5]]$Dataset), 1)
plotdelta(dataset = M3a_sliding_window[[5]]$Dataset)

## Check the minimum Fish Year for samples from Northwestern Australia
## Ensure the annual growth band is explained by IOD conditions from the prior two years (Oct 1981 - Sept 1983).
min(Ljohnii_Australia$FishYear)

## Group 1
env_Australia1 <- env_data %>% 
  filter(Site == "Kimberley") %>% 
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>%
  subset(select = c(ymd, Site, Year, Month, Day, DMI, s.dmi, group, group1))

env_Australia1 <- env_Australia1 %>% 
  group_by(Site) %>% 
  mutate(group1 = consecutive_id(2 * floor(group1 / 2))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>% 
  group_by(Site, group1) %>% 
  mutate(group2 = row_number() - 1)

env_Australia_name1 <- env_Australia1 %>% 
  group_by(Site, group1) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/"))

env_Australia1 <- env_Australia1 %>% 
  filter(group2 == "9" | group2 == "8" | group2 == "7" | group2 == "6" | group2 == "5") %>%
  group_by(Site, group1) %>%
  summarise(mean_s.dmi = mean(s.dmi), mean_DMI = mean(DMI)) 

env_Australia1 <- merge(env_Australia1, env_Australia_name1)
env_Australia1 <- env_Australia1 %>%
  filter(!group1 == "20") 
env_Australia1

## Group 2 
env_Australia2 <- env_data %>% 
  filter(Site == "Kimberley") %>% 
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>%
  subset(select = c(ymd, Site, Year, Month, Day, DMI, s.dmi, group, group1))

env_Australia2 <- env_Australia2 %>% 
  group_by(Site) %>% 
  mutate(group2 = consecutive_id(2 * floor((group1 + 1) / 2))) %>% 
  mutate(group2 = as.numeric(as.factor(group2)) - 1) %>%
  group_by(Site, group2) %>% 
  mutate(group3 = row_number() - 1)
# Note the difference in code between group 1 and 2.

env_Australia_name2 <- env_Australia2 %>% 
  group_by(Site, group2) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/")) 

env_Australia2 <- env_Australia2 %>% 
  filter(group3 == "9" | group3 == "8" | group3 == "7" | group3 == "6" | group3 == "5") %>%
  group_by(Site, group2) %>%
  summarise(mean_s.dmi = mean(s.dmi), mean_DMI = mean(DMI)) 

env_Australia2 <- merge(env_Australia2, env_Australia_name2)
env_Australia2 <- env_Australia2 %>%
  filter(!group2 == "0") %>%
  rename(group1 = group2)
env_Australia2

env_Australia <- rbind(env_Australia1, env_Australia2)
env_Australia <- env_Australia %>%
  arrange(s_Year) %>%
  subset(select = -c(group1))
env_Australia <- env_Australia %>% 
  tidyr::separate(col="s_Year", sep=c(4), into=c("FishYear", "delim1"), remove=FALSE) %>%
  subset(select = -c(delim1))
env_Australia

### Construct the optimal growth model for Lutjanus johnii ###
Ljohnii_Australia_final <- merge(Ljohnii_Australia, env_Australia, by = c('Site','FishYear'))

M3a_final <- lmer(log(Width) ~ sAge + sAAC + 
                    mean_s.dmi + 
                    (sAge|SampleID) + (1|FishYear), 
                  Ljohnii_Australia_final, REML= T)

summary(M3a_final)
anova(M3a_final)
rsquared(M3a_final)

M3a_conf <- sim(M3a_final, 5000)
M3a_conf <- apply(M3a_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3a_conf

p <- predictorEffect("mean_s.dmi", (M3a_final)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto")) 

### Model fitted with local and regional environmental variables ###
## Linear DMI 
M3a_dmi <- as.data.frame (Effect (c('mean_s.dmi'), M3a_final, xlevels= list (mean_s.dmi = seq(-2.66, 2.11, by= 0.01)) )) 
M3a_mean_dmi <- mean(Ljohnii_Australia_final$mean_DMI) ; M3a_mean_dmi
M3a_sd_dmi <- sd(Ljohnii_Australia_final$mean_DMI) ; M3a_sd_dmi
M3a_dmi$dmi <- (M3a_dmi$mean_s.dmi * M3a_sd_dmi + M3a_mean_dmi)

M3a_dmi$transfit <- exp(M3a_dmi$fit)
M3a_dmi$transupper <- exp(M3a_dmi$upper)
M3a_dmi$translower <- exp(M3a_dmi$lower)

min(M3a_dmi$dmi) ; max(M3a_dmi$dmi)
min(Ljohnii_Australia_final$mean_DMI) ; max(Ljohnii_Australia_final$mean_DMI)

## Calculate the effect of the Indian Ocean Dipole on Lutjanus johnii growth
dmi_min <- M3a_dmi %>% filter(dmi == min(dmi)) %>% dplyr::select(transfit) 
dmi_max <- M3a_dmi %>% filter(dmi == max(dmi)) %>% dplyr::select(transfit)
dmi_unit <- (M3a_dmi %>% filter(dmi == max(dmi)) %>% dplyr::select(dmi)) - (M3a_dmi %>% filter(dmi == min(dmi)) %>% dplyr::select(dmi))
dmi_rate <- ((dmi_max - dmi_min) / (dmi_unit)) * 100
dmi_rate

### Plot the effect of DMI on Lutjanus johnii growth 
plot11 <- M3a_dmi %>%
  ggplot(aes(y = transfit, x = dmi)) +
  geom_line(linewidth=0.75) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.3, alpha = 0.3, colour = 0) + 
  scale_x_continuous(limits = c(-0.606, 0.400), breaks = c(-0.6, -0.4, -0.2, 0.0, 0.2, 0.4), labels = c('-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4')) +
  labs(y = 'Annual otolith growth (mm)', x = 'Mean IOD intensity') +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        legend.position.inside = c(0.8,0.125),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot11

ggpubr::ggarrange(plot9, plot11, labels = c("A", "B"))

#################################################################################################################################
#################################################################################################################################

### Mature Lutjanus johnii from Australia ####
Ljohnii_Australia_adult <- Ljohnii_Australia
Ljohnii_Australia_adult <- Ljohnii_Australia_adult %>% 
  filter(MaturityStage == "Mature")

M1a_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID), 
            Ljohnii_Australia_adult, REML= T)

M1b_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID), 
            Ljohnii_Australia_adult, REML= T)

M1c_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|FishYear), 
            Ljohnii_Australia_adult, REML= T)

M1d_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Australia_adult, REML= T)

M1e_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|Cohort), 
            Ljohnii_Australia_adult, REML= T)

M1f_adult <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Ljohnii_Australia_adult, REML= T)

bbmle::AICctab(M1a_adult,M1b_adult,M1c_adult,M1d_adult,M1e_adult,M1f_adult, base=T,logLik=T,weights=T) 

M1_adult <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
             Ljohnii_Australia_adult, REML= F)

options(na.action=na.fail)
M1_adult_dredge<-dredge(M1_adult, trace=2) 

subset(M1_adult_dredge)

M1_adult <- lmer(log(Width) ~ sAge + 
             (sAge|SampleID) + (1|FishYear), 
             Ljohnii_Australia_adult, REML= T)

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3b <- lmer(log(Width) ~ sAge + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Australia_adult, REML= F)

## Applying the sliding window to the environmental data ####
M3b_sliding_window <- slidingwin(baseline = M3b, 
                                 xvar = list(sst = env_data$s.temp,
                                             rain = env_data$s.precip,
                                             pdo = env_data$s.pdo,
                                             enso = env_data$s.enso,
                                             dmi = env_data$s.dmi),  
                                 type = "absolute",
                                 range = c(23, 0),    
                                 stat = c("mean"), 
                                 func = c("lin"),  
                                 refday = c(1, 9),
                                 cinterval = "month",  
                                 cdate = env_data$ymd,   
                                 bdate = Ljohnii_Australia_adult$ymd,
                                 cohort = Ljohnii_Australia_adult$FishYear,
                                 spatial = list(Ljohnii_Australia_adult$Site, env_data$Site))

M3b_sliding_window$combos 

M3b_random_window <- randwin(baseline = M3b, 
                             repeats = 1000, 
                             xvar = list(sst = env_data$s.temp,
                                         rain = env_data$s.precip,
                                         pdo = env_data$s.pdo,
                                         enso = env_data$s.enso,
                                         dmi = env_data$s.dmi), 
                             type = "absolute",                       
                             range = c(23, 0),                      
                             stat = c("mean"),                          
                             func = c("lin"),                         
                             refday=c(1, 9),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Ljohnii_Australia_adult$ymd,
                             cohort = Ljohnii_Australia_adult$FishYear,
                             spatial = list(Ljohnii_Australia_adult$Site, env_data$Site),
                             window= "sliding")

## Calculate the p value to establish significance 
pvalue(datasetrand = M3b_random_window[[1]], dataset = M3b_sliding_window[[1]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3b_random_window[[2]], dataset = M3b_sliding_window[[2]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3b_random_window[[3]], dataset = M3b_sliding_window[[3]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3b_random_window[[4]], dataset = M3b_sliding_window[[4]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3b_random_window[[5]], dataset = M3b_sliding_window[[5]]$Dataset, metric = "AIC") 

## Site-specific sea surface temperature
plotall(datasetrand = M3b_random_window[[1]],
        dataset = M3b_sliding_window[[1]]$Dataset, 
        bestmodel = M3b_sliding_window[[1]]$BestModel,
        bestmodeldata = M3b_sliding_window[[1]]$BestModelData)

## Site-specific precipitation  
plotall(datasetrand = M3b_random_window[[2]],
        dataset = M3b_sliding_window[[2]]$Dataset, 
        bestmodel = M3b_sliding_window[[2]]$BestModel,
        bestmodeldata = M3b_sliding_window[[2]]$BestModelData)

## Regional Pacific Decadal Oscillation
plotall(datasetrand = M3b_random_window[[3]],
        dataset = M3b_sliding_window[[3]]$Dataset, 
        bestmodel = M3b_sliding_window[[3]]$BestModel,
        bestmodeldata = M3b_sliding_window[[3]]$BestModelData)

## Regional El Niño Southern Oscillation
plotall(datasetrand = M3b_random_window[[4]],
        dataset = M3b_sliding_window[[4]]$Dataset, 
        bestmodel = M3b_sliding_window[[4]]$BestModel,
        bestmodeldata = M3b_sliding_window[[4]]$BestModelData)

## Regional Dipole Mode Index
plotall(datasetrand = M3b_random_window[[5]],
        dataset = M3b_sliding_window[[5]]$Dataset, 
        bestmodel = M3b_sliding_window[[5]]$BestModel,
        bestmodeldata = M3b_sliding_window[[5]]$BestModelData)

### Construct the optimal growth model for Lutjanus johnii ###
M3b_final <- lmer(log(Width) ~ sAge + 
                    (sAge|SampleID) + (1|FishYear), 
                  Ljohnii_Australia_adult, REML= T)

summary(M3b_final)
anova(M3b_final)
rsquared(M3b_final)

M3b_conf <- sim(M3b_final, 5000)
M3b_conf <- apply(M3b_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3b_conf

#################################################################################################################################
#################################################################################################################################

### Immature Lutjanus johnii #### 
Ljohnii_Australia_juvenile <- Ljohnii_Australia
Ljohnii_Australia_juvenile <- Ljohnii_Australia_juvenile %>% 
  filter(MaturityStage == "Immature")

M1a_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID), 
                  Ljohnii_Australia_juvenile, REML= T)

M1b_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID), 
                  Ljohnii_Australia_juvenile, REML= T)

M1c_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID) + (1|FishYear), 
                  Ljohnii_Australia_juvenile, REML= T)

M1d_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|FishYear), 
                  Ljohnii_Australia_juvenile, REML= T)

M1e_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID) + (1|Cohort), 
                  Ljohnii_Australia_juvenile, REML= T)

M1f_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|Cohort), 
                  Ljohnii_Australia_juvenile, REML= T)

bbmle::AICctab(M1a_juvenile,M1b_juvenile,M1c_juvenile,M1d_juvenile,M1e_juvenile,M1f_juvenile, base=T,logLik=T,weights=T) 

M1_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                      (sAge|SampleID) + (1|Cohort), 
                 Ljohnii_Australia_juvenile, REML= F)

options(na.action=na.fail)
M1_juvenile_dredge<-dredge(M1_juvenile, trace=2) 

subset(M1_juvenile_dredge)

M1_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                   (sAge|SampleID) + (1|Cohort), 
                 Ljohnii_Australia_juvenile, REML= T)

plot(M1_juvenile) 
qqnorm(residuals(M1_juvenile)) ; qqline(residuals(M1_juvenile)) 

summary(M1_juvenile) 
anova(M1_juvenile)
rsquared(M1_juvenile)

M1_juvenile_conf <- sim(M1_juvenile, 5000)
M1_juvenile_conf <- apply(M1_juvenile_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M1_juvenile_conf

p <- predictorEffect("sAge", (M1_juvenile)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))
p <- predictorEffect("sAAC", (M1_juvenile)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto"))

#### Cohort model for Australian samples ###
M1_juvenile_year <- ranef(M1_juvenile)$Cohort[,1]
M1_juvenile_year_se <-sqrt (attr(ranef(M1_juvenile, postVar=TRUE) [["Cohort"]], "postVar")[1,1,])

modelcohort <- data.frame(y=M1_juvenile_year)
modelcohort$upper <- (modelcohort$y + M1_juvenile_year_se)
modelcohort$lower <- (modelcohort$y - M1_juvenile_year_se)
modelcohort$cohort <- rownames(ranef(M1_juvenile)$Cohort)
modelcohort$cohort <- as.integer(modelcohort$cohort)
sample_depth_3 <- Ljohnii_Australia %>% 
  group_by(Cohort) %>% 
  summarise(N=n())
sample_depth_3$Cohort <- as.integer(as.character(sample_depth_3$Cohort))
colnames(sample_depth_3)[1] <- "cohort"
modelcohort <- modelcohort %>% 
  left_join(sample_depth_3) 
modelcohort

#### Plot relative annual otolith growth of Australian samples over time ###
p <- modelcohort %>% 
  ggplot(aes(y=y, x=cohort))+
  geom_point(size=1)+ 
  geom_line(linewidth=0.5)+
  geom_ribbon(aes(y=y, ymin=lower, ymax=upper), alpha=0.4)+
  scale_x_continuous(limits=c(1975,2023), breaks=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020), labels=c(1975,1980,1985,1990,1995,2000,2005,2010,2015,2020))+
  labs(y="Relative annual otolith growth", x="Cohort")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
p

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3c <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Ljohnii_Australia_juvenile, REML= F)

## Applying the sliding window to the environmental data ####
M3c_sliding_window <- slidingwin(baseline = M3c, 
                                 xvar = list(sst = env_data$s.temp,
                                             rain = env_data$s.precip,
                                             pdo = env_data$s.pdo,
                                             enso = env_data$s.enso,
                                             dmi = env_data$s.dmi),  
                                 type = "absolute",
                                 range = c(23, 0),    
                                 stat = c("mean"), 
                                 func = c("lin"),  
                                 refday = c(1, 9),
                                 cinterval = "month",  
                                 cdate = env_data$ymd,   
                                 bdate = Ljohnii_Australia_juvenile$ymd,
                                 cohort = Ljohnii_Australia_juvenile$FishYear,
                                 spatial = list(Ljohnii_Australia_juvenile$Site, env_data$Site))

M3c_sliding_window$combos 

M3c_random_window <- randwin(baseline = M3c, 
                             repeats = 1000, 
                             xvar = list(sst = env_data$s.temp,
                                         rain = env_data$s.precip,
                                         pdo = env_data$s.pdo,
                                         enso = env_data$s.enso,
                                         dmi = env_data$s.dmi), 
                             type = "absolute",                       
                             range = c(23, 0),                      
                             stat = c("mean"),                          
                             func = c("lin"),                         
                             refday=c(1, 9),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Ljohnii_Australia_juvenile$ymd,
                             cohort = Ljohnii_Australia_juvenile$FishYear,
                             spatial = list(Ljohnii_Australia_juvenile$Site, env_data$Site),
                             window= "sliding")

## Calculate the p value to establish significance 
pvalue(datasetrand = M3c_random_window[[1]], dataset = M3c_sliding_window[[1]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3c_random_window[[2]], dataset = M3c_sliding_window[[2]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3c_random_window[[3]], dataset = M3c_sliding_window[[3]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3c_random_window[[4]], dataset = M3c_sliding_window[[4]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3c_random_window[[5]], dataset = M3c_sliding_window[[5]]$Dataset, metric = "AIC") 

## Site-specific sea surface temperature
plotall(datasetrand = M3c_random_window[[1]],
        dataset = M3c_sliding_window[[1]]$Dataset, 
        bestmodel = M3c_sliding_window[[1]]$BestModel,
        bestmodeldata = M3c_sliding_window[[1]]$BestModelData)

## Site-specific precipitation  
plotall(datasetrand = M3c_random_window[[2]],
        dataset = M3c_sliding_window[[2]]$Dataset, 
        bestmodel = M3c_sliding_window[[2]]$BestModel,
        bestmodeldata = M3c_sliding_window[[2]]$BestModelData)

## Regional Pacific Decadal Oscillation
plotall(datasetrand = M3c_random_window[[3]],
        dataset = M3c_sliding_window[[3]]$Dataset, 
        bestmodel = M3c_sliding_window[[3]]$BestModel,
        bestmodeldata = M3c_sliding_window[[3]]$BestModelData)

## Regional El Niño Southern Oscillation
plotall(datasetrand = M3c_random_window[[4]],
        dataset = M3c_sliding_window[[4]]$Dataset, 
        bestmodel = M3c_sliding_window[[4]]$BestModel,
        bestmodeldata = M3c_sliding_window[[4]]$BestModelData)

## Regional Dipole Mode Index
plotall(datasetrand = M3c_random_window[[5]],
        dataset = M3c_sliding_window[[5]]$Dataset, 
        bestmodel = M3c_sliding_window[[5]]$BestModel,
        bestmodeldata = M3c_sliding_window[[5]]$BestModelData)

### Sliding window for rainfall data #### 
M3c_plot <- as.data.frame(M3c_sliding_window[[2]]$Dataset)

plot12 <- M3c_plot %>% 
  ggplot(aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+
  geom_tile(aes(fill=deltaAICc))+
  geom_tile(data = M3c_plot[c(1),], fill = NA, color = "white", linewidth = 1.5) +
  scale_fill_viridis(option = 'magma')+
  scale_x_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec\n(t-1)','Nov\n(t-1)','Oct\n(t-1)','Sep\n(t-1)','Aug\n(t-1)','Jul\n(t-1)','Jun\n(t-1)','May\n(t-1)','Apr\n(t-1)','Mar\n(t-1)','Feb\n(t-1)','Jan\n(t-1)','Dec\n(t-2)','Nov\n(t-2)','Oct\n(t-2)'), breaks=c(0:23), expand=c(0,0))+
  scale_y_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec (t-1)','Nov (t-1)','Oct (t-1)','Sep (t-1)','Aug (t-1)','Jul (t-1)','Jun (t-1)','May (t-1)','Apr (t-1)','Mar (t-1)','Feb (t-1)','Jan (t-1)','Dec (t-2)','Nov (t-2)','Oct (t-2)'), breaks=c(0:23), expand=c(0,0))+
  labs(x="Window Close", y="Window Open")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot12

### Distribution of AICc values plotted in a histogram
M3c_hist <- as.data.frame(M3c_random_window[[2]])

plot13 <- M3c_hist %>%
  ggplot(aes(x= deltaAICc))+
  geom_histogram(fill='grey20', bins = 50)+
  geom_hline(yintercept = 0, colour = 'grey80')+
  geom_vline(xintercept = -8.93, colour = 'red', linetype = 'dashed')+
  labs(x="Delta AICc", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot13

#################################################################################################################################
#################################################################################################################################

### Construct the optimal growth model for juvenile Lutjanus johnii from Northwestern Australia ####
## Our results indicate that juvenile growth was best explained by the rainfall patterns based on the p values.

## Site-specific rainfall
head(arrange(M3c_sliding_window[[2]]$Dataset), 1)
plotdelta(dataset = M3c_sliding_window[[2]]$Dataset)

## Check the minimum Fish Year for samples from Northwestern Australia
## Ensure the annual growth band is explained by rainfall conditions from the prior two years (Oct 1981 - Sept 1983).
min(Ljohnii_Australia_juvenile$FishYear)

## Group 1
env_Australia_juvenile1 <- env_data %>% 
  filter(Site == "Kimberley") %>% 
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>%
  subset(select = c(ymd, Site, Year, Month, Day, precip, s.precip, group, group1))

env_Australia_juvenile1 <- env_Australia_juvenile1 %>% 
  group_by(Site) %>% 
  mutate(group1 = consecutive_id(2 * floor(group1 / 2))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>% 
  group_by(Site, group1) %>% 
  mutate(group2 = row_number() - 1)

env_Australia_juvenile_name1 <- env_Australia_juvenile1 %>% 
  group_by(Site, group1) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/"))

env_Australia_juvenile1 <- env_Australia_juvenile1 %>% 
  filter(group2 == "7" | group2 == "6") %>%
  group_by(Site, group1) %>%
  summarise(mean_s.precip = mean(s.precip), mean_precip = mean(precip)) 

env_Australia_juvenile1 <- merge(env_Australia_juvenile1, env_Australia_juvenile_name1)
env_Australia_juvenile1 <- env_Australia_juvenile1 %>%
  filter(!group1 == "20") 
env_Australia_juvenile1

## Group 2 
env_Australia_juvenile2 <- env_data %>% 
  filter(Site == "Kimberley") %>% 
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>%
  subset(select = c(ymd, Site, Year, Month, Day, precip, s.precip, group, group1))

env_Australia_juvenile2 <- env_Australia_juvenile2 %>% 
  group_by(Site) %>% 
  mutate(group2 = consecutive_id(2 * floor((group1 + 1) / 2))) %>% 
  mutate(group2 = as.numeric(as.factor(group2)) - 1) %>%
  group_by(Site, group2) %>% 
  mutate(group3 = row_number() - 1)
# Note the difference in code between group 1 and 2.

env_Australia_juvenile_name2 <- env_Australia_juvenile2 %>% 
  group_by(Site, group2) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/")) 

env_Australia_juvenile2 <- env_Australia_juvenile2 %>% 
  filter(group3 == "7" | group3 == "6") %>%
  group_by(Site, group2) %>%
  summarise(mean_s.precip = mean(s.precip), mean_precip = mean(precip)) 

env_Australia_juvenile2 <- merge(env_Australia_juvenile2, env_Australia_juvenile_name2)
env_Australia_juvenile2 <- env_Australia_juvenile2 %>%
  filter(!group2 == "0") %>%
  rename(group1 = group2)
env_Australia_juvenile2

env_Australia_juvenile <- rbind(env_Australia_juvenile1, env_Australia_juvenile2)
env_Australia_juvenile <- env_Australia_juvenile %>%
  arrange(s_Year) %>%
  subset(select = -c(group1))
env_Australia_juvenile <- env_Australia_juvenile %>% 
  tidyr::separate(col="s_Year", sep=c(4), into=c("FishYear", "delim1"), remove=FALSE) %>%
  subset(select = -c(delim1))
env_Australia_juvenile

### Construct the optimal growth model for Lutjanus johnii ###
Ljohnii_Australia_juvenile_final <- merge(Ljohnii_Australia_juvenile, env_Australia_juvenile, by = c('Site','FishYear'))

M3c_final <- lmer(log(Width) ~ sAge + sAAC + 
                    mean_s.precip + 
                    (sAge|SampleID) + (1|Cohort), 
                  Ljohnii_Australia_juvenile_final, REML= T)

summary(M3c_final)
anova(M3c_final)
rsquared(M3c_final)

M3c_conf <- sim(M3c_final, 5000)
M3c_conf <- apply(M3c_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3c_conf

p <- predictorEffect("mean_s.precip", (M3c_final)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto")) 

### Model fitted with local and regional environmental variables ###
## Linear site-specific precipitation
M3c_precip <- as.data.frame (Effect (c('mean_s.precip'), M3c_final, xlevels= list (mean_s.precip = seq(-1.173, 1.78, by= 0.01)) )) 
M3c_mean_precip <- mean(Ljohnii_Australia_juvenile_final$mean_precip) ; M3c_mean_precip
M3c_sd_precip <- sd(Ljohnii_Australia_juvenile_final$mean_precip) ; M3c_sd_precip
M3c_precip$precip <- (M3c_precip$mean_s.precip * M3c_sd_precip + M3c_mean_precip)

M3c_precip$transfit <- exp(M3c_precip$fit)
M3c_precip$transupper <- exp(M3c_precip$upper)
M3c_precip$translower <- exp(M3c_precip$lower)

min(M3c_precip$precip) ; max(M3c_precip$precip)
min(Ljohnii_Australia_juvenile_final$mean_precip) ; max(Ljohnii_Australia_juvenile_final$mean_precip)

## Calculate the effect of precipitation on Lutjanus johnii growth
precip_min <- M3c_precip %>% filter(precip == min(precip)) %>% dplyr::select(transfit) 
precip_max <- M3c_precip %>% filter(precip == max(precip)) %>% dplyr::select(transfit)
precip_unit <- (M3c_precip %>% filter(precip == max(precip)) %>% dplyr::select(precip)) - (M3c_precip %>% filter(precip == min(precip)) %>% dplyr::select(precip))
precip_rate <- ((precip_max - precip_min) / (precip_unit)) * 100
precip_rate

### Plot the effect of site-specific precipitation on Lutjanus johnii growth 
plot14 <- M3c_precip %>%
  ggplot(aes(y = transfit, x = precip)) +
  geom_line(linewidth=0.75) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.3, alpha = 0.3, colour = 0) + 
  scale_x_continuous(limits = c(0, 80), breaks = c(0.0, 20.0, 40.0, 60.0, 80.0), labels = c('0.0', '20.0', '40.0', '60.0', '80.0')) +
  scale_y_continuous(limits = c(0.26, 0.44), breaks = c(0.26, 0.30, 0.34, 0.38, 0.42), labels = c('0.26', '0.30', '0.34', '0.38', '0.42')) +
  labs(y = 'Annual otolith growth (mm)', x = 'Mean site rainfall (mm)') +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        legend.position.inside = c(0.8,0.125),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot14

ggpubr::ggarrange(plot12, plot14, labels = c("A", "B"))

#################################################################################################################################
#################################################################################################################################

### Prepare the data for the exploratory sliding window analysis
Ljohnii_Asia <- Ljohnii_Asia %>%
  dplyr::rename(capture_year = Year, capture_month = Month, capture_day = Day)

## Dummy variable for the alignment of environmental variables
Ljohnii_Asia$Month <- "9" 
Ljohnii_Asia$Month <- as.numeric(Ljohnii_Asia$Month)
Ljohnii_Asia <- Ljohnii_Asia %>% mutate(ymd = make_date(FishYear, Month)) 
Ljohnii_Asia 

#### Base model ####
M4a <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Ljohnii_Asia, REML= F)

## Applying the sliding window to the environmental data ####
M4a_sliding_window <- slidingwin(baseline = M4a, 
                                 xvar = list(sst = env_data$s.temp,
                                             rain = env_data$s.precip,
                                             pdo = env_data$s.pdo,
                                             enso = env_data$s.enso,
                                             dmi = env_data$s.dmi),  
                                 type = "absolute",                                 # "absolute" or "relative", relative climate window=number of days before each biological record is measured; absolute climate window=number of days before a set point (refday) in time.
                                 range = c(23, 0),                                  # the furthest and closest number of time intervals (set by cinterval)
                                 stat = c("mean"),                                  # aggregate statistics used
                                 func = c("lin"),                                   # functions used to fit the climate variable. Can be linear ("lin"), quadratic ("quad"), cubic ("cub"), inverse ("inv") or log ("log")
                                 refday = c(1, 9),                                 # day and month respectively of the year from which the absolute window analysis will start
                                 cinterval = "month",                               # resolution at which climate window analysis will be conducted. 
                                 cdate = env_data$ymd,                              # climate date variable (dd/mm/yyyy)
                                 bdate = Ljohnii_Asia$ymd,
                                 cohort = Ljohnii_Asia$FishYear,       # cohort is required because our environmental data spans across two years (see Advanced Vignette for more information)
                                 spatial = list(Ljohnii_Asia$Site, env_data$Site))

### Sliding window data ####
## The 'combos' function provides an overview of the best fitted climate window, the delta AICc value, the start and end window date, and slope of the best window. 
M4a_sliding_window
M4a_sliding_window$combos 

### Applying the random window approach 
## The random window is repeated 1000 times to estimate the probability that a given result represents a false positive. 
## It randomises a given dataset (i.e., removes any climate signal) and conducts a sliding window analysis to extract a value of delta AIC. 
M4a_random_window <- randwin(baseline = M4a, 
                             repeats = 1000, # Randomisation of 1000 times 
                             xvar = list(sst = env_data$s.temp,
                                         rain = env_data$s.precip,
                                         pdo = env_data$s.pdo,
                                         enso = env_data$s.enso,
                                         dmi = env_data$s.dmi), 
                             type = "absolute",                       
                             range = c(23, 0),                      
                             stat = c("mean"),                          
                             func = c("lin"),                         
                             refday=c(1, 9),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Ljohnii_Asia$ymd,
                             cohort = Ljohnii_Asia$FishYear,
                             spatial = list(Ljohnii_Asia$Site, env_data$Site),
                             window= "sliding")

## There are warning messages when running the random window. 
## Fixed effect model matrix is rank deficient so dropping 1 column / coefficient. 

## As the reference date does not matter in this scenario, the absolute sliding window indicates that annual growth is dependent on the average annual temperature at depth conditions across the Gascoyne region over the Fish Year. 

## Calculate the p value to establish significance 
pvalue(datasetrand = M4a_random_window[[1]], dataset = M4a_sliding_window[[1]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M4a_random_window[[2]], dataset = M4a_sliding_window[[2]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M4a_random_window[[3]], dataset = M4a_sliding_window[[3]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M4a_random_window[[4]], dataset = M4a_sliding_window[[4]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M4a_random_window[[5]], dataset = M4a_sliding_window[[5]]$Dataset, metric = "AIC") 

### Exploratory plots ####
## Note based on the warning that Pc will be overly conservative when the sample size is greater than 47. 
## Site-specific sea surface temperature
plotall(datasetrand = M4a_random_window[[1]],
        dataset = M4a_sliding_window[[1]]$Dataset, 
        bestmodel = M4a_sliding_window[[1]]$BestModel,
        bestmodeldata = M4a_sliding_window[[1]]$BestModelData)

## Site-specific precipitation  
plotall(datasetrand = M4a_random_window[[2]],
        dataset = M4a_sliding_window[[2]]$Dataset, 
        bestmodel = M4a_sliding_window[[2]]$BestModel,
        bestmodeldata = M4a_sliding_window[[2]]$BestModelData)

## Regional Pacific Decadal Oscillation
plotall(datasetrand = M4a_random_window[[3]],
        dataset = M4a_sliding_window[[3]]$Dataset, 
        bestmodel = M4a_sliding_window[[3]]$BestModel,
        bestmodeldata = M4a_sliding_window[[3]]$BestModelData)

## Regional El Niño Southern Oscillation
plotall(datasetrand = M4a_random_window[[4]],
        dataset = M4a_sliding_window[[4]]$Dataset, 
        bestmodel = M4a_sliding_window[[4]]$BestModel,
        bestmodeldata = M4a_sliding_window[[4]]$BestModelData)

## Regional Dipole Mode Index
plotall(datasetrand = M4a_random_window[[5]],
        dataset = M4a_sliding_window[[5]]$Dataset, 
        bestmodel = M4a_sliding_window[[5]]$BestModel,
        bestmodeldata = M4a_sliding_window[[5]]$BestModelData)

### Sliding window for temperature data #### 
M4a_plot <- as.data.frame(M4a_sliding_window[[1]]$Dataset)

plot15 <- M4a_plot %>% 
  ggplot(aes(x=WindowClose, y=WindowOpen, fill=deltaAICc))+
  geom_tile(aes(fill=deltaAICc))+
  geom_tile(data = M4a_plot[c(1),], fill = NA, color = "white", linewidth = 1.5) +
  scale_fill_viridis(option = 'magma')+
  scale_x_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec\n(t-1)','Nov\n(t-1)','Oct\n(t-1)','Sep\n(t-1)','Aug\n(t-1)','Jul\n(t-1)','Jun\n(t-1)','May\n(t-1)','Apr\n(t-1)','Mar\n(t-1)','Feb\n(t-1)','Jan\n(t-1)','Dec\n(t-2)','Nov\n(t-2)','Oct\n(t-2)'), breaks=c(0:23), expand=c(0,0))+
  scale_y_continuous(labels=c('Sep','Aug','Jul','Jun','May','Apr','Mar','Feb','Jan','Dec (t-1)','Nov (t-1)','Oct (t-1)','Sep (t-1)','Aug (t-1)','Jul (t-1)','Jun (t-1)','May (t-1)','Apr (t-1)','Mar (t-1)','Feb (t-1)','Jan (t-1)','Dec (t-2)','Nov (t-2)','Oct (t-2)'), breaks=c(0:23), expand=c(0,0))+
  labs(x="Window Close", y="Window Open")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=10),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot15

### Distribution of AICc values plotted in a histogram
M4a_hist <- as.data.frame(M4a_random_window[[1]])

plot16 <- M4a_hist %>%
  ggplot(aes(x= deltaAICc))+
  geom_histogram(fill='grey20', bins = 50)+
  geom_hline(yintercept = 0, colour = 'grey80')+
  geom_vline(xintercept = -17.09, colour = 'red', linetype = 'dashed')+
  labs(x="Delta AICc", y="Count")+
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot16

#################################################################################################################################
#################################################################################################################################

### Construct the optimal growth model for Lutjanus johnii from Southeast Asia ####
## Our results indicate that growth was best explained by SST based on the p values.

## Site-specific sea surface temperature
head(arrange(M4a_sliding_window[[1]]$Dataset), 1)
plotdelta(dataset = M4a_sliding_window[[1]]$Dataset)

## Check the minimum Fish Year for samples from Southeast Asia
## Ensure the annual growth band is explained by SST conditions from the prior two years (Oct 1993 - Sept 1995).
min(Ljohnii_Asia$FishYear)

## Group 1
env_Asia1 <- env_data %>% 
  filter(!Site == "Kimberley") %>% 
  filter(!ymd < "1992-10-01" & !ymd > "2022-09-01") %>%
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) 

env_Asia1 <- env_Asia1 %>% 
  group_by(Site) %>% 
  mutate(group1 = consecutive_id(2 * floor(group1 / 2))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) %>% 
  group_by(Site, group1) %>% mutate(group2 = row_number() - 1)

env_Asia_name1 <- env_Asia1 %>% 
  group_by(Site, group1) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/"))

env_Asia1 <- env_Asia1 %>% 
  filter(!group2 == "0" | !group2 == "1") %>%
  group_by(Site, group2) %>%
  summarise(mean_s.temp = mean(s.temp), mean_temp = mean(temp)) %>%
  rename(group1 = group2)

env_Asia1 <- merge(env_Asia1, env_Asia_name1)
env_Asia1

## Group 2 
env_Asia2 <- env_data %>% 
  filter(!Site == "Kimberley") %>% 
  filter(!ymd < "1992-10-01" & !ymd > "2022-09-01") %>%
  arrange(desc(ymd)) %>% 
  group_by(Site) %>% 
  mutate(group = (1:n()) - 1) %>%
  mutate(group1 = consecutive_id(12 * floor(group / 12))) %>% 
  mutate(group1 = as.numeric(as.factor(group1)) - 1) 

env_Asia2 <- env_Asia2 %>% 
  group_by(Site) %>% 
  mutate(group2 = consecutive_id(2 * floor((group1 + 1) / 2))) %>% 
  mutate(group2 = as.numeric(as.factor(group2)) - 1) %>%
  group_by(Site, group2) %>% mutate(group3 = row_number() - 1)
# Note the difference in code between group 1 and 2.

env_Asia_name2 <- env_Asia2 %>% 
  group_by(Site, group2) %>% 
  summarise(s_Year = str_c(unique(Year), collapse = "/")) 

env_Asia2 <- env_Asia2 %>% 
  filter(!group2 == "0" | !group2 == "1") %>%
  group_by(Site, group2) %>%
  summarise(mean_s.temp = mean(s.temp), mean_temp = mean(temp)) 
  
env_Asia2 <- merge(env_Asia2, env_Asia_name2)
env_Asia2 <- env_Asia2 %>%
  filter(!group2 == "0" & !group2 == "15") %>%
  rename(group1 = group2)
env_Asia2

env_Asia <- rbind(env_Asia1, env_Asia2)
env_Asia <- env_Asia %>%
  arrange(s_Year) %>%
  subset(select = -c(group1))
env_Asia <- env_Asia %>% 
  tidyr::separate(col="s_Year", sep=c(4), into=c("FishYear", "delim1"), remove=FALSE) %>%
  subset(select = -c(delim1))
env_Asia

### Construct the optimal growth model for Lutjanus johnii ###
Ljohnii_Asia_final <- merge(Ljohnii_Asia, env_Asia, by = c('Site','FishYear'))

M4a_final <- lmer(log(Width) ~ sAge + sAAC + 
                    mean_s.temp + 
                    (sAge|SampleID) + (1|FishYear), 
                  Ljohnii_Asia_final, REML= T)

summary(M4a_final)
anova(M4a_final)
rsquared(M4a_final)

M4a_conf <- sim(M4a_final, 5000)
M4a_conf <- apply(M4a_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M4a_conf

p <- predictorEffect("mean_s.temp", (M4a_final)) ; plot(p, lines=list(multiline=TRUE), confint=list(style="auto")) 

### Model fitted with local and regional environmental variables ###
## Linear site-specific temperature 
M4a_temp <- as.data.frame (Effect (c('mean_s.temp'), M4a_final, xlevels= list (mean_s.temp = seq(-3.0, 1.7, by= 0.01)) )) 
M4a_mean_temp <- mean(Ljohnii_Asia_final$mean_temp) ; M4a_mean_temp
M4a_sd_temp <- sd(Ljohnii_Asia_final$mean_temp) ; M4a_sd_temp
M4a_temp$temp <- (M4a_temp$mean_s.temp * M4a_sd_temp + M4a_mean_temp)

M4a_temp$transfit <- exp(M4a_temp$fit)
M4a_temp$transupper <- exp(M4a_temp$upper)
M4a_temp$translower <- exp(M4a_temp$lower)

min(M4a_temp$temp) ; max(M4a_temp$temp)
min(Ljohnii_Asia_final$mean_temp) ; max(Ljohnii_Asia_final$mean_temp)

## Calculate the effect of temperature on Lutjanus johnii growth
temp_min <- M4a_temp %>% filter(temp == min(temp)) %>% dplyr::select(transfit) 
temp_max <- M4a_temp %>% filter(temp == max(temp)) %>% dplyr::select(transfit)
temp_unit <- (M4a_temp %>% filter(temp == max(temp)) %>% dplyr::select(temp)) - (M4a_temp %>% filter(temp == min(temp)) %>% dplyr::select(temp))
temp_rate <- ((temp_max - temp_min) / (temp_unit)) * 100
temp_rate

### Plot the effect of temperature on Lutjanus johnii growth 
plot17 <- M4a_temp %>%
  ggplot(aes(y = transfit, x = temp)) +
  geom_line(linewidth=0.75) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.3, alpha = 0.3, colour = 0) + 
  scale_x_continuous(limits = c(27.4, 30.5), breaks = c(27.5, 28.0, 28.5, 29.0, 29.5, 30.0, 30.5), labels = c('27.5', '28.0', '28.5','29.0', '29.5', '30.0', '30.5')) +
  scale_y_continuous(limits = c(0.090, 0.235), breaks = c(0.100, 0.125, 0.150, 0.175, 0.200, 0.225), labels = c('0.100', '0.125', '0.150', '0.175', '0.200', '0.225')) +
  labs(y = 'Annual otolith growth (mm)', x = 'Mean site temperature (Celsius)') +
  theme(plot.title = element_text(hjust = 0.5, size=22),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 14),
        axis.line = element_line(colour = "grey15", linewidth=0.3), 
        panel.border = element_blank(),
        panel.background = element_rect(colour="grey50", linewidth=0.3, fill="white"),
        panel.grid.major = element_line(colour = "grey85", linewidth=0.15),
        panel.grid.minor = element_line(colour = "grey85", linewidth=0.15),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        legend.position.inside = c(0.8,0.125),
        legend.background = element_rect(fill = "white", color = "black"),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot17

ggpubr::ggarrange(plot15, plot17, labels = c("A", "B"))

ggpubr::ggarrange(plot9, plot10, plot12, plot13, plot15, plot16, labels = c("A", "D", "B", "E", "C", "F"), ncol = 2, nrow = 3)

ggpubr::ggarrange(plot11, plot14, plot17, labels = c("A", "B", "C"), ncol = 1, nrow = 3)

#################################################################################################################################
#################################################################################################################################

#### Study region ####
## New update as of 12 November 2023, Stadia Maps requires an API Key.
## I will be using Google Maps instead. 
remotes::install_github("dkahle/ggmap")
library(ggmap)
register_google(key = "AIzaSyBBfJPr7DR7GYl6FxeslPLBqQE7cia55j4", write = TRUE)

lutjanus_johnii_study_region <- get_map(location = c(110, 14.25), maptype='satellite', source="google", zoom=3)

eez_shp <- readOGR("../../Climate/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp") %>% spTransform("+proj=longlat +ellps=WGS84")
eez_shp <- ggplot2::fortify(eez_shp)

ggmap(lutjanus_johnii_study_region) +
  #geom_path(data=eez_shp, aes(long, lat, group=group), linewidth= 0.075)+ # Exclusive economic zones
  xlab("Longitude") + ylab("Latitude")+
  scale_x_continuous(limits = c(80, 140), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-25, 10), expand = c(0, 0)) +
  
  annotate('point', y= 3.5, x= 100.0, color='firebrick3', size=20, alpha=0.4)+
  annotate('point', y= 2.2, x= 107.0, size= 25, color='firebrick3', alpha= 0.5) + 
  #annotate('point', y= -18.5, x= 117.0, size= 25, color='firebrick3', alpha= 0.5) + 
  annotate('point', y= -14.0, x= 124.5, size= 25, color='firebrick3', alpha= 0.5) + 
  
  annotate('text', y=  3.5, x= 100.0, color='grey98', size= 3, label='Malacca\nStraits', fontface = 2) +
  annotate('text', y=  2.2, x= 107.0, color='grey98', size= 3, label='Riau\nArchipelago', fontface= 2) +
  #annotate('text', y= -18.5, x= 117.0, color='grey98', size= 3, label='Pilbara', fontface= 2) +
  annotate('text', y= -14.0, x= 124.5, color='grey98', size= 3, label='Kimberley', fontface= 2) +
  
  theme(plot.title = element_text(hjust=0.5, size=18), 
        axis.title=element_text(size=18),
        axis.text = element_text(size=16),
        legend.position='none')

#################################################################################################################################
#################################################################################################################################

## Calculation of IAPE and ACV for Lutjanus johnii
SEA_age <- readxl::read_xlsx("Julio_Snapper_Readings_SEA.xlsx", 1, col_names = F)
SEA_age <- SEA_age %>%
  subset(select = -c(2,3,5,7,8:10)) %>%
  dplyr::rename("SampleID" = "...1", "Reader1" = "...4", "Reader2" = "...6") 
SEA_age = SEA_age[-1,]
SEA_age$Reader1 <- as.numeric(SEA_age$Reader1)
SEA_age$Reader2 <- as.numeric(SEA_age$Reader2)
SEA_age$diff <- SEA_age$Reader1 - SEA_age$Reader2

WA_age <- readxl::read_xlsx("Corey_Snapper_Readings_WA.xlsx", 1, col_names = F)
WA_age <- WA_age %>%
  subset(select = -c(3,4,6:7)) %>%
  dplyr::rename("SampleID" = "...1", "Reader1" = "...2","Reader2" = "...5" ) 
WA_age = WA_age[-1,]
WA_age$Reader1 <- as.numeric(WA_age$Reader1)
WA_age$Reader2 <- as.numeric(WA_age$Reader2)
WA_age$diff <- WA_age$Reader1 - WA_age$Reader2
WA_age <- as.data.frame(WA_age)

SEA_precision <- agePrecision(~Reader1+Reader2, data=SEA_age)
summary(SEA_precision, what="precision")
SEA_precision

WA_precision <- agePrecision(~Reader1+Reader2, data=WA_age)
summary(WA_precision, what="precision")
WA_precision

#################################################################################################################################
#################################################################################################################################

#save.image("lutjanus_johnii.RData")

#################################################################################################################################
#################################################################################################################################


