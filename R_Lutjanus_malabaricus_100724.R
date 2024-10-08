#### Chapter 2: Growth dynamics of Tropical Marine Snappers from the Indo-Pacific 
## Author: Clement Ng 

#### Lutjanus malabaricus ####

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
library(lattice)
library(viridis)
library(rgdal)
library(sf)
library(terra)

cbPalette <- c("#D55E00", "#E69F00", "#009E73", "#56B4E9", "#0072B2", "#F0E442", "#CC79A7", "#999999") 

c. <- function (x) scale(x, scale = FALSE) 

#load("lutjanus_malabaricus.RData") 

#################################################################################################################################
#################################################################################################################################

#### Data coding ####
## Lutjanus malabaricus from Southeast Asia 
SEA_lutjanus_malabaricus_txt_files <- list.files(path = ".", recursive = TRUE,
                                                 pattern = "\\.txt$", 
                                                 full.names = TRUE)
SEA_lutjanus_malabaricus_increment_data <- rbindlist(sapply(SEA_lutjanus_malabaricus_txt_files, fread, simplify = FALSE),
                                                     use.names = TRUE, fill = TRUE, idcol = "TxtFileName")

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  subset(select = -c(ypos, xpos)) %>%
  dplyr::rename("CalYear" = "IncNum v1.3c") %>%
  rename("Width" = "Thickness (mm)")

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  tidyr::separate(col="TxtFileName", sep=c(2, 7), into=c("delim1", "FishID", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim1, delim2))

SEA_lutjanus_malabaricus_increment_data$Width[SEA_lutjanus_malabaricus_increment_data$Width == 0.000] <- NA
SEA_lutjanus_malabaricus_increment_data <- na.omit(SEA_lutjanus_malabaricus_increment_data) # remove points on the edge

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>%
  arrange(CalYear, .by_group = TRUE) %>% 
  mutate(Age = 1:n()+1) %>% # Add one point because the first ring was not included in the growth measurements
  ungroup()

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Age == max(Age) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>% # Extract the edge as a new column in the data frame
  ungroup()

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(EdgeRatio = case_when(Age == (max(Age) - 1) ~ Width)) %>% 
  fill(EdgeRatio, .direction = "downup") %>% # Establish the ratio between the edge and the thickness of the previous ring to identify if rings have been missed
  mutate(EdgeRatio = Edge / EdgeRatio) %>%
  mutate(EdgeType = case_when(EdgeRatio < 0.2 ~ "Narrow", 
                              EdgeRatio >= 0.2 & EdgeRatio <= 0.7 ~ "Intermediate",
                              EdgeRatio > 0.7 ~ "Wide")) %>% 
  ungroup()

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge_Delete = case_when(Age == (max(Age)) ~ "Edge")) 

SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data[- grep("Edge", SEA_lutjanus_malabaricus_increment_data$Edge_Delete),]
SEA_lutjanus_malabaricus_increment_data <- subset(SEA_lutjanus_malabaricus_increment_data, select = -c(Edge_Delete))
SEA_lutjanus_malabaricus_increment_data <- SEA_lutjanus_malabaricus_increment_data %>% rename("SampleID" = "FishID") 

#### Inputting fish data into the increment data frame ####
SEA_lutjanus_malabaricus_fish_data <- read.csv("lutjanus_malabaricus_fish_data.csv")
colnames(SEA_lutjanus_malabaricus_fish_data)
SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  rename(Year = Capture.Year, Month = Capture.Month, Day = Capture.Day, Method = Harvest.Method, 
         TL = Total.Length..mm., FL = Fork.Length..mm., SL = Standard.Length..mm., Weight = Wet.weight..g.,
         LOtoWt = L.Otolith.Weight..g., ROtoWt = R.Otolith.Weight..g.) 

unique(SEA_lutjanus_malabaricus_fish_data$Harvest.Source) 
unique(SEA_lutjanus_malabaricus_fish_data$State.of.Left.Otolith) 
unique(SEA_lutjanus_malabaricus_fish_data$State.of.Right.Otolith) 

SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>% 
  filter(Harvest.Source == "Commercial") %>%
  mutate(LOtoWt = case_when(State.of.Left.Otolith == "Chipped" | State.of.Left.Otolith == "Broken" ~ NA, TRUE ~ LOtoWt)) %>% # remove otolith weights if they are chipped
  mutate(ROtoWt = case_when(State.of.Right.Otolith == "Chipped" | State.of.Right.Otolith == "Chipped " | State.of.Right.Otolith == "Broken" ~ NA, TRUE ~ ROtoWt)) %>% 
  subset(select = c(1, 3:5, 12, 17:23, 29:30)) %>%
  rename(Stage = `Sex.Stage`)
SEA_lutjanus_malabaricus_fish_data$Stage <- as.character(SEA_lutjanus_malabaricus_fish_data$Stage)

unique(SEA_lutjanus_malabaricus_fish_data$Location)
SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  mutate(Region = case_when(Location == "Eastern Makassar" ~ "Sulawesi", 
                            
                            Location == "Melacca Straits" ~ "Malacca Straits", 
                            Location == "Malacca" ~ "Malacca Straits", 
                            Location == "Medan" ~ "Malacca Straits", 
                            
                            Location == "Anambas" ~ "Riau Archipelago", 
                            Location == "Batam" ~ "Riau Archipelago", 
                            Location == "Kijang" ~ "Riau Archipelago", 
                            Location == "South-Eastern Sumatra" ~ "Riau Archipelago", 
                            
                            Location == "Unknown" ~ "",
                            TRUE ~ Location))
unique(SEA_lutjanus_malabaricus_fish_data$Region)

## Reclassify the names into grouping categories
SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  mutate(Location = case_when(Location == "Malacca" ~ "Malacca Straits", 
                              Location == "Melacca Straits" ~ "Malacca Straits", 
                              Location == "Medan" ~ "Malacca Straits", 
                              TRUE ~ Location))

SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  mutate(Location = case_when(Location == "Anambas" ~ "Riau Archipelago", 
                              Location == "Batam" ~ "Riau Archipelago", 
                              Location == "Kijang" ~ "Riau Archipelago", 
                              Location == "South-Eastern Sumatra" ~ "Riau Archipelago", 
                              TRUE ~ Location))

SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  mutate(Location = case_when(Location == "Unknown" ~ "", 
                              TRUE ~ Location))

SEA_lutjanus_malabaricus_fish_data <- SEA_lutjanus_malabaricus_fish_data %>%
  mutate(Location = case_when(Location == "" ~ "", 
                              TRUE ~ Location))
unique(SEA_lutjanus_malabaricus_fish_data$Location)

#### Merge to form a dataset ####
SEA_lutjanus_malabaricus_full_data <- merge(SEA_lutjanus_malabaricus_increment_data, SEA_lutjanus_malabaricus_fish_data, by = "SampleID")
SEA_lutjanus_malabaricus_full_data <- SEA_lutjanus_malabaricus_full_data %>% rename(Site = Location)
SEA_lutjanus_malabaricus_full_data$birth_month = 10 # Lutjanus malabaricus annuli were completed in late Sept based on Cappo et al., (2000)
SEA_lutjanus_malabaricus_full_data <- as.data.frame(SEA_lutjanus_malabaricus_full_data) 

SEA_lutjanus_malabaricus_full_data$TL <- as.integer(SEA_lutjanus_malabaricus_full_data$TL)
SEA_lutjanus_malabaricus_full_data$FL <- as.integer(SEA_lutjanus_malabaricus_full_data$FL)
SEA_lutjanus_malabaricus_full_data$SL <- as.integer(SEA_lutjanus_malabaricus_full_data$SL)

SEA_lutjanus_malabaricus_full_data$Year <- as.character(SEA_lutjanus_malabaricus_full_data$Year)
SEA_lutjanus_malabaricus_full_data$Month <- as.character(SEA_lutjanus_malabaricus_full_data$Month)
SEA_lutjanus_malabaricus_full_data$Day <- as.character(SEA_lutjanus_malabaricus_full_data$Day)
SEA_lutjanus_malabaricus_full_data$birth_month <- as.character(SEA_lutjanus_malabaricus_full_data$birth_month)

head(SEA_lutjanus_malabaricus_full_data, 3) ; str(SEA_lutjanus_malabaricus_full_data) # looks good
SEA_lutjanus_malabaricus_full_data %>% distinct(SampleID) # Here, we see that we have 108 unique individuals that were aged

#################################################################################################################################
#################################################################################################################################

### Input data into the model 
## Lutjanus malabaricus from Australia
WA_lutjanus_malabaricus_txt_files <- list.files(path = "../../2022 WA DPIRD/2022_DPIRD_Images/Lutjanus malabaricus/ImageJ", recursive = TRUE,
                                                pattern = "\\.txt$", 
                                                full.names = TRUE)
WA_lutjanus_malabaricus_increment_data <- rbindlist(sapply(WA_lutjanus_malabaricus_txt_files, fread, simplify = FALSE),
                                                    use.names = TRUE, fill = TRUE, idcol = "TxtFileName")

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  subset(select = -c(V1, ypos, xpos)) %>%
  rename("CalYear" = "IncNum v1.3c") %>%
  rename("Width" = "Thickness (µm)")

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  tidyr::separate(col="TxtFileName", sep=c(66), into=c("delim1", "FishID"), remove=TRUE) %>% 
  subset(select = -c(delim1)) %>%
  tidyr::separate(col="FishID", sep=c(-19), into=c("FishID", "delim2"), remove=TRUE) %>% 
  subset(select = -c(delim2)) 

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>%
  mutate(FishID = str_replace(FishID, "malabaricus_", "M")) %>% # Some are named in a different format, you can replace Lmalabaricus with LM
  mutate(FishID = str_replace(FishID, "_", "-")) 

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>%
  mutate(FishID = str_replace(FishID, "LM15", "S15")) # Samples that start with the number LM15xxx should be labelled S15xxx. 

WA_lutjanus_malabaricus_increment_data$Width[WA_lutjanus_malabaricus_increment_data$Width == 0.000] <- NA
WA_lutjanus_malabaricus_increment_data <- na.omit(WA_lutjanus_malabaricus_increment_data)

WA_lutjanus_malabaricus_increment_data$Width <- WA_lutjanus_malabaricus_increment_data$Width / 1000

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>%
  arrange(CalYear, .by_group = TRUE) %>% 
  mutate(Age = 1:n()+1) %>%
  ungroup()

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge = case_when(Age == max(Age) ~ Width)) %>%
  fill(Edge, .direction = "downup") %>%
  ungroup()

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(EdgeRatio = case_when(Age == (max(Age) - 1) ~ Width)) %>% 
  fill(EdgeRatio, .direction = "downup") %>%
  mutate(EdgeRatio = Edge / EdgeRatio) %>%
  mutate(EdgeType = case_when(EdgeRatio < 0.2 ~ "Narrow",
                              EdgeRatio >= 0.2 & EdgeRatio <= 0.7 ~ "Intermediate",
                              EdgeRatio > 0.7 ~ "Wide")) %>%
  ungroup()

WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% 
  group_by(FishID) %>% 
  mutate(Edge_Delete = case_when(Age == (max(Age)) ~ "Edge"))
WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data[- grep("Edge", WA_lutjanus_malabaricus_increment_data$Edge_Delete),]
WA_lutjanus_malabaricus_increment_data <- subset(WA_lutjanus_malabaricus_increment_data, select = -c(Edge_Delete))
WA_lutjanus_malabaricus_increment_data <- WA_lutjanus_malabaricus_increment_data %>% rename("SampleID" = "FishID") 

#### Inputting fish data into the increment data frame ####
## There are two data sheets containing the data -- sheet 1 contains sectioned otolith metadata and sheet 4 includes the unsectioned otolith metadata
## Sectioned otolith metadata
WA_lutjanus_malabaricus_fish_data <- readxl::read_xlsx("../../2022 WA DPIRD/DPIRD_Snapper_Biological_Data/DPIRD_Lutjanus_malabaricus.xlsx", 1) 
WA_lutjanus_malabaricus_fish_data <- WA_lutjanus_malabaricus_fish_data %>%
  rename(SampleID = `SlideNum`,Capture_Year = Year, Capture_Month = Month, Capture_Day = Day, Weight = TotalWt) %>%
  subset(select = -c(1, 3)) 

head(WA_lutjanus_malabaricus_fish_data, 3) ; str(WA_lutjanus_malabaricus_fish_data)
WA_lutjanus_malabaricus_fish_data$SampleID <- as.character(WA_lutjanus_malabaricus_fish_data$SampleID)
WA_lutjanus_malabaricus_fish_data$Capture_Year <- as.character(WA_lutjanus_malabaricus_fish_data$Capture_Year)
WA_lutjanus_malabaricus_fish_data$Capture_Month <- as.character(WA_lutjanus_malabaricus_fish_data$Capture_Month)
WA_lutjanus_malabaricus_fish_data$Capture_Day <- as.character(WA_lutjanus_malabaricus_fish_data$Capture_Day)
WA_lutjanus_malabaricus_fish_data$Stage <- as.character(WA_lutjanus_malabaricus_fish_data$Stage)

WA_lutjanus_malabaricus_fish_data <-  as.data.frame(WA_lutjanus_malabaricus_fish_data)
colnames(WA_lutjanus_malabaricus_fish_data)
setcolorder(WA_lutjanus_malabaricus_fish_data, c("SampleID", "Region", "Location", 
                                                 "Capture_Year", "Capture_Month", "Capture_Day",
                                                 "TL",  "FL", "SL", "Weight", "Sex"))

#### Merge to form a dataset ####
WA_lutjanus_malabaricus_full_data <- merge(WA_lutjanus_malabaricus_increment_data, WA_lutjanus_malabaricus_fish_data, by = "SampleID")
WA_lutjanus_malabaricus_full_data <- WA_lutjanus_malabaricus_full_data %>% rename(Site = Location)
sum(is.na(WA_lutjanus_malabaricus_full_data$Site)) 

WA_lutjanus_malabaricus_full_data$birth_month = 10 # Lutjanus malabaricus annuli were completed in late Sept based on Cappo et al., (2000)
WA_lutjanus_malabaricus_full_data <- as.data.frame(WA_lutjanus_malabaricus_full_data)

unique(WA_lutjanus_malabaricus_full_data$Region)
WA_lutjanus_malabaricus_full_data <- WA_lutjanus_malabaricus_full_data %>% 
  mutate(Region = case_when(Region == "PIlbara" ~ "Pilbara", 
                            TRUE ~ Region))

WA_lutjanus_malabaricus_full_data <- WA_lutjanus_malabaricus_full_data %>% 
  mutate(Region = case_when(Region == "Kimberley NDSF" ~ "Kimberley", 
                            TRUE ~ Region))

## Discussions with DPIRD suggested that I should combine the Regions due to the lack of genetic differences between populations. 
WA_lutjanus_malabaricus_full_data$Site <- WA_lutjanus_malabaricus_full_data$Region
unique(WA_lutjanus_malabaricus_full_data$Site)

WA_lutjanus_malabaricus_full_data$Region <- "NW Australia"
unique(WA_lutjanus_malabaricus_full_data$Region)

WA_lutjanus_malabaricus_full_data <- WA_lutjanus_malabaricus_full_data %>% rename(LOtoWt = LeftOtoWt)
WA_lutjanus_malabaricus_full_data$TL <- as.integer(WA_lutjanus_malabaricus_full_data$TL)
WA_lutjanus_malabaricus_full_data$FL <- as.integer(WA_lutjanus_malabaricus_full_data$FL)
WA_lutjanus_malabaricus_full_data$SL <- as.integer(WA_lutjanus_malabaricus_full_data$SL)
WA_lutjanus_malabaricus_full_data$Weight <- as.integer(WA_lutjanus_malabaricus_full_data$Weight)
WA_lutjanus_malabaricus_full_data <- WA_lutjanus_malabaricus_full_data %>% 
  rename(Year = Capture_Year, Month = Capture_Month, Day = Capture_Day) 
WA_lutjanus_malabaricus_full_data <- as.data.frame(WA_lutjanus_malabaricus_full_data)

head(WA_lutjanus_malabaricus_full_data, 3) ; str(WA_lutjanus_malabaricus_full_data) # looks good
WA_lutjanus_malabaricus_full_data %>% distinct(SampleID) # Here, we see that we have 340 unique individuals that were aged

## Additional check for missing metadata 
## The merge function is specified to be all = T, which allows for NAs in the data frame unlike in the earlier line of code, where it is all = F. 
checkMissing <- WA_lutjanus_malabaricus_full_data
checkMissing <- checkMissing %>%
  filter(is.na(Site) & !is.na(Width)) %>% # Filter Regions without any information but has increment Width data
  distinct(SampleID) 

checkMissing # 37 samples without data on the Location of capture
#data.table::fwrite(checkMissing, file = "missingMetadata.csv") 

#################################################################################################################################
#################################################################################################################################

#### Combine the Southeast Asia and Western Australia data sets ####
lutjanus_malabaricus_full_data <- merge(SEA_lutjanus_malabaricus_full_data, WA_lutjanus_malabaricus_full_data, all=TRUE)
lutjanus_malabaricus_full_data <- as.data.frame(lutjanus_malabaricus_full_data)
head(lutjanus_malabaricus_full_data, 3) ; str(lutjanus_malabaricus_full_data)

lutjanus_malabaricus_full_data %>% distinct(SampleID) 

lutjanus_malabaricus_full_data <- lutjanus_malabaricus_full_data %>% 
  group_by(SampleID) %>% 
  mutate(AAC = max(Age)) # fully formed rings

lutjanus_malabaricus_full_data <- lutjanus_malabaricus_full_data %>%
  mutate(Cohort = CalYear-Age)

lutjanus_malabaricus_full_data <- lutjanus_malabaricus_full_data %>% rename(FishYear = CalYear)

lutjanus_malabaricus_full_data <- lutjanus_malabaricus_full_data %>%
  mutate(MaturityStage = case_when(Age < 5 ~ "Immature", 
                                   Age >= 5 ~ "Mature"))

### Create grouping variables
lutjanus_malabaricus_full_data <- as.data.frame(lutjanus_malabaricus_full_data)
lutjanus_malabaricus_full_data$SampleID <- as.factor(lutjanus_malabaricus_full_data$SampleID)
lutjanus_malabaricus_full_data$Region <- as.factor(lutjanus_malabaricus_full_data$Region)
lutjanus_malabaricus_full_data$Site <- as.factor(lutjanus_malabaricus_full_data$Site)
lutjanus_malabaricus_full_data$Cohort <- as.factor(lutjanus_malabaricus_full_data$Cohort)

lutjanus_malabaricus_full_data$region_year <- as.factor(paste(lutjanus_malabaricus_full_data$Region, lutjanus_malabaricus_full_data$FishYear, sep = "_")) 
lutjanus_malabaricus_full_data$site_year <- as.factor(paste(lutjanus_malabaricus_full_data$Site, lutjanus_malabaricus_full_data$FishYear, sep = "_")) 
lutjanus_malabaricus_full_data$cohort_year <- as.factor(paste(lutjanus_malabaricus_full_data$Cohort, lutjanus_malabaricus_full_data$FishYear, sep = "_")) 
lutjanus_malabaricus_full_data$region_cohort <- as.factor(paste(lutjanus_malabaricus_full_data$Region, lutjanus_malabaricus_full_data$Cohort, sep = "_")) 
lutjanus_malabaricus_full_data$site_cohort <- as.factor(paste(lutjanus_malabaricus_full_data$Site, lutjanus_malabaricus_full_data$Cohort, sep = "_")) 

lutjanus_malabaricus_full_data$sAge <- c.(log(lutjanus_malabaricus_full_data$Age))
lutjanus_malabaricus_full_data$sAAC <- c.(log(lutjanus_malabaricus_full_data$AAC))

colnames(lutjanus_malabaricus_full_data)
setcolorder(lutjanus_malabaricus_full_data, c("SampleID", "FishYear", 
                                              "Year", "Month", "Day", "birth_month", 
                                              "Method", "Region", "region_year", "Site", "site_year", 
                                              "Cohort", "cohort_year", "region_cohort", "site_cohort",
                                              "TL", "FL", "SL", "Weight", "Sex", "Stage","MaturityStage","GonadWt",
                                              "Fishery","Source","Depth",
                                              "LOtoWt", "ROtoWt", 
                                              "Width", "Age", "sAge", "AAC", "sAAC",
                                              "Edge", "EdgeRatio", "EdgeType"))

#################################################################################################################################
#################################################################################################################################

## Intrinsic variables
## Note that this data set is based on increments, so we need to split it based on distinct individuals based on the function, summarise(n_distinct(SampleID))

## Age at capture table and distribution
lutjanus_malabaricus_full_data %>%
  group_by(AAC) %>% summarise(n_distinct(SampleID)) %>% # summarise based on long format
  rename(n = `n_distinct(SampleID)`) %>%
  #spread(AAC, n) %>% # change to long format, if required
  print(n=Inf)

lutjanus_malabaricus_full_data %>% 
  group_by(AAC) %>% summarise(n=n_distinct(SampleID)) %>%
  ggplot(., aes(x=AAC, y=n))+ # histogram to highlight age distribution
  geom_bar(stat='identity',color='grey5',fill='grey45')+
  labs(y= 'Count')+ theme_bw()

## Region
lutjanus_malabaricus_full_data %>% 
  group_by(Region) %>% summarise(n=n_distinct(SampleID)) %>%
  print(n=Inf)

lutjanus_malabaricus_full_data %>% 
  filter(Region != "NA") %>%
  group_by(AAC, Region) %>% summarise(n=n_distinct(SampleID)) %>%
  ggplot(., aes(x=AAC, y=n))+ # histogram to highlight AAC distribution based on Region
  geom_bar(stat='identity',color='grey5',fill='grey45')+
  scale_y_continuous(limits=c(0, 80))+
  scale_x_continuous(limits=c(0, 30))+
  facet_wrap(.~Region, ncol = 2)+
  labs(y= 'Count', x= 'Age-at-Capture')+ theme_bw()
# We can see that majority of values are from NW Australia and there are some values from the Riau Archipelago. 

## Total length
## Replace with Standard and Fork lengths respectively. 
lutjanus_malabaricus_full_data %>% 
  group_by(TL) %>% summarise(n=n_distinct(SampleID)) %>%
  print(n=Inf)

lutjanus_malabaricus_full_data %>% 
  group_by(TL, Region) %>% summarise(n=n_distinct(SampleID)) %>% drop_na(TL) %>%
  ggplot(., aes(y=as.numeric(TL), x=Region, group=Region))+ # histogram to highlight age distribution
  geom_boxplot(color='grey5',fill='grey45')+ geom_jitter()+
  #facet_wrap(.~Region)+
  labs(y= 'Total Length (mm)')+ theme_bw()

#################################################################################################################################
#################################################################################################################################

#### Check missing values ####
colSums(is.na(lutjanus_malabaricus_full_data)) # note that this is repeated for all increments, not unique values
# There are a large number of missing values for the left and right otolith weights, as they were not recorded in Australia. 
# There are also missing methods of capture, FL and SL measurements. This indicates that we should ideally use TL measurements. 

## Create an additional table
checkTable <- lutjanus_malabaricus_full_data %>% 
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
checkTable <- checkTable %>% drop_na(TL, LOtoWt)
cor.test(checkTable$TL, checkTable$LOtoWt, method = "pearson") 

## Assess the relationship between left and right otoliths
checkTable <- checkTable %>% drop_na(LOtoWt, ROtoWt)
t.test(checkTable$LOtoWt, checkTable$ROtoWt) 

### Check for outliers in the response variable ####
## Examine the relationship between increment width and fish age for each individual
xyplot(Width ~ Age, group=SampleID, lutjanus_malabaricus_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus malabaricus)", ylab = "Width (mm)") # Overall pattern

xyplot(Width ~ Age | Region, group=SampleID, lutjanus_malabaricus_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus malabaricus)", ylab = "Width (mm)") # Regional pattern
# There may potentially be three outliers in the data frame, present in Riau Archipelago, and two data points in the Pilbara. 

xyplot(log(Width) ~ log(Age), group=SampleID, lutjanus_malabaricus_full_data, type=c('l','p'), 
       main = "Increment width by Age for individual fish (Lutjanus malabaricus)", ylab = "Width (mm)") # Log transformed

### Check for dependency 
lutjanus_malabaricus_full_data %>% 
  ggplot()+ 
  geom_point(aes(y=Width, x=Region, group=SampleID), position=position_dodge(0.8))+
  theme(text = element_text(size = 12), 
        legend.position = "none", 
        axis.text.x = element_text(size = 7, angle = 45, hjust = 0.7))

## Create additional grouping variables 
lutjanus_malabaricus_full_data$region_stage_year <- as.factor(paste(lutjanus_malabaricus_full_data$Region, lutjanus_malabaricus_full_data$MaturityStage, lutjanus_malabaricus_full_data$FishYear, sep = "_")) 
lutjanus_malabaricus_full_data$region_stage_cohort <- as.factor(paste(lutjanus_malabaricus_full_data$Region, lutjanus_malabaricus_full_data$MaturityStage, lutjanus_malabaricus_full_data$Cohort, sep = "_")) 

#################################################################################################################################
#################################################################################################################################

lutjanus_malabaricus_data <- lutjanus_malabaricus_full_data

unique(lutjanus_malabaricus_data$Method)
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>%
  mutate(Method = case_when(Method == "" ~ NA, 
                            Method == "Angling" ~ "Line",
                            TRUE ~ Method))

unique(lutjanus_malabaricus_data$Region)
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>% subset(select = -c(GonadWt))
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>% subset(select = -c(Depth))
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>% subset(select = -c(Fishery))
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>% subset(select = -c(Source))

sum(is.na(lutjanus_malabaricus_data$Site))
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>% drop_na(Site) 

table(lutjanus_malabaricus_data$Site)

## For reading 1 
#AAC_Measurements <- lutjanus_malabaricus_data %>% distinct(SampleID, Region, AAC) 
#write.csv(AAC_Measurements, 'LM_AAC_Measurements_3.csv')

### Split the data frame into the different geographic regions 
unique(lutjanus_malabaricus_data$Region)
lutjanus_malabaricus_data <- lutjanus_malabaricus_data %>%
  mutate(Region = case_when(Region == "NW Australia" ~ "North Western Australia", 
                            TRUE ~ Region))
lutjanus_malabaricus_data

#################################################################################################################################
#################################################################################################################################

## Following our discussion on 31 May 2024, we thought it would be best to split my models into region specific ones. 
## This is because we have unequal sample sizes for Kimberley versus Malacca Straits versus Riau Archipelago.
## The underlying driver of growth would then be explained by Kimberley, which has a much longer time series, compared to the Southeast Asian samples. 

Lmalabaricus_Australia <- lutjanus_malabaricus_data %>% filter(Region == "North Western Australia")
Lmalabaricus_Asia <- lutjanus_malabaricus_data %>% filter(Region == "Riau Archipelago")

### Construct the basic intrinsic linear mixed models 
## Australian samples
M1a <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID), 
            Lmalabaricus_Australia, REML= T)

M1b <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID), 
            Lmalabaricus_Australia, REML= T)

M1c <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia, REML= T)

M1d <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia, REML= T)

M1e <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|Cohort), 
            Lmalabaricus_Australia, REML= T)

M1f <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Lmalabaricus_Australia, REML= T)

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
           Lmalabaricus_Australia, REML= F)

options(na.action=na.fail)
M1_dredge<-dredge(M1, trace=2) 

subset(M1_dredge)

## Optimal intrinsic mixed model structure comprises of
M1 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Lmalabaricus_Australia, REML= T)

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
sample_depth_1 <- Lmalabaricus_Australia %>% 
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

#### Predict how Lutjanus malabaricus growth from Australia varies with age ####
Lmalabaricus_Australia %>% 
  dplyr::slice(which.max(Age)) %>% 
  pull(Age)

## Age plot
Lmalabaricus_Australia_Age <- as.data.frame (Effect (c('sAge'), M1, xlevels = list(sAge=seq(-1.0326, 4.54, by=0.01)) )) 
mean_Australia_age <- mean(Lmalabaricus_Australia$Age) ; mean_Australia_age
sd_Australia_age <- sd(Lmalabaricus_Australia$Age) ; sd_Australia_age
Lmalabaricus_Australia_Age$age <- (Lmalabaricus_Australia_Age$sAge * sd_Australia_age + mean_Australia_age)

Lmalabaricus_Australia_Age$transfit<-exp(Lmalabaricus_Australia_Age$fit)
Lmalabaricus_Australia_Age$transupper<-exp(Lmalabaricus_Australia_Age$upper)
Lmalabaricus_Australia_Age$translower<-exp(Lmalabaricus_Australia_Age$lower)

plot2 <- Lmalabaricus_Australia_Age %>% 
  filter(age <= 28.0) %>%
  ggplot(aes(y=transfit, x=age)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age", y = "Predicted growth (mm)", fill = "Site", color = "Site")+ 
  scale_x_continuous(limits=c(0,30), breaks = c(0,5,10,15,20,25,30), labels = c('0','5','10','15','20','25','30'))+
  scale_y_continuous(limits=c(0.0,0.45), breaks = c(0.000, 0.200, 0.400), labels = c('0.0','0.2','0.4'))+
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
Lmalabaricus_Australia_AAC <- as.data.frame (Effect (c('sAAC'), M1, xlevels = list(sAAC=seq(-1.578, 2.69 ,by=0.01)) )) 
mean_Australia_AAC <- mean(Lmalabaricus_Australia$AAC) ; mean_Australia_AAC
sd_Australia_AAC <- sd(Lmalabaricus_Australia$AAC) ; sd_Australia_AAC
Lmalabaricus_Australia_AAC$AAC <- (Lmalabaricus_Australia_AAC$sAAC * sd_Australia_AAC + mean_Australia_AAC)

Lmalabaricus_Australia_AAC$transfit<-exp(Lmalabaricus_Australia_AAC$fit)
Lmalabaricus_Australia_AAC$transupper<-exp(Lmalabaricus_Australia_AAC$upper)
Lmalabaricus_Australia_AAC$translower<-exp(Lmalabaricus_Australia_AAC$lower)

plot3 <- Lmalabaricus_Australia_AAC %>% 
  filter(AAC <= 28.0) %>%
  ggplot(aes(y=transfit, x=AAC)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age-at-Capture", y = "Predicted growth (mm)")+ 
  scale_x_continuous(limits=c(0,30), breaks = c(0,5,10,15,20,25,30), labels = c('0','5','10','15','20','25','30'))+
  scale_y_continuous(limits=c(0.0,0.3), breaks = c(0.000, 0.100, 0.200, 0.300), labels = c('0.0','0.1','0.2','0.3'))+
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
            Lmalabaricus_Asia, REML= T)

M2b <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID), 
            Lmalabaricus_Asia, REML= T)

M2c <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|FishYear), 
            Lmalabaricus_Asia, REML= T)

M2d <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Asia, REML= T)

M2e <- lmer(log(Width) ~ sAge + sAAC + 
              (1|SampleID) + (1|Cohort), 
            Lmalabaricus_Asia, REML= T)

M2f <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Lmalabaricus_Asia, REML= T)

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
           Lmalabaricus_Asia, REML= F)

options(na.action=na.fail)
M2_dredge<-dredge(M2, trace=2) 

subset(M2_dredge)

## Optimal intrinsic mixed model structure comprises of
M2 <- lmer(log(Width) ~ sAge + sAAC + 
             (sAge|SampleID) + (1|FishYear), 
           Lmalabaricus_Asia, REML= T)

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
sample_depth_2 <- Lmalabaricus_Asia %>% 
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
  scale_y_continuous(limits=c(-0.4,0.4))+
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

#### Predict how Lutjanus malabaricus growth from Asia varies with age ####
Lmalabaricus_Asia %>% 
  dplyr::slice(which.max(Age)) %>% 
  pull(Age)

## Age plot
Lmalabaricus_Asia_Age <- as.data.frame (Effect (c('sAge'), M2, xlevels = list(sAge=seq(-0.805, 4.74 ,by=0.01)) )) 
mean_Asia_age <- mean(Lmalabaricus_Asia$Age) ; mean_Asia_age
sd_Asia_age <- sd(Lmalabaricus_Asia$Age) ; sd_Asia_age
Lmalabaricus_Asia_Age$age <- (Lmalabaricus_Asia_Age$sAge * sd_Asia_age + mean_Asia_age)

Lmalabaricus_Asia_Age$transfit<-exp(Lmalabaricus_Asia_Age$fit)
Lmalabaricus_Asia_Age$transupper<-exp(Lmalabaricus_Asia_Age$upper)
Lmalabaricus_Asia_Age$translower<-exp(Lmalabaricus_Asia_Age$lower)

plot5 <- Lmalabaricus_Asia_Age %>% 
  filter(age <= 11.0) %>%
  ggplot(aes(y=transfit, x=age)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age", y = "Predicted growth (mm)", fill = "Site", color = "Site")+ 
  scale_x_continuous(limits=c(0,12), breaks = c(0,2,4,6,8,10,12), labels = c('0','2','4','6','8','10','12'))+
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
        legend.position.inside = c(1,1),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
plot5

## AAC plot
Lmalabaricus_Asia_AAC <- as.data.frame (Effect (c('sAAC'), M2, xlevels = list(sAAC=seq(-1.185, 2.895 ,by=0.01)) )) 
mean_Asia_AAC <- mean(Lmalabaricus_Asia$AAC) ; mean_Asia_AAC
sd_Asia_AAC <- sd(Lmalabaricus_Asia$AAC) ; sd_Asia_AAC
Lmalabaricus_Asia_AAC$AAC <- (Lmalabaricus_Asia_AAC$sAAC * sd_Asia_AAC + mean_Asia_AAC)

Lmalabaricus_Asia_AAC$transfit<-exp(Lmalabaricus_Asia_AAC$fit)
Lmalabaricus_Asia_AAC$transupper<-exp(Lmalabaricus_Asia_AAC$upper)
Lmalabaricus_Asia_AAC$translower<-exp(Lmalabaricus_Asia_AAC$lower)

plot6 <- Lmalabaricus_Asia_AAC %>% 
  filter(AAC <= 11.0) %>%
  ggplot(aes(y=transfit, x=AAC)) +
  geom_line(linewidth=0.5) +
  geom_ribbon(aes(ymin=translower, ymax=transupper), linewidth=0.4, alpha=0.3, color=0)+
  labs(x = "Age-at-Capture", y = "Predicted growth (mm)")+ 
  scale_x_continuous(limits=c(0,12), breaks = c(0,2,4,6,8,10,12), labels = c('0','2','4','6','8','10','12'))+
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
        legend.position = "none",
        legend.title=element_text(size=18), 
        legend.text=element_text(size=14)
  )
plot6

ggpubr::ggarrange(plot5, plot6, labels = c("A", "B"))

plot7 <- ggplot() + 
  geom_line(aes(y = transfit, x = age, colour="Australia"), linewidth=0.5, data = Lmalabaricus_Australia_Age) + 
  geom_ribbon(aes(y = transfit, x = age, ymin = translower, ymax = transupper, fill="Australia"), linewidth=0.4, alpha=0.3, colour=0, data = Lmalabaricus_Australia_Age) +
  geom_line(aes(y = transfit, x = age, colour="Southeast Asia"), linewidth=0.5,data = Lmalabaricus_Asia_Age) + 
  geom_ribbon(aes(y = transfit, x = age, ymin = translower, ymax = transupper, fill="Southeast Asia"), linewidth=0.4, alpha=0.3, colour=0, data = Lmalabaricus_Asia_Age) + 
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
  geom_line(aes(y = transfit, x = AAC, colour="Australia"), linewidth=0.5, data = Lmalabaricus_Australia_AAC) + 
  geom_ribbon(aes(y = transfit, x = AAC, ymin = translower, ymax = transupper, fill="Australia"), linewidth=0.4, alpha=0.3, colour=0, data = Lmalabaricus_Australia_AAC) +
  geom_line(aes(y = transfit, x = AAC, colour="Southeast Asia"), linewidth=0.5,data = Lmalabaricus_Asia_AAC) + 
  geom_ribbon(aes(y = transfit, x = AAC, ymin = translower, ymax = transupper, fill="Southeast Asia"), linewidth=0.4, alpha=0.3, colour=0, data = Lmalabaricus_Asia_AAC) + 
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
                                   
                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Pilbara",
                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Pilbara", 
                                   
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
                                   
                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Pilbara",
                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Pilbara", 
                                   
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
#                                   latitude <= 6.0 & latitude >= 3.0 & longitude >= 97.5 & longitude <= 101.25 ~ "Pilbara",
#                                   latitude <= 3.0 & latitude >= 1.0 & longitude >= 99.5 & longitude <= 103.0 ~ "Pilbara", 
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

##save.image("lutjanus_malabaricus.RData")
##save.image("lutjanus_malabaricus_backup.RData")

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

oisst_data$FishYear <- as.character(as.factor(cut(oisst_data$date, breaks=as.Date(paste(years,"-10-01",sep="")), labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))))
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

precip_data$FishYear <- as.character(as.factor(cut(precip_data$date, breaks=as.Date(paste(years,"-10-01",sep="")), labels=paste(years[-length(years)], years[-length(years)]+1,sep="/"))))
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

#chla_data$FishYear <- as.character(as.factor(cut(chla_data$date, breaks=as.Date(paste(years,"-10-01",sep="")), labels=paste(years[-length(years)], years[-length(years)]+1,sep="/"))))
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
  filter(Site == "Riau Archipelago" | Site == "Pilbara" | Site == "Kimberley") %>% 
  rename(ymd = date) 
oisst_month_data[, temp := mean(temp) , by = .(Year, Month, Site)][] 
oisst_month_data <- oisst_month_data %>%
  distinct(temp, Year, Month, .keep_all = TRUE) 

precip_month_data <- precip_month_data %>% 
  filter(Site == "Riau Archipelago" | Site == "Pilbara" | Site == "Kimberley") %>% 
  rename(ymd = date) 
precip_month_data[, precip := mean(precip) , by = .(Year, Month, Site)][] 
precip_month_data <- precip_month_data %>%
  distinct(precip, Year, Month, .keep_all = TRUE) 

env_data <- merge(env_data, oisst_month_data[,c('ymd','Site','temp')], by = c('ymd')) 
env_data <- merge(env_data, precip_month_data[,c('ymd','Site','precip')], by = c('ymd','Site')) 

env_data$ymd <- ymd(env_data$ymd)
env_data$FishYear <- as.character(as.factor(cut(env_data$ymd, breaks=as.Date(paste(years,"-10-01",sep="")), labels=paste(years[-length(years)],years[-length(years)]+1,sep="/"))))
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

## Dipole Mode Index
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

Australia_env_data <- env_data %>% filter(Site == "Kimberley" | Site == "Pilbara")
Asia_env_data <- env_data %>% filter(Site == "Riau Archipelago")

corvif(Australia_env_data[, c(7:11)]) 
corvif(Asia_env_data[, c(7:11)]) 

## Check the environmental data frame
env_data <- as.data.frame(env_data)
env_data

#################################################################################################################################
#################################################################################################################################

### Prepare the data for the exploratory sliding window analysis
Lmalabaricus_Australia <- Lmalabaricus_Australia %>%
  rename(capture_year = Year, capture_month = Month, capture_day = Day)

## Dummy variable for the alignment of environmental variables
Lmalabaricus_Australia$Month <- "10" 
Lmalabaricus_Australia$Month <- as.numeric(Lmalabaricus_Australia$Month)
Lmalabaricus_Australia <- Lmalabaricus_Australia %>% mutate(ymd = make_date(FishYear, Month)) 

## Preliminary analyses of the data indicated that we have biological data which starts in 1976. However, our environmental data only starts in September 1981. 
## Truncate the data to 1982 to ensure that we have environmental data coverage for all increment years. 
Lmalabaricus_Australia <- Lmalabaricus_Australia %>% filter(ymd >= '1983-10-01')
Lmalabaricus_Australia # Note the new data frame name
## Note that we truncated the data to 1983 to allow us to explore if growth was affected by environmental conditions in 1981 (two years).

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3a <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia, REML= F)

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
                                 refday = c(1, 10),                                  # day and month respectively of the year from which the absolute window analysis will start
                                 cinterval = "month",                               # resolution at which climate window analysis will be conducted. 
                                 cdate = env_data$ymd,                              # climate date variable (dd/mm/yyyy)
                                 bdate = Lmalabaricus_Australia$ymd,
                                 cohort = Lmalabaricus_Australia$FishYear,       # cohort is required because our environmental data spans across two years (see Advanced Vignette for more information)
                                 spatial = list(Lmalabaricus_Australia$Site, env_data$Site))

## As there is little temporal variation as to when the increment was completed, I used the absolute sliding window to analyse the data. 
## This means that there is a fixed reference time for all climate windows. 
## The climate year is important for the alignment of biological data, but the month is irrelevant for the biological data, as we assumed that the opaque zone completion period is consistent every year.
## On the other hand, bdate represents the biological date variable for the alignment of cdate. 
## The day and month of bdate is irrelevant and is only required for the alignment of cdate. 

## Note: Delta AICc is generated in comparison with the base null model M3. 

### Sliding window data ####
## The 'combos' function provides an overview of the best fitted climate window, the delta AICc value, the start and end window date, and slope of the best window. 
M3a_sliding_window
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
                             refday = c(1, 10),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Lmalabaricus_Australia$ymd,
                             cohort = Lmalabaricus_Australia$FishYear,
                             spatial = list(Lmalabaricus_Australia$Site, env_data$Site),
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

M3a_final <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia, REML= T)

summary(M3a_final)
anova(M3a_final)
rsquared(M3a_final)

M3a_conf <- sim(M3a_final, 5000)
M3a_conf <- apply(M3a_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3a_conf

#################################################################################################################################
#################################################################################################################################

### Mature Lutjanus malabaricus from Australia ####
Lmalabaricus_Australia_adult <- Lmalabaricus_Australia
Lmalabaricus_Australia_adult <- Lmalabaricus_Australia_adult %>% 
  filter(MaturityStage == "Mature")

M1a_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID), 
                  Lmalabaricus_Australia_adult, REML= T)

M1b_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID), 
                  Lmalabaricus_Australia_adult, REML= T)

M1c_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID) + (1|FishYear), 
                  Lmalabaricus_Australia_adult, REML= T)

M1d_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|FishYear), 
                  Lmalabaricus_Australia_adult, REML= T)

M1e_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (1|SampleID) + (1|Cohort), 
                  Lmalabaricus_Australia_adult, REML= T)

M1f_adult <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|Cohort), 
                  Lmalabaricus_Australia_adult, REML= T)

bbmle::AICctab(M1a_adult,M1b_adult,M1c_adult,M1d_adult,M1e_adult,M1f_adult, base=T,logLik=T,weights=T) 

M1_adult <- lmer(log(Width) ~ sAge + sAAC + 
                   (sAge|SampleID) + (1|FishYear), 
                 Lmalabaricus_Australia_adult, REML= F)

options(na.action=na.fail)
M1_adult_dredge<-dredge(M1_adult, trace=2) 

subset(M1_adult_dredge)

M1_adult <- lmer(log(Width) ~ sAge + sAAC + 
                   (sAge|SampleID) + (1|FishYear), 
                 Lmalabaricus_Australia_adult, REML= T)

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3b <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia_adult, REML= F)

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
                                 refday = c(1, 10),
                                 cinterval = "month",  
                                 cdate = env_data$ymd,   
                                 bdate = Lmalabaricus_Australia_adult$ymd,
                                 cohort = Lmalabaricus_Australia_adult$FishYear,
                                 spatial = list(Lmalabaricus_Australia_adult$Site, env_data$Site))

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
                             refday = c(1, 10),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Lmalabaricus_Australia_adult$ymd,
                             cohort = Lmalabaricus_Australia_adult$FishYear,
                             spatial = list(Lmalabaricus_Australia_adult$Site, env_data$Site),
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

M3b_final <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|FishYear), 
                  Lmalabaricus_Australia_adult, REML= T)

summary(M3b_final)
anova(M3b_final)
rsquared(M3b_final)

M3b_conf <- sim(M3b_final, 5000)
M3b_conf <- apply(M3b_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3b_conf

#################################################################################################################################
#################################################################################################################################

### Immature Lutjanus malabaricus #### 
Lmalabaricus_Australia_juvenile <- Lmalabaricus_Australia
Lmalabaricus_Australia_juvenile <- Lmalabaricus_Australia_juvenile %>% 
  filter(MaturityStage == "Immature")

M1a_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (1|SampleID), 
                     Lmalabaricus_Australia_juvenile, REML= T)

M1b_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (sAge|SampleID), 
                     Lmalabaricus_Australia_juvenile, REML= T)

M1c_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (1|SampleID) + (1|FishYear), 
                     Lmalabaricus_Australia_juvenile, REML= T)

M1d_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (sAge|SampleID) + (1|FishYear), 
                     Lmalabaricus_Australia_juvenile, REML= T)

M1e_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (1|SampleID) + (1|Cohort), 
                     Lmalabaricus_Australia_juvenile, REML= T)

M1f_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                       (sAge|SampleID) + (1|Cohort), 
                     Lmalabaricus_Australia_juvenile, REML= T)

bbmle::AICctab(M1a_juvenile,M1b_juvenile,M1c_juvenile,M1d_juvenile,M1e_juvenile,M1f_juvenile, base=T,logLik=T,weights=T) 

M1_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                      (sAge|SampleID) + (1|Cohort), 
                    Lmalabaricus_Australia_juvenile, REML= F)

options(na.action=na.fail)
M1_juvenile_dredge<-dredge(M1_juvenile, trace=2) 

subset(M1_juvenile_dredge)

M1_juvenile <- lmer(log(Width) ~ sAge + sAAC + 
                      (sAge|SampleID) + (1|Cohort), 
                    Lmalabaricus_Australia_juvenile, REML= T)

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3c <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|Cohort), 
            Lmalabaricus_Australia_juvenile, REML= F)

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
                                 refday = c(1, 10),
                                 cinterval = "month",  
                                 cdate = env_data$ymd,   
                                 bdate = Lmalabaricus_Australia_juvenile$ymd,
                                 cohort = Lmalabaricus_Australia_juvenile$FishYear,
                                 spatial = list(Lmalabaricus_Australia_juvenile$Site, env_data$Site))

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
                             refday = c(1, 10),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Lmalabaricus_Australia_juvenile$ymd,
                             cohort = Lmalabaricus_Australia_juvenile$FishYear,
                             spatial = list(Lmalabaricus_Australia_juvenile$Site, env_data$Site),
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

M3c_final <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|Cohort), 
                  Lmalabaricus_Australia_juvenile, REML= T)

summary(M3c_final)
anova(M3c_final)
rsquared(M3c_final)

M3c_conf <- sim(M3c_final, 5000)
M3c_conf <- apply(M3c_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M3c_conf

#################################################################################################################################
#################################################################################################################################

#### Exploratory species-specific analysis using mean linear environmental variables ####
#### Base model ####
M3d <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Australia_juvenile, REML= F)

## Applying the sliding window to the environmental data ####
M3d_sliding_window <- slidingwin(baseline = M3d, 
                                 xvar = list(sst = env_data$s.temp,
                                             rain = env_data$s.precip,
                                             pdo = env_data$s.pdo,
                                             enso = env_data$s.enso,
                                             dmi = env_data$s.dmi),  
                                 type = "absolute",
                                 range = c(23, 0),    
                                 stat = c("mean"), 
                                 func = c("lin"),  
                                 refday = c(1, 10),
                                 cinterval = "month",  
                                 cdate = env_data$ymd,   
                                 bdate = Lmalabaricus_Australia_juvenile$ymd,
                                 cohort = Lmalabaricus_Australia_juvenile$FishYear,
                                 spatial = list(Lmalabaricus_Australia_juvenile$Site, env_data$Site))

M3d_sliding_window$combos

M3d_random_window <- randwin(baseline = M3d, 
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
                             refday = c(1, 10),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Lmalabaricus_Australia_juvenile$ymd,
                             cohort = Lmalabaricus_Australia_juvenile$FishYear,
                             spatial = list(Lmalabaricus_Australia_juvenile$Site, env_data$Site),
                             window= "sliding")

## Calculate the p value to establish significance 
pvalue(datasetrand = M3d_random_window[[1]], dataset = M3d_sliding_window[[1]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3d_random_window[[2]], dataset = M3d_sliding_window[[2]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3d_random_window[[3]], dataset = M3d_sliding_window[[3]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3d_random_window[[4]], dataset = M3d_sliding_window[[4]]$Dataset, metric = "AIC") 
pvalue(datasetrand = M3d_random_window[[5]], dataset = M3d_sliding_window[[5]]$Dataset, metric = "AIC") 

#################################################################################################################################
#################################################################################################################################

### Prepare the data for the exploratory sliding window analysis
Lmalabaricus_Asia <- Lmalabaricus_Asia %>%
  dplyr::rename(capture_year = Year, capture_month = Month, capture_day = Day)

## Dummy variable for the alignment of environmental variables
Lmalabaricus_Asia$Month <- "10" 
Lmalabaricus_Asia$Month <- as.numeric(Lmalabaricus_Asia$Month)
Lmalabaricus_Asia <- Lmalabaricus_Asia %>% mutate(ymd = make_date(FishYear, Month)) 
Lmalabaricus_Asia 

#### Base model ####
M4a <- lmer(log(Width) ~ sAge + sAAC + 
              (sAge|SampleID) + (1|FishYear), 
            Lmalabaricus_Asia, REML= F)

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
                                 refday = c(1, 10),                                 # day and month respectively of the year from which the absolute window analysis will start
                                 cinterval = "month",                               # resolution at which climate window analysis will be conducted. 
                                 cdate = env_data$ymd,                              # climate date variable (dd/mm/yyyy)
                                 bdate = Lmalabaricus_Asia$ymd,
                                 cohort = Lmalabaricus_Asia$FishYear,       # cohort is required because our environmental data spans across two years (see Advanced Vignette for more information)
                                 spatial = list(Lmalabaricus_Asia$Site, env_data$Site))

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
                             refday = c(1, 10),
                             cinterval = "month",
                             cdate = env_data$ymd, 
                             bdate = Lmalabaricus_Asia$ymd,
                             cohort = Lmalabaricus_Asia$FishYear,
                             spatial = list(Lmalabaricus_Asia$Site, env_data$Site),
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

M4a_final <- lmer(log(Width) ~ sAge + sAAC + 
                    (sAge|SampleID) + (1|FishYear), 
                  Lmalabaricus_Asia, REML= T)

summary(M4a_final)
anova(M4a_final)
rsquared(M4a_final)

M4a_conf <- sim(M4a_final, 5000)
M4a_conf <- apply(M4a_conf@fixef,2, function(x) quantile(x,probs=c(0.5, 0.025, 0.975))) ; M4a_conf

#################################################################################################################################
#################################################################################################################################

#### Study region ####
## New update as of 12 November 2023, Stadia Maps requires an API Key.
## I will be using Google Maps instead. 
remotes::install_github("dkahle/ggmap")
library(ggmap)
register_google(key = "AIzaSyBBfJPr7DR7GYl6FxeslPLBqQE7cia55j4", write = TRUE)
lutjanus_malabaricus_study_region <- get_map(location = c(110, 14.25), maptype='satellite', source="google", zoom=3)

eez_shp <- readOGR("../../Climate/EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp") %>% spTransform("+proj=longlat +ellps=WGS84")
eez_shp <- ggplot2::fortify(eez_shp)

ggmap(lutjanus_malabaricus_study_region) +
  #geom_path(data=eez_shp, aes(long, lat, group=group), linewidth= 0.075)+ # Exclusive economic zones
  xlab("Longitude") + ylab("Latitude")+
  scale_x_continuous(limits = c(80, 140), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-25, 10), expand = c(0, 0)) +
  
  #annotate('point', y= 3.5, x= 100.0, color='firebrick3', size=20, alpha=0.4)+
  annotate('point', y= 2.2, x= 107.0, size= 25, color='firebrick3', alpha= 0.5) + 
  annotate('point', y= -18.5, x= 117.0, size= 25, color='firebrick3', alpha= 0.5) + 
  annotate('point', y= -14.0, x= 124.5, size= 25, color='firebrick3', alpha= 0.5) + 
  
  #annotate('text', y=  3.5, x= 100.0, color='grey98', size= 3, label='Malacca\nStraits', fontface = 2) +
  annotate('text', y=  2.2, x= 107.0, color='grey98', size= 3, label='Riau\nArchipelago', fontface= 2) +
  annotate('text', y= -18.5, x= 117.0, color='grey98', size= 3, label='Pilbara', fontface= 2) +
  annotate('text', y= -14.0, x= 124.5, color='grey98', size= 3, label='Kimberley', fontface= 2) +
  
  theme(plot.title = element_text(hjust=0.5, size=18), 
        axis.title=element_text(size=18),
        axis.text = element_text(size=16),
        legend.position='none')

indo_pacific <- get_map(location = c(90, 1), maptype='satellite', source="google", zoom=2)

ggmap(indo_pacific) + 
  xlab("Longitude") + ylab("Latitude") +
  scale_x_continuous(limits = c(20, 200), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-45, 45), expand = c(0, 0)) 

#################################################################################################################################
#################################################################################################################################

## Calculation of IAPE and ACV for Lutjanus malabaricus
SEA_age <- readxl::read_xlsx("Julio_Snapper_Readings_SEA.xlsx", 1, col_names = F)
SEA_age <- SEA_age %>%
  subset(select = -c(2,3,5,7)) %>%
  dplyr::rename("SampleID" = "...1", "Reader1" = "...4", "Reader2" = "...6") 
SEA_age = SEA_age[-1,]
SEA_age$Reader1 <- as.numeric(SEA_age$Reader1)
SEA_age$Reader2 <- as.numeric(SEA_age$Reader2)
SEA_age$diff <- SEA_age$Reader1 - SEA_age$Reader2

WA_age <- readxl::read_xlsx("Corey_Snapper_Readings_WA.xlsx", 1, col_names = F)
WA_age <- WA_age %>%
  subset(select = -c(3)) %>%
  dplyr::rename("SampleID" = "...1", "Reader1" = "...2","Reader2" = "...4" ) 
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

#save.image("lutjanus_malabaricus.RData")

#################################################################################################################################
#################################################################################################################################
