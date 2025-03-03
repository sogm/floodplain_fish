setwd("..../JSDM_SMDB/")# to be defined by the user

library(Hmsc)
library(readr)
library(tidyverse)
library(readxl)
library(writexl)
library(corrplot)
library(ape)
library(coda)
library(sf)
library(caret)

localDir = "."
dataDir = file.path(localDir, "data")# the path to a folder named data
modelDir = file.path(localDir, "models")# The path to a folder name models
#if(!dir.exists(modelDir)) dir.create(modelDir)

# Exploratory analysis
spp_env <- read_csv("data/fish_env_data_jsdm.csv")
TrData<- read_xlsx("data/fi_traits.xlsx")
taxonomy <- read_xlsx("data/fi_phylo.xlsx")

#Filter the data to accommodate our research interest only the Southern Basin Catchment SBC
south_basin_catch_mdb_cewovall <-c("Lower Murray", "Wimmera", "Avoca",
                                   "Loddon", "Campaspe", "Goulburn", 
                                   "Broken", "Ovens", "Central Murray", 
                                   "Lower Darling", "Edward Wakool",
                                   "Murrumbidgee", "Mitta Mitta", "Kiewa", "Upper Murray")

wetland_type <- c("Floodplain", "Palustrine", "Lacustrine")

spp_env_floodplains <- spp_env %>% filter(SystmTy %in% wetland_type) %>% 
                       filter(CEWOVll %in% south_basin_catch_mdb_cewovall) # 83927 Only floodplain wetlands in the SMDB

#Select based on the list wetlands in floodplains in a 2.5 km buffer within the SMDB
floodplains_mask_5000 <- read_sf("/fish_in_buffer_5000.shp")# Buffer of the area of interest for the analysis
d1 <- spp_env_floodplains %>% filter(UID %in% floodplains_mask_5000$UID)#66411

#Selecting variables of interest
d1 <- d1 %>% dplyr::select(UID, lon, lat, CEWOVll,# spatial random effects, study design. 
                                nearest_mapped_dist, n_wetlands_pwet, area_wetlands_pwet,  # spatially relevant variables for metapopulation/patch approach
                                wet_prop_mean, area_mean_max, STRANNTEMP, STRELEMEAN, STRANNRAIN, RUNANNMEAN, VALLEYSLOPE,# environmental variables
                                'Perca fluviatilis', 
                                'Philypnodon grandiceps',
                                'Cyprinus carpio',            
                                'Gambusia holbrooki',
                                'Nannoperca australis',
                                 'Retropinna semoni', 
                                'Craterocephalus fluviatilis',
                                'Gadopsis marmoratus', 
                                'Mogurnda adspersa', 
                                'Hypseleotris klunzingeri', 
                                'Galaxias rostratus',         
                                'Misgurnus anguillicaudatus')# species records 

d1 <- d1 %>% mutate(n_wetlands_pwet = if_else(is.na(n_wetlands_pwet), 0, n_wetlands_pwet),
                    area_wetlands_pwet = if_else(is.na(area_wetlands_pwet), 0, area_wetlands_pwet))

#Exploratory
f_data_long <- d1 [,c(5:14)]%>%
  pivot_longer(cols = everything(),
               names_to ="variable",
               values_to = "value")

# Create a histogram for all numeric variables in one plot
f_histograms <- ggplot(f_data_long, aes(x = value)) +
  geom_histogram(bins = 30, color = "black", fill = "lightblue") +
  facet_wrap(~variable, scales = "free", ncol = 5, nrow= 2) +
  labs(title = "Histograms of Numeric Variables in Ovens",
       x = "Value",
       y = "Frequency") +
  theme_minimal()
f_histograms

#Correlations
enviro_cor<- d1[,c(5:14)]
enviro_cor2 <- na.omit(enviro_cor)#59416
c <-cor(enviro_cor2)
corrplot(c, method = 'number', order = 'alphabet')

my_cor <- findCorrelation(
  c,
  cutoff = 0.6,
  verbose = TRUE,
  names = TRUE,
  exact = TRUE # Whether average correlation should be recomputed at each step
)

my_cor# remove elevation, precipitation and valley slope

#make this example reproducible
set.seed(1)
d1$id <- 1:nrow(d1) # 66411 values
d2 <- na.omit(d1)# generates 1032 observations with all the variables in the MDB.
#The more environmental variables in the model the lower the number of samples, since there are many NA that can not be used for modelling.

# Organize the community data in the matrix Y
Y <- d2[,c(15:21, 23:26)] %>% 
     mutate(across(.cols = c(1:11), .fns = function(x) ifelse( x >  0, 1, 0)))# IDEM

Y <- Y %>%  rename(Perca_fluviatilis = 'Perca fluviatilis', 
                   Philypnodon_grandiceps  = 'Philypnodon grandiceps' ,
                   Cyprinus_carpio = 'Cyprinus carpio', 
                   Gambusia_holbrooki = 'Gambusia holbrooki', 
                   Nannoperca_australis = 'Nannoperca australis', 
                   Retropinna_semoni = 'Retropinna semoni', 
                   Craterocephalus_fluviatilis = 'Craterocephalus fluviatilis', 
                   Mogurnda_adspersa = 'Mogurnda adspersa', 
                   Hypseleotris_spp = 'Hypseleotris klunzingeri', 
                   Galaxias_rostratus = 'Galaxias rostratus',
                   Misgurnus_anguillicaudatus = 'Misgurnus anguillicaudatus')

wet<-d2 %>% filter(`Galaxias rostratus`> 0)

mean(wet$area_mean_max)

#species prevalence
fish_prevalence <- data.frame(prevalence = colMeans(Y), n_wetlands = colSums(Y))
fish_prevalence$spp_names <-row.names(fish_prevalence)

fish_prevalence$spp_names <- factor(fish_prevalence$spp_names, levels = c("Cyprinus_carpio",
                                                                          "Gambusia_holbrooki",
                                                                          "Perca_fluviatilis",
                                                                          "Misgurnus_anguillicaudatus",
                                                                          "Retropinna_semoni",
                                                                          "Philypnodon_grandiceps",
                                                                          "Hypseleotris_spp",
                                                                          "Nannoperca_australis", 
                                                                          "Craterocephalus_fluviatilis",
                                                                          "Galaxias_rostratus",
                                                                          "Mogurnda_adspersa"))

#spp_prevalence plot
ggplot(fish_prevalence, aes(fill = spp_names, y= prevalence, x= spp_names)) + 
  geom_bar(stat="identity", colour="black", width =0.7, position=position_dodge(width =.85))+
  scale_y_continuous(breaks=seq(0,1,0.25), limits=c(0,1))+
  geom_text(aes(label = n_wetlands, vjust = 0.5, hjust = -1))+
  ylab("\nSpecies")+
  xlab("Prevalence\n")+
  coord_flip()+
  
  theme(legend.position = "none",# similar conditions as theme classic.
        axis.title.x = element_text(color="black", size=16, face=1),# color of the text for the axes x
        axis.title.y = element_text(color="black", size=16, face=1),# color of the text for the axes y
        axis.text.x = element_text(colour = "black",size=14, face=1),# color of the text for the labels of the data in axes x
        axis.text.y = element_text(colour = "black",size=14, face=1), # color of the text for the labels of the data in axes y
        panel.background = element_blank(),# removes the gray background on the panel and leaves it white
        axis.line = element_line(colour = "black"),# adds black solid lines around the panel
        panel.grid.major = element_blank(),#provides white background in each element of the panel
        panel.grid.minor = element_blank(),
        strip.background = element_blank())# removes the gray background on the name of the reefs

# Organize the environmental data into a dataframe XData
XData = d2 %>% dplyr::select(nearest_mapped_dist, n_wetlands_pwet, area_wetlands_pwet, wet_prop_mean, area_mean_max, STRANNTEMP, RUNANNMEAN) %>% 
  rename (mean_wet_time =  wet_prop_mean,
                           wetland_area = area_mean_max,
                           temperature = STRANNTEMP,
                           runoff = RUNANNMEAN,
                           distance_to_river = nearest_mapped_dist,
                           number_permanent_wetlands = n_wetlands_pwet,
                           area_permanent_wetlands = area_wetlands_pwet)

# Define the environmental model through XFormula
 XFormula = ~ mean_wet_time + wetland_area + temperature + runoff + distance_to_river + number_permanent_wetlands + area_permanent_wetlands 

library(summarytools)
summary_table_1 <- XData %>% summarytools::descr()
names <- rownames(summary_table_1)
summary_table_1 <- as_data_frame(summary_table_1) 
summary_table_1 <-cbind(names, summary_table_1)

# Organize the trait data into a dataframe TrData
TrData <- TrData %>% dplyr:: select (1:14) %>% 
            mutate(resilience =  factor(resilience),
            water_column = factor(water_column),
            residency = factor(residency)) %>% 
            data.frame()

rownames(TrData) = colnames(Y)

TrFormula = ~ length + resilience + residency + length_maturity + age_maturity + fecundity_max + temp_m 

# Set up a phylogenetic (or taxonomic tree) as myTree
taxonomy <- taxonomy %>%  mutate(order = factor(order),
                                 family = factor (family),
                                 genus = factor (genus),
                                 species = factor (species))

myTree = as.phylo(~order/family/genus/species, data = taxonomy, collapse = FALSE)
myTree$edge.length = rep(1, length(myTree$edge))
myTree = myTree
plot(myTree)

# Define the studyDesign as a dataframe 

#First, set up coordinates
xy = as.matrix(cbind(d2[,2], d2[, 3]))
rownames(xy) = unique(d2$UID)# rows names must match the names of sampling units in this case wetlands
colnames(xy) = c("x", "y")

studyDesign = data.frame(catchment = as.factor(d2$CEWOVll),#catchments
                         location = as.factor(d2$UID))#locatiion

# Set up the random effects
rL.location = Hmsc::HmscRandomLevel(sData = xy)
rL.catchment = Hmsc::HmscRandomLevel(units = levels(studyDesign$catchment))

# Use the Hmsc model constructor to define a model

m_full = Hmsc(Y=Y,
         distr="probit",
         XData = XData,  XFormula=XFormula,
         TrData = TrData, TrFormula = TrFormula,
         phyloTree = myTree,
         studyDesign = studyDesign,
         ranLevels=list(location = rL.location, catchment = rL.catchment))

 models = list(m_full)
 names(models) = c("m_full")
 save(models, file = file.path(modelDir, "unfitted_models.RData"))

for(i in 1:length(models)){
 print(i)
 sampleMcmc(models[[i]],samples=2)
}