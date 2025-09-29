## Will Patrone
## UW SAFS
## Metier Analysis for West Coast Groundfish

## Set path for data
#Will's Path
#setwd("C:/Users/Will Patrone/Documents/Master's Thesis/Metier Analysis")
#Andre's Path
#setwd()


#load potentially useful packages



vmstoolsPackages <-  c("data.table","doBy","lubridate","sf","mixtools","segmented", 
                       "cluster", "maps", "tidyverse", "devtools", "reshape2", "SOAR",
                       "FactoMineR", "writexl", "clustree", "mapdata", "PBSmapping", 
                       "writexl", "readxl")
for(i in vmstoolsPackages)       try(install.packages(pkgs=i,repos=getOption("repos")))
install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
install.packages("C:/Users/Will Patrone/Documents/test metier/vmstools_0.76.zip")




#Once you've installed all secondary packages install VMStools using

#In office


#At home
#install.packages("C:/Users/Will Patrone/Documents/Master's Thesis/Metier Analysis/vmstools_0.76")

library(tidyverse)
library(devtools)
library(vmstools)

# Loop through each package and load it
for (pkg in vmstoolsPackages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}





## load data
#load("PacFIN.GRND.CompFT.10.Jul.2024.Rdata")





#Get list of species in analysis


data = catch.pacfin[c("LANDING_YEAR", "LANDING_MONTH", "GEAR_CODE", "GEAR_NAME",
                      "PACFIN_GEAR_DESCRIPTION", "PACFIN_CATCH_AREA_DESCRIPTION",
                      "PACFIN_GROUP_CATCH_AREA_CODE", "COUNCIL_CODE", "PORT_CODE",
                      "PACFIN_PORT_NAME", "COUNTY_NAME", "COUNTY_STATE", 
                      "SUBREGION_NAME", "REGION_NAME", 
                      "PACFIN_SPECIES_CODE", "PACFIN_SPECIES_SCIENTIFIC_NAME",
                      "PACFIN_SPECIES_COMMON_NAME", "LANDED_WEIGHT_LBS", 
                      "LANDED_WEIGHT_MTONS", "PRICE_PER_POUND", "EXVESSEL_REVENUE")]

dummy = catch.pacfin[c("PACFIN_SPECIES_SCIENTIFIC_NAME",
          "PACFIN_SPECIES_COMMON_NAME",
          "LANDED_WEIGHT_MTONS")]



dummy = filter(dummy, !str_detect(PACFIN_SPECIES_SCIENTIFIC_NAME, regex("\\bSPP", ignore_case = TRUE)))
dummy = filter(dummy, !str_detect(PACFIN_SPECIES_SCIENTIFIC_NAME, "N/A"))

dummy = dummy[c("PACFIN_SPECIES_COMMON_NAME",
                       "LANDED_WEIGHT_MTONS")]

#Condense dummy to relevant species

allcatchmt= aggregate(.~PACFIN_SPECIES_COMMON_NAME, data = dummy, FUN = sum)

allcatchmt = allcatchmt[order(allcatchmt$LANDED_WEIGHT_MTONS, decreasing = TRUE),]

choosespp = filter(allcatchmt, LANDED_WEIGHT_MTONS > 5000)

sppname = unique(choosespp$PACFIN_SPECIES_COMMON_NAME)

data$PACFIN_GEAR_DESCRIPTION <- ifelse(data$PACFIN_GEAR_DESCRIPTION %in% c("SHRIMP TRAWL, SINGLE RIGGED",
                                                                           "SHRIMP TRAWL, SINGLE OR DOUBLE RIG", 
                                                                           "DANISH/SCOTTISH SEINE (TRAWL)", 
                                                                           "BEAM TRAWL",
                                                                           "PAIR TRAWL"), "OTHER TRAWL GEAR", 
                                       data$PACFIN_GEAR_DESCRIPTION)

dummy = catch.pacfin[c("PACFIN_GEAR_DESCRIPTION", "LANDED_WEIGHT_MTONS")]

test = aggregate(.~PACFIN_GEAR_DESCRIPTION, data = dummy, FUN = sum)
test = test[order(test$LANDED_WEIGHT_MTONS, decreasing = TRUE),]
choosegear = filter(test, LANDED_WEIGHT_MTONS > 25000) 


data = filter(data, PACFIN_GEAR_DESCRIPTION %in% choosegear$PACFIN_GEAR_DESCRIPTION)
data = filter(data, PACFIN_SPECIES_COMMON_NAME %in% choosespp$PACFIN_SPECIES_COMMON_NAME)
data$shotID = paste(data$LANDING_MONTH, data$PACFIN_PORT_NAME, data$PACFIN_GEAR_DESCRIPTION, sep = "_")
data = relocate(data, shotID)


plzwork = data %>% dplyr::select(shotID, PACFIN_SPECIES_COMMON_NAME, PRICE_PER_POUND, EXVESSEL_REVENUE)

plzwork$product = plzwork$PRICE_PER_POUND * plzwork$EXVESSEL_REVENUE

useme_df = matrix(0, nrow = length(unique(plzwork$shotID)), ncol = length(unique(plzwork$PACFIN_SPECIES_COMMON_NAME)),
                  dimnames = list(unique(plzwork$shotID), unique(plzwork$PACFIN_SPECIES_COMMON_NAME)))

for (i in 1:nrow(useme_df)) {
  row = plzwork$shotID[i]
  col = plzwork$PACFIN_SPECIES_COMMON_NAME[i]
  value = plzwork$product[i]
  useme_df[row, col] = value
}

casting = dcast(formula = shotID ~ PACFIN_SPECIES_COMMON_NAME, data = plzwork, 
                fun.aggregate = sum, na.rm=TRUE)

Step1 <- extractTableMainSpecies(dat = casting, sppname,
                                 paramTotal=100, paramLogevent = 100)


head(Step1)
## Get PCA
rowNamesSave=row.names(Step1)
row.names(Step1)=1:nrow(Step1)



analysisName <- paste(data$LANDING_MONTH, data$PACFIN_PORT_NAME, data$COUNTY_STATE, data$PACFIN_GEAR_DESCRIPTION, sep="_")


critPca="NO_PCA"

Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)

Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL)



originalcluster = Step3$clustering

###################################################
summarized_data <- casting %>%
  group_by(clusternumber) %>%                             # Group by the grouping variable
  summarise(across("ARROWTOOTH FLOUNDER":"YELLOWTAIL ROCKFISH", 
                   sum, na.rm = TRUE)) 

work <- sweep(summarized_data[, -1], 2, colSums(summarized_data[, -1]), FUN = "/")

work$clusternumber = 1:17

work = relocate(work, clusternumber)

#write_xlsx(summarized_data, "Patrone Metier Clusters Landed Weight By Cluster.xlsx")




tidy_data <- work %>%
  pivot_longer(cols = -clusternumber, names_to = "Category", values_to = "Value")

plot1 = ggplot(tidy_data, aes(x = clusternumber, y = Value, fill = clusternumber)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = seq(1, 17, by = 1))+
  facet_wrap(~ Category, scales = "free_y") +
  theme_minimal() +
  labs(x = "Cluster Number", y = "Proportion of Catch in Each Cluster", title = "Proportion of Catch by Species in Each Cluster") +
  theme( 
    strip.text = element_text(size = 10),  # Reduce the size of facet titles
    legend.position = "none"              # Hide the legend if not needed
  )


ggsave("Proportion of Catch by Species in Each Cluster.jpg", 
       plot = plot1, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)   # Resolution in DPI



tidy_data$SpeciesNumber <- as.numeric(factor(tidy_data$Category))  # Assign numbers to species

# Create the plot
plot2 = ggplot(tidy_data, aes(x = SpeciesNumber, y = Value, fill = factor(SpeciesNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clusternumber, scales = "free_x") +
  scale_x_continuous(
    breaks = 1:length(unique(tidy_data$Category)),  # Set numeric breaks
    labels = 1:length(unique(tidy_data$Category))  # Use numbers as labels
  ) +
  scale_fill_manual(
    values = scales::hue_pal()(length(unique(tidy_data$Category))),  # Use default ggplot colors
    labels = paste(1:length(unique(tidy_data$Category)), unique(tidy_data$Category), sep = ": ")  # Combine number and species in the legend
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),  # Customize x-axis text size
    legend.title = element_text(size = 10)
  ) +
  labs(
    x = "Species (Numbered)",
    y = "Proprotion of Catch in Each Metier",
    title = "Proportion of Catch in Each Metier by Species",
    fill = "Species Key"  # Legend title
  )


ggsave("Proportion of Catch by Species for Each Cluster.jpg", 
       
       plot = plot2, 
       width = 30,  # Width in inches
       height = 20, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)   # Resolution in DPI


#########TEsting


data_long <- casting %>%
  pivot_longer(
    cols = -c(clusternumber, shotID),
    names_to = "species",
    values_to = "catch"
  )

# Calculate proportions
data_long <- data_long %>%
  group_by(clusternumber) %>%
  mutate(proportion = catch / sum(catch))

# Plot the stacked bar chart
stackbar = ggplot(data_long, aes(x = factor(clusternumber), y = proportion, fill = species)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Proportion of Catch by Speices for Each Cluster",
    x = "Cluster Number",
    y = "Proportion of Total Catch In Metier",
    fill = "Species"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()

ggsave("Proportion of Catch by Speices for Each Cluster.jpg",  
       plot = stackbar, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)





#####################################################


write_xlsx(casting, "Patrone Metier Clusters.xlsx")


Step3$medoids



sciname = unique(data$PACFIN_SPECIES_COMMON_NAME)
comname = unique(data$PACFIN_SPECIES_COMMON_NAME)

# propor <- postclustdata %>%
#   group_by(Port, GEAR, clusternumber) %>%
#   summarise(Count = n()) %>%
#   mutate(Proportion = Count / sum(Count))
# 
# print(propor)
# 
# write_xlsx(propor, "Patrone Metier Proportions Test.xlsx")


########### read in excel file separated ############

geardata = read_excel("Charts by gear Data.xlsx")


propor <- geardata %>%
  group_by(clusternumber, Gear) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(propor)


gearbar = ggplot(propor, aes(x = factor(clusternumber), y = Proportion, fill = Gear)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Proportion of Gear Used in Each Cluster",
    x = "Cluster Number",
    y = "Proportion of Gear In Metier",
    fill = "Gear"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()

ggsave("Proportion of Gear Used in Each Cluster.jpg",  
       plot = gearbar, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)


#################
portdata = read_excel("cluster By Port for graph.xlsx")
names(portdata)


filterdata = data[c("PACFIN_PORT_NAME",
                            "COUNTY_STATE")]

unqiuedf <- unique(filterdata)

dataframe1 <- merge(portdata, unqiuedf, by = "PACFIN_PORT_NAME", all.x = TRUE)

# dataframe1 <- portdata %>%
#   left_join(filterdata %>% select(PACFIN_PORT_NAME, COUNTY_STATE), by = "PACFIN_PORT_NAME")
# head(dataframe1)

dataframe1 = dataframe1 %>% drop_na()

propor <- dataframe1 %>%
  group_by(clusternumber, COUNTY_STATE) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(propor)


statebar = ggplot(propor, aes(x = factor(clusternumber), y = Proportion, fill = COUNTY_STATE)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Proportion of Cacth for Each State Each Cluster",
    x = "Cluster Number",
    y = "Proportion of Catch in State in Metier",
    fill = "State"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()


ggsave("Proportion of Catch By State Each Cluster.jpg",  
       plot = statebar, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)





donedata = read_excel("Patrone Metier Clusters Separated withour state.xlsx")
names(donedata)

uniquegear = unique(donedata$Gear)

geardesc = c("Line or Pole", "Line or Pole", "Bottom Trawl", "Set Net and Fish Pot", "Bottom Trawl", 
             "Midwater Trawl", "Bottom Trawl", "Bottom Trawl", "Bottom Trawl", "Bottom Trawl", "Bottom Trawl", 
             "Set Net and Fish Pot")


gearmatrix = matrix(c(uniquegear, geardesc), ncol = 2)
geardf = as.data.frame(gearmatrix)
names(geardf) = c("Gear", "GearDescription")

dataframe2 <- merge(donedata, geardf, by = "Gear", all.x = TRUE)
dataframe2$Sum = NULL

summarizeddata <- dataframe2 %>%
  dplyr::select(!c(shotID, Gear)) %>% 
  group_by(GearDescription, clusternumber) %>%
  summarise(across("ARROWTOOTH FLOUNDER":"YELLOWTAIL ROCKFISH", 
                   sum, na.rm = TRUE))
    .groups = "drop" # Ensure the result is ungrouped
  )

summarizeddata[22,24] = summarizeddata[22,24]*0.01

df_long <- summarizeddata %>%
  pivot_longer(
    cols = 3:33,                
    names_to = "Species",      
    values_to = "Catch"         
  )


geargroupbar = ggplot(df_long, aes(x = Species, y = Catch, fill = as.factor(clusternumber))) +
  geom_bar(stat = "identity") +  # `position = "fill"` makes it 100% stacked
  facet_wrap(~ GearDescription) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  labs(
    x = "Species",
    y = "Proportion",
    fill = "Cluster Number",
    title = "Proportion of Catch By Gear Tyoe and Cluster for Each Species"
  )


ggsave("Proportion of Catch By Gear Tyoe and Cluster for Each Species.jpg",  
       plot = geargroupbar, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)
