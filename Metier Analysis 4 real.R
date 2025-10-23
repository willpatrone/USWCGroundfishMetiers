## Will Patrone
## UW SAFS
## Metier Analysis for West Coast Groundfish

## Set path for data
#Will's Path
#setwd("C:/Users/Will Patrone/Documents/Master's Thesis/Metier Analysis")
#Andre's Path
#setwd()


#load potentially useful packages



vmstoolsPackages=c("data.table","doBy","lubridate","sf","mixtools","segmented", 
                       "cluster", "maps", "tidyverse", "devtools", "reshape2", "SOAR",
                       "FactoMineR", "writexl", "clustree", "mapdata", "PBSmapping", 
                       "writexl", "readxl", "patchwork")
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

data$PACFIN_GEAR_DESCRIPTION=ifelse(data$PACFIN_GEAR_DESCRIPTION %in% c("SHRIMP TRAWL, SINGLE RIGGED",
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

Step1=extractTableMainSpecies(dat = casting, sppname,
                                 paramTotal=100, paramLogevent = 100)


head(Step1)
## Get PCA
rowNamesSave=row.names(Step1)
row.names(Step1)=1:nrow(Step1)



analysisName=paste(data$LANDING_MONTH, data$PACFIN_PORT_NAME, data$COUNTY_STATE, data$PACFIN_GEAR_DESCRIPTION, sep="_")


critPca="NO_PCA"

Step2=getTableAfterPCA(Step1,analysisName,pcaYesNo="nopca",criterion=NULL)

Step3=getMetierClusters(Step1,Step2,analysisName=analysisName,methMetier="clara",param1="euclidean",param2=NULL)



casting$clusternumber = Step3$clustering



clusterstometiers=c( "1" = "A", "2" = "B","3" = "B", "4" = "C", "5" = "D", "6" = "B", 
                        "7" = "A","8" = "C", "9" = "C", "10" = "D", "11" = "B", "12" = "B",
                        "13" = "B", "14" = "C", "15" = "B", "16" = "B", "17" = "A" )

metiercasting = casting

metiercasting$metier=clusterstometiers[as.character(metiercasting$clusternumber)]

metiercasting = metiercasting[,c(1:32, 34)]

clusterpropor = c(0.32, 0.07, 0.20, 5.62, 24.1, 10.9, 3.40, 0.36, 17.9, 3.49, 
                  0.20, 0.01, 0.003, 0.59, 0.05, 0.01, 33.0)
clusterdf = data.frame(cluster = 1:17, proportion = clusterpropor)

clusterproporbar = ggplot(clusterdf, aes(x = factor(cluster), y = proportion)) +
  geom_col(fill = "gray30") +  
  scale_y_reverse(limits = c(35, 0)) +  
  scale_x_discrete(breaks = 1:17) +
  labs(x = "", y = "Proportion of EVR (%)", title = "Proportion of Total EVR in the Fishery") +
  theme_minimal()+
  theme(strip.text = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),      
        axis.title.x = element_text(size = 20),     
        axis.title.y = element_text(size = 17),     
        legend.position = "none")
  
metierpropor = c(36.7, 11.4, 24.5, 27.6)
metierdf = data.frame(metier = c("A", "B", "C", "D"), proportion = metierpropor)

metierproporbar = ggplot(metierdf, aes(x = factor(metier), y = proportion)) +
  geom_col(fill = "gray30") +  
  scale_y_reverse(limits = c(40, 0)) +  
  scale_x_discrete() +
  labs(x = "", y = "Proportion of EVR (%)", title = "Proportion of Total EVR in the Fishery") +
  theme_minimal()+
  theme(strip.text = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),      
        axis.title.x = element_text(size = 20),     
        axis.title.y = element_text(size = 17),     
        legend.position = "none")
metierproporbar

###################################################
summarized_data=casting %>%
  group_by(clusternumber) %>% 
  summarise(across("ARROWTOOTH FLOUNDER":"YELLOWTAIL ROCKFISH", 
                   sum, na.rm = TRUE)) 

work=sweep(summarized_data[, -1], 2, colSums(summarized_data[, -1]), FUN = "/")

work$clusternumber = 1:17

work = relocate(work, clusternumber)




tidy_data=work %>%
  pivot_longer(cols = -clusternumber, names_to = "Category", values_to = "Value")



plot1 = ggplot(tidy_data, aes(x = clusternumber, y = Value)) +
  geom_bar(stat = "identity", fill = "black") +
  scale_x_continuous(breaks = seq(1, 17, by = 1))+
  facet_wrap(~ Category, scales = "free_y", ncol = 4) +
  theme_minimal() +
  labs(x = "Cluster Number", y = "Proportion of EVR in Each Cluster") +
  theme_minimal()+
  theme(strip.text = element_text(size = 18), 
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 17),      
    axis.title.x = element_text(size = 24),     
    axis.title.y = element_text(size = 24),     
    legend.position = "none")

plot1

ggsave("Figure A3 - Proportion of EVR by Species in Each Cluster.jpg", 
       plot = plot1, 
       width = 27.5, 
       height = 15, 
       units = "in", 
       dpi = 300) 


colors_31 = c("#640b0b", "#15226E", "#000000", "#324899", "#3B3B3B", "#525252", "#696969",
              "#08519c", "#a50f15", "#B53F3F", "#FF8800", "#B400FF", "#6F84F7", "#FFB700",
              "#75B2FF", "#CF8888", "#CCBC43", "#F50000", "#ABABAB", "#D1114D", "#0073FF",
              "#FFEA00", "#F558C6", "#0026FF", "#C4C4C4", "#520075", "#00EDFF", "#F500C4",
              "#3f301d","#795c34","#9a7b4f")

tidy_data$SpeciesNumber=as.numeric(factor(tidy_data$Category))

plot2 = ggplot(tidy_data, aes(x = SpeciesNumber, y = Value, fill = factor(SpeciesNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clusternumber, scales = "free_x", ncol = 3) +
  scale_x_continuous(
    breaks = 1:length(unique(tidy_data$Category)),  
    labels = 1:length(unique(tidy_data$Category))) +
  scale_fill_manual(
    values = colors_31, 
    labels = paste(1:length(unique(tidy_data$Category)), unique(tidy_data$Category), sep = ": ")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal()+
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24)) +
  labs(x = "Species (Numbered)",
    y = "Proprotion of EVR in Each Cluster",
    fill = "Species Key")
  
plot2
ggsave("Figure A2 - Proportion of EVR by Species for Each Cluster.jpg",
       plot = plot2, 
       width = 40, 
       height = 15, 
       units = "in", 
       dpi = 300)   


summarized_data_m=metiercasting %>%
  group_by(metier) %>% 
  summarise(across("ARROWTOOTH FLOUNDER":"YELLOWTAIL ROCKFISH", 
                   sum, na.rm = TRUE)) 

workm=sweep(summarized_data_m[, -1], 2, colSums(summarized_data_m[, -1]), FUN = "/")

workm$metiers = c("A", "B", "C", "D")

workm = relocate(workm, metiers)

tidy_data_m=workm %>%
  pivot_longer(cols = -metiers, names_to = "Category", values_to = "Value")


tidy_data_m$SpeciesNumber=as.numeric(factor(tidy_data_m$Category)) 



plot2 = ggplot(tidy_data_m, aes(x = SpeciesNumber, y = Value, fill = factor(SpeciesNumber))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~metiers, scales = "free_x", ncol = 2) +
  scale_x_continuous(breaks = 1:length(unique(tidy_data_m$Category)),  
    labels = 1:length(unique(tidy_data_m$Category)) ) +
  scale_fill_manual(values = colors_31, 
    labels = paste(1:length(unique(tidy_data_m$Category)), unique(tidy_data_m$Category), sep = ": ")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_minimal()+
  theme(legend.title = element_text(size = 30),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24)) +
  labs(x = "Species (Numbered)",
    y = "Proprotion of EVR in Each Metier",
    fill = "Species Key")

plot2

ggsave("Figure A7 - Proportion of EVR by Species for Each Metier.jpg",
       plot = plot2, 
       width = 35, 
       height = 15, 
       units = "in", 
       dpi = 300)   

#########Testing


data_long=casting %>%
  pivot_longer(
    cols = -c(clusternumber, shotID),
    names_to = "species",
    values_to = "catch"
  )

data_long=data_long %>%
  group_by(clusternumber) %>%
  mutate(proportion = catch / sum(catch))


stackbartest = ggplot(data_long, aes(x = factor(clusternumber), y = proportion, fill = species)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster Number",
    y = "Proportion of Total Catch In Cluster",
    fill = "Species Key") +
  scale_fill_manual(values = colors_31, 
    labels = paste(1:length(unique(tidy_data$Category)), unique(tidy_data$Category), sep = ": ")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(breaks = 1:17) +
  theme_minimal() + 
  theme(legend.title = element_text(size = 10),
                          strip.text = element_text(size = 30),
                          axis.text.y = element_text(size = 20),
                          axis.text.x = element_text(size = 24),    
                          axis.title.x = element_text(size = 24),     
                          axis.title.y = element_text(size = 24))



A4 = clusterproporbar / stackbartest + 
  plot_layout(heights = c(1, 3))

A4

ggsave("Figure A4 - Proportion of EVR by Speices for Each Cluster.jpg",  
       plot = A4, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)


data_longm=metiercasting %>%
  pivot_longer(
    cols = -c(metier, shotID),
    names_to = "species",
    values_to = "catch"
  )

data_longm=data_longm %>%
  group_by(metier) %>%
  mutate(proportion = catch / sum(catch))

stackbartestm = ggplot(data_longm, aes(x = factor(metier), y = proportion, fill = species)) +
  geom_bar(stat = "identity") +
  labs(x = "Metier",
    y = "Proportion of Total Catch In Metier",
    fill = "Species Key") +
  scale_fill_manual(values = colors_31, 
    labels = paste(1:length(unique(tidy_data$Category)), unique(tidy_data$Category), sep = ": ")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 10),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24))


F2 = metierproporbar / stackbartestm + 
  plot_layout(heights = c(1, 3))

F2

ggsave("Figure 2 - Proportion of EVR by Speices for Each Metier.jpg",  
       plot = F2, 
       width = 15,  
       height = 10, 
       units = "in", 
       dpi = 300)





sciname = unique(data$PACFIN_SPECIES_COMMON_NAME)
comname = unique(data$PACFIN_SPECIES_COMMON_NAME)


########### read in excel file separated ############

geardata = read_excel("2011-2023 Gears - 17 clusters.xlsx")


propor=geardata %>%
  group_by(clusternumber, Gear) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(propor)


gearbar = ggplot(propor, aes(x = factor(clusternumber), y = Proportion, fill = Gear)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster Number",
    y = "Proportion of Gear In Metier",
    fill = "Gear Key") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 10),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24))

A5 = clusterproporbar / gearbar + 
  plot_layout(heights = c(1, 3))

A5

ggsave("Figure A5- Proportion of Gear Used in Each Cluster.jpg",  
       plot = A5, 
       width = 15,  
       height = 10, 
       units = "in",
       dpi = 300)





metiergeardata = geardata

metiergeardata$metier=clusterstometiers[as.character(metiergeardata$clusternumber)]

metiergeardata = metiergeardata[c(1, 3:4)]


proporm=metiergeardata %>%
  group_by(metier, Gear) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))


gearbarm = ggplot(proporm, aes(x = factor(metier), y = Proportion, fill = Gear)) +
  geom_bar(stat = "identity") +
  labs(x = "Metier",
    y = "Proportion of Gear In Metier",
    fill = "Gear Key") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 10),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24))

F3 = metierproporbar / gearbarm + 
  plot_layout(heights = c(1, 3))

F3

ggsave("Figure 3- Proportion of Gear Used in Each Metier.jpg",  
       plot = F3, 
       width = 15,  # Width in inches
       height = 10, # Height in inches
       units = "in", # Units can be 'in', 'cm', or 'mm'
       dpi = 300)

#################
portdata = read_excel("cluster By Port for graph.xlsx")
names(portdata)


filterdata = data[c("PACFIN_PORT_NAME","COUNTY_STATE")]

unqiuedf=unique(filterdata)

dataframe1=merge(portdata, unqiuedf, by = "PACFIN_PORT_NAME", all.x = TRUE)

dataframe1 = dataframe1 %>% drop_na()

propor=dataframe1 %>%
  group_by(clusternumber, COUNTY_STATE) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(propor)


statebar = ggplot(propor, aes(x = factor(clusternumber), y = Proportion, fill = COUNTY_STATE)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster Number",
    y = "Proportion of Catch in State in Cluster",
    fill = "State Key") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24))


A6 = clusterproporbar / statebar + 
  plot_layout(heights = c(1, 3))

A6

ggsave("Figure A6- Proportion of Catch By State Each Cluster.jpg",  
       plot = A6, 
       width = 15,  
       height = 10,
       units = "in", 
       dpi = 300)


metierportdata = portdata

metierportdata$metier=clusterstometiers[as.character(metierportdata$clusternumber)]

metierportdata = metierportdata[c(1:2, 4)]

dataframe1m=merge(metierportdata, unqiuedf, by = "PACFIN_PORT_NAME", all.x = TRUE)

dataframe1m = dataframe1m %>% drop_na()

propor=dataframe1m %>%
  group_by(metier, COUNTY_STATE) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(propor)


statebarm = ggplot(propor, aes(x = factor(metier), y = Proportion, fill = COUNTY_STATE)) +
  geom_bar(stat = "identity") +
  labs(    x = "Metier",
    y = "Proportion of EVR in State in Metier",
    fill = "State") +
  scale_fill_viridis_d() +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 20),
                strip.text = element_text(size = 30),
                axis.text.y = element_text(size = 24),
                axis.text.x = element_text(size = 24),    
                axis.title.x = element_text(size = 24),     
                axis.title.y = element_text(size = 24))

F4 = metierproporbar / statebarm + 
  plot_layout(heights = c(1, 3))

F4

ggsave("Figure 4 - Proportion of EVR By State Each Metier.jpg",  
       plot = F4, 
       width = 15,  
       height = 10,
       units = "in", 
       dpi = 300)



silhouettescores = c(0, 0, 0, 0.3235494, 0.3908636, 0.3991081, 0.4037481, 
                       0.4207898, 0.4314634, 0.4316692, 0.4308888, 0.4283984, 
                       0.421593, 0.4365652, 0.436694, 0.4437081, 0.4555859, 
                       0.4449974, 0.4475046)
df = data.frame(nclusters = 1:length(silhouettescores),
  silhouette = silhouettescores)

silhouetteplot = ggplot(df, aes(x = nclusters, y = silhouette)) +
  geom_point(size = 4) +
  labs(x = "Number of Clusters",
    y = "Silhouette score") +
  theme_minimal()+ 
  theme(legend.title = element_text(size = 10),
        strip.text = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),    
        axis.title.x = element_text(size = 24),     
        axis.title.y = element_text(size = 24))


silhouetteplot


ggsave("Figure A1- Silhouette Classification Plot for CLARA Clustering.jpg",  
       plot = silhouetteplot, 
       width = 15,  
       height = 10,
       units = "in", 
       dpi = 300)

