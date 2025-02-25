library(MCMCglmm)
library(tidyverse)
library(lme4)

library(MCMCglmm)
library(MCMCpack)
library(tidyverse)
library(pedtricks)
library(brms)
library(AGHmatrix)
library(tidybayes)
source("Functions/Estimates from MCMCglmm.R")
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf

data_observed<-read.table("C:/Users/torhf/OneDrive - NTNU/Density dependent fitness variance/R/Density_dependent_fitness_variance/Cleaned data for density dependence analyses/cleaned_density_dependence_data.txt", sep=";", header=TRUE)
ped.hamish <-read.delim("C:/Users/torhf/Dropbox/data/New pedigree/helge_ped_err0.0027_adj_04-03-2024.txt", sep ="", stringsAsFactors = F)
morph <-read.delim("C:/Users/torhf/Dropbox/data/sql/Raw_Data/morphology.txt", sep =";", stringsAsFactors = F)
head(morph)

#From course notes:
#If the parent(s) of an individual are unknown then a missing value (NA)
#should be assigned in the relevant column. All individuals appearing as dams or
#sires need to have their own record, even if both of their parents are unknown.
#Often the number of individuals in a pedigree will be greater than the number of
#individuals for which phenotypic data exist. MCMCglmm can handle this, as long
#as all the individuals appearing in the data frame passed to data also appear in
#the pedigree.

#remove blood sample information from ID columns by keeping the 1st to 7th
#character in the ID. Should work with all versions of the pedigree.
ped.hamish<-ped.hamish %>% 
  mutate(id = substr(id, 1, 7)) %>% 
  mutate(sire = substr(sire, 1, 7)) %>% 
  mutate(dam = substr(dam, 1, 7)) %>% 
  dplyr::select(id, dam, sire)

#fix_ped:Returns a pedigree in which all individuals that exist 
#in the dam and sire columns are represented by their own record 
#lines, occurring before the records of their first offspring. If data 
#are supplied, then fix_ped will return a dataframe, the first three 
#columns are the 'fixed' pedigree, and the following columns of which 
#contain appropriately reordered data.
pedigree_data<-fix_ped(ped.hamish)
str(pedigree_data)

pedigree_data<-pedigree_data %>% 
  rename(animal = id)

data_observed<-data_observed %>% 
  mutate(Least_age = Least_age-1) %>% 
  mutate(Least_age_sq = Least_age^2)

data_observed<-data_observed %>% 
  mutate(Experienced_breeder = ifelse(Least_age == 0, 0, 1))

mean_expbreed<-mean(data_observed$Experienced_breeder)

data_observed<-data_observed %>% 
  mutate(mc_Experienced_breeder = Experienced_breeder-mean_expbreed)


#create an animal column
data_observed$animal <- data_observed$ID

data_wped <- data_observed
#This is probably better to do with the hestman dataframe
#data_wped<-data_wped %>% 
# filter(animal %in% pedigree_data$animal) #there are only 402 rows out of 10961 not in pedigree
#data_wped<-data_wped %>% 
#  left_join(pedigree_data)

data_wped<-data_wped %>% 
  left_join(pedigree_data)
nrow(data_wped)

data_wped<-data.frame(data_wped)
pedigree_data<-data.frame(pedigree_data)

str(pedigree_data)
class(data_wped)

#Hestmannøy####
#The poisson still gets the overdispersion with units

#retain only individuals in the pedigree. We only lose 45.
hestman<-data_wped %>% 
  filter(Location == "hestmannøy") %>% 
  filter(ID %in% pedigree_data$animal)

#prunePed: Creates a subset of a pedigree by retaining the 
#ancestors of a specified subset of individuals
#make.base=TRUE removes uninformative individuals and should have no impact
#on the posterior distributions. See ?prunePed
hestped<-prunePed(pedigree_data,hestman$animal, make.base = T)

write_delim(hestped, file = "Generated data/pedigree_hestmannoy.txt", delim =";")

sd(hestman$N) #we should do 50!

hestman<-hestman %>% 
  mutate(Nc = N/50)

hestman<-hestman %>% 
  mutate(Experienced_breeder = ifelse(Least_age == 0, 0, 1))

mean_expbreed<-mean(hestman$Experienced_breeder)

hestman<-hestman %>% 
  mutate(mc_Experienced_breeder = Experienced_breeder-mean_expbreed)

hestman<-data.frame(hestman)
hestped<-data.frame(hestped)


##uninformative improper flat ####
#From Hadfield course notes: Although inverse-Wishart distributions with negative degree of belief parameters
#are not defined, the resulting posterior distribution can be defined if there is
#sufficient replication. Specifying V=0 and n=-1 is equivalent to a uniform prior
#for the standard deviation on the the interval (0;1], and specifying V=0 and
#n=-2 is non-informative for a variance component.

#V= 1e-16
#nu = -2

nitt = 2130000
burnin = 30000
nt=2000

prior_flat <- list(
  G = list(
    G1 = list(V = diag(2) * 1e-16, nu = -2),  # For `us(1 + n):animal` (2x2 matrix)
    G2 = list(V = diag(2) * 1e-16, nu = -2),# For `us(1 + n):ID` (2x2 matrix)
    G3 = list(V = 1e-16, nu = -2)      # For `Year` (1x1 matrix)
  ),
  R = list(
    V = 1e-16, nu = -2                    # Residual variance (1x1 matrix)
  )
)

hestmod_Nc_flat <- MCMCglmm(
  fixed = r2s ~ Nc + mc_Experienced_breeder*sex,  
  random = ~ us(1 + Nc):animal +us(1 + Nc):ID + Year,
  rcov = ~ units,
  family = "poisson",
  pedigree = hestped,
  pr=TRUE, #to store BLUPs!
  prior = prior_flat, 
  data = hestman,
  nitt = nitt,
  burnin = burnin,
  thin = nt
)

saveRDS(hestmod_Nc_flat, "Workspace backup/hestmod_Nc_flat.rds")
hestmod_Nc_flat<-readRDS("Workspace backup/hestmod_Nc_flat.rds")
summary(hestmod_Nc_flat)
plot(hestmod_Nc_flat[["VCV"]])
MCMC2table(hestmod_Nc_flat)
autocorr.diag(hestmod_Nc_flat[["VCV"]])

varpost_uninf<-as.data.frame(hestmod_Nc_flat[["VCV"]])

median_qi(varpost_uninf$`(Intercept):(Intercept).animal`)
median_qi(varpost_uninf$`Nc:Nc.animal`)
median_qi(varpost_uninf$`Nc:(Intercept).animal`)
hist(varpost_uninf$`(Intercept):(Intercept).animal`, breaks = 100)
hist(varpost_uninf$`(Intercept):Nc.animal`, breaks = 100)
hist(varpost_uninf$`Nc:Nc.animal`, breaks = 100, col = "red",
     main="Variance distributions of random individual dens.reg function slopes",
     xlab = "Variance")


mycol2 <- rgb(0, 255, 0, max = 255, alpha = 100, names = "hej")
hist(varpost_uninf$`Nc:Nc.ID`, breaks = 100,col = mycol2,border=mycol2, add = T) 

#compare with the p model
hestmod_p_Nc<-readRDS("Workspace backup/hestmod_p_Nc.rds")
MCMC2table(hestmod_p_Nc, width = 0.95)
MCMC2table(hestmod_p_Nc, width = 0.95)[10,]
#the sum of the animal model slopes and permanent environment slopes should be the same as the phenotypic ID slopes
animal_ID_sum<-median(varpost_uninf$`Nc:Nc.ID` + varpost_uninf$`Nc:Nc.animal`)
#compare with the lme4 model
hestlme4.pois_bin_age <- readRDS("Workspace backup/hestlme4.pois_bin_age.rds")
VarCorr(hestlme4.pois_bin_age)$ID
VarCorr(hestlme4.pois_bin_age)$ID[4]

abline(v=(median(varpost_uninf$`Nc:Nc.animal`)), lwd=2, col = "red")
abline(v=(median(varpost_uninf$`Nc:Nc.ID`)), lwd=2, col = "green")
abline(v=animal_ID_sum, lwd=2, col = "purple")
ID_p_MCMC_slope<-MCMC2table(hestmod_p_Nc, width = 0.95)[10,2]
abline(v=ID_p_MCMC_slope, lwd=2, col = "orange")
abline(v=VarCorr(hestlme4.pois_bin_age)$ID[4], lwd=2, col = "yellow")

leg.txt <- c("animal", "permanent_env", "sum_animal_perm_env","ID_p_MCMC", "ID_lme4")
cols <- c("red", "green", "purple","orange", "yellow")
legend("topright", leg.txt, lty = 1,col = cols, lwd=2)

#same with intercept variance

hist(varpost_uninf$`(Intercept):(Intercept).animal`, breaks = 100, col = "red",
     main="Variance distributions of random individual dens.reg function intercepts",
     xlab = "Variance")

mycol2 <- rgb(0, 255, 0, max = 255, alpha = 100, names = "hej")
hist(varpost_uninf$`(Intercept):(Intercept).ID`, breaks = 100,col = mycol2,border=mycol2, add = T) 

#compare with the p model
hestmod_p_Nc<-readRDS("Workspace backup/hestmod_p_Nc.rds")
MCMC2table(hestmod_p_Nc, width = 0.95)
MCMC2table(hestmod_p_Nc, width = 0.95)[7,]

#the sum of the animal model intercepts and permanent environment intercepts should be the same as the phenotypic ID intercepts
animal_ID_sum_int<-median(varpost_uninf$`(Intercept):(Intercept).ID` + varpost_uninf$`(Intercept):(Intercept).animal`)

#compare with the lme4 model
hestlme4.pois_bin_age <- readRDS("Workspace backup/hestlme4.pois_bin_age.rds")
VarCorr(hestlme4.pois_bin_age)$ID
VarCorr(hestlme4.pois_bin_age)$ID[1]

abline(v=(median(varpost_uninf$`(Intercept):(Intercept).animal`)), lwd=2, col = "red")
abline(v=(median(varpost_uninf$`(Intercept):(Intercept).ID`)), lwd=2, col = "green")
abline(v=animal_ID_sum_int, lwd=2, col = "purple")
ID_p_MCMC_int<-MCMC2table(hestmod_p_Nc, width = 0.95)[7,2]
abline(v=ID_p_MCMC_int, lwd=2, col = "orange")
abline(v=VarCorr(hestlme4.pois_bin_age)$ID[1], lwd=2, col = "yellow")

leg.txt <- c("animal", "permanent_env", "sum_animal_perm_env","ID_p_MCMC", "ID_lme4")
cols <- c("red", "green", "purple","orange", "yellow")
legend("topright", leg.txt, lty = 1,col = cols, lwd=2)
