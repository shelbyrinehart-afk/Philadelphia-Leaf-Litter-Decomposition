### Data Analysis and Visualization for Leaf Litter Field Incubation ###
### Associated Manuscript: Heat vulnerability and leaf traits affect litter decomposition across an urban heat island in Philadelphia, Pennsylvania.   
### Incubation time frame: May 2024 - November 2024.
### Factors: Species, Site, (Mesh) - excluding the Mesh factor and focusing on small mesh samples only due to lack of balanced sampling design. 
###Independent variable: Mass loss. 
### Analysis conducted by: Dana Frankenstein, Shelby Rinehart ###

library(car) #levene's test, 
library(ggplot2)#plotting 
library(dplyr) #mutate
library(DescTools) #DunnTest
library(multcomp) #multcompLetters2
library(multcompView)
library(ggpmisc)#stat_poly_eq
library(lme4)#linear mixed models
library(lmerTest)#linear mixed models 
library(emmeans) #estimated marginal means for tea bags 


data1 <- read.csv("litter_data.csv")
head(data1)
table(data1$Species_Name)
table(data1$Site)
table(data1$Mesh)
table(data1$Common_Name, data1$Site, data1$Mesh)

#Determine decay rate based on exact deployment and retrieval dates and mass loss data. 

data2 <-mutate(data1, date_deployed = case_when(Site == "BA4" ~ "2024/05/02", Site == "CC1" ~ "2024/05/06", Site == "CG3" ~ "2024/05/07", Site == "W2"~ "2024/05/02"))
data3 <-mutate(data2, date_retrieved = "2024/11/05")
data3$diff_in_days<- difftime(data3$date_retrieved ,data3$date_deployed, units = c("days"), tz='EST')
data3$diff_in_days<-as.numeric(data3$diff_in_days)
data4<- mutate(data3, decay_rate_g_day = Total_wt_lost_g / diff_in_days)
data5<-mutate(data4, decay_rate_mg_day = decay_rate_g_day * 1000)

#subset the data to include only small mesh samples, since we do not have a balanced design if we include large as well (we are missing large sycamore bags at all sites).

data6<-subset(data5, Mesh == "Small")
table(data6$Common_Name, data6$Site)


### STATISTICAL ANALYSIS ###

#### Main Factors: Decay rate by Site and Species ####


#ANOVA analysis with Tukey post hoc. 

mod101<-aov(data=data6, decay_rate_mg_day ~ Site * Common_Name)
summary(mod101)
anova(mod101)

posthoc101a<-TukeyHSD(x=mod101, 'Common_Name', conf.level=0.95)
posthoc101b<-TukeyHSD(x=mod101, 'Site', conf.level=0.95)
posthoc101c<-TukeyHSD(x=mod101, 'Site:Common_Name', conf.level=0.95)


#letter differences for plotting purposes. 
cld101a <-multcompLetters2(decay_rate_mg_day ~ Common_Name, posthoc101a$Common_Name[,"p adj"], data6)
cld101b <-multcompLetters2(decay_rate_mg_day ~ Site, posthoc101b$Site[,"p adj"], data6)
cld101c <-multcompLetters2(decay_rate_mg_day ~ 'Site:Common_Name', posthoc101c$'Site:Common_Name'[,"p adj"], data6)

#check assumptions
shapiro.test(resid(mod101))
hist(resid(mod101))
leveneTest(mod101)
plot(mod101, which=1)
plot(mod101, which=2)
table(data6$Common_Name, data6$Site)


#Data is independent. Group sample sizes are nearly equal. Assumptions were investigated both with descriptive statistics and visual inspection to provide a holistic overview. Variance is heterogeneous and residuals are non-normal, but ANOVA should be robust enough to proceed. For added justification, see non-parametric Kruskal-wallis and DunnTest output below, which mirrors ANOVA and Tukey results. 

kruskal.test(data=data6, decay_rate_mg_day ~ Site)
kruskal.test(data=data6, decay_rate_mg_day ~ Common_Name)

DunnTest(decay_rate_mg_day ~ Site, data=data6, method='bonferroni')
DunnTest(decay_rate_mg_day ~ Common_Name, data=data6, method='bonferroni')

#### end ####

#### Regression: Temperature by Site ####

#Since we collected temperature data at each of our sites, we produce a linear regression between Site and temperature, to provide further justification for our site selection (based on heat-stress) and suggest reasons for the trends observed above. 

BA_sitedata<-read.csv("BA_sitedata.csv")
CC_sitedata<-read.csv("CC_sitedata.csv")
G_sitedata<-read.csv("G_sitedata.csv")
W_sitedata<-read.csv("W_sitedata.csv")

#calculate mean temperature at each site for the entire incubation. 

W_sitedata2<-W_sitedata %>% summarize (temp_mean = mean(temp))
W_sitedata3<-mutate(W_sitedata2, Site = "W2")

G_sitedata2<-G_sitedata %>% summarize (temp_mean = mean(temp))
G_sitedata3<-mutate(G_sitedata2, Site = "CG3")

CC_sitedata2<-CC_sitedata %>% summarize (temp_mean = mean(temp))
CC_sitedata3<-mutate(CC_sitedata2, Site = "CC1")

BA_sitedata2<-BA_sitedata %>% summarize (temp_mean = mean(temp))
BA_sitedata3<-mutate(BA_sitedata2, Site = "BA4")

site_data<-rbind(BA_sitedata3, CC_sitedata3)
site_data2<-rbind(G_sitedata3, W_sitedata3)
site_data3<-rbind(site_data, site_data2)

#merge with data frame that shows mean decay rate at each site. 

data7<-data6 %>% group_by(Site) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       stdev_decay_rate_mg_day = sd(decay_rate_mg_day))

site_data4<-merge(data7, site_data3, by =c("Site"))

cor(site_data4$mean_decay_rate_mg_day, site_data4$temp_mean) #strong negative correlation. 



#Binning sites by less or more urbanized to later calculate percent difference in mean decay.
data7.1<-mutate(data6, urban_grp = case_when(Site == "CC1" ~ "less urbanized", Site == "W2" ~ "less urbanized", Site == "CG3" ~ 'more urbanized', Site == "BA4" ~ "more urbanized"))

data7.2<-data7.1 %>% group_by(urban_grp) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       stdev_decay_rate_mg_day = sd(decay_rate_mg_day))

#### end ####

##### Carbon Source: labile vs. recalcitrant #### 

data8<-mutate(data6, C_source = case_when(Common_Name == "Gingko" ~ "L", Common_Name == "Sweetgum" ~ "L", Common_Name == "Red Oak" ~ 'R', Common_Name == "Sycamore" ~ "R"))

#mean by carbon source. 
data8.1<-data8 %>% group_by(C_source) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))

#mean by carbon source and site. 
data8.2 <-data8 %>% group_by(C_source, Site) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))


mod102<-aov(data=data8, decay_rate_mg_day ~ C_source *Site)
anova(mod102) #both C source and site are significant. 

leveneTest(mod102) #significant.
shapiro.test(resid(mod102)) #Not normal. 
table(data8$Site, data8$C_source)
plot(mod102, which =2)
hist(resid(mod102))

#anova is robust to these slight deviations and visual inspections largely pass. 

posthoc102a<-TukeyHSD(x=mod102, 'Site', conf.level=0.95) #confirming our previous results. 
posthoc102b<-TukeyHSD(x=mod102, "C_source:Site", conf.level=0.95)

#non-parametric test for added clarity and corroboration of our findings. 
kruskal.test(data=data8, decay_rate_mg_day ~ C_source)
DunnTest(decay_rate_mg_day ~ C_source, data=data8, method='bonferroni')

#### end ####

#### Stoichiometry ####

###Organic Matter ###

#INITIAL
om_initial <-read.csv('om_initial.csv')
table(om_initial$species)
om_initial1<-om_initial%>% group_by(species) %>% summarize (mean_om = mean(om_perc),                                                       sd_om= sd(om_perc))


#FINAL

om_final<-read.csv('om_final.csv')
#remove duplicate sample 
om_final1<-om_final[-c(13),]
#remove large mesh 
om_final2<-subset(om_final1, mesh == "Small")
table(om_final2$site, om_final2$species) #good 
class(om_final2$species)
om_final2$site<-factor(om_final2$site, levels=c("CC","W","G",'BA'))
om_final2$species<-factor(om_final2$species, levels=c("Gingko","Sweetgum","Red Oak ","Sycamore"))



#COMPARISON 

om_initial2<-mutate(om_initial, time='initial')
om_initial3<-dplyr::select(om_initial2, "species",'time','om_perc')

om_final3<-mutate(om_final2, time='final')
om_final4<-mutate(om_final3, om_perc = om_percent)
om_final5<-dplyr::select(om_final4, 'species', 'time','om_perc')

om1<-rbind(om_initial3, om_final5)
table(om1$species, om1$time)

om_gk<-subset(om1, species == "Gingko")
om_sg<-subset(om1, species == "Sweetgum")
om_ro<-subset(om1, species == 'Red Oak ')
om_sy<-subset(om1, species == "Sycamore")


wilcox.test(om_perc ~ time, data=om_gk)#significant 
wilcox.test(om_perc ~ time, data=om_sg)#significant 
wilcox.test(om_perc ~ time, data=om_ro)#significant 
wilcox.test(om_perc ~ time, data=om_sy)#significant 


#Next, calculate percent difference between final and initial and then do anova on those values by species. 

om_final6<-mutate(om_final5, om_perc_fin =om_perc)
om_final6.1<-om_final6[,-c(2,3)]


om2<-merge(om_initial1, om_final6.1, by='species')
om3<-mutate(om2, perc_diff = ((abs(mean_om-om_perc_fin))/((mean_om+om_perc_fin)/2))*100
)

mod329<-aov(perc_diff ~ species, data=om3)
anova(mod329) #n.s.
summary(mod329)

shapiro.test(resid(mod329))
leveneTest(mod329)
hist(resid(mod329))
plot(mod329)
#no significant differences in percent difference in organic matter by species. 


### C:N ###

###INITIAL###

cn_initial <- read.csv('cn_initial.csv')
cn_initial1<-mutate(cn_initial, cn_ratio =percC / percN)
#removing samples with 0 nitrogen and thus NA for C:N ratio (likely instrument error). 
cn_initial2 <- cn_initial1[-c(3,4,13), ]
table(cn_initial2$species)
cn_initial3 <-cn_initial2%>% group_by(species) %>% summarize (mean_cn_in = mean(cn_ratio),                                                       sd_cn_in= sd(cn_ratio))

#cn_initial2$species<-factor(cn_initial2$species, levels=c('Gingko ', 'Sweetgum',"Red Oak","Sycamore"))


### FINAL ###
cn_final<-read.csv("cn_final.csv")
table(cn_final$species, cn_final$site, cn_final$mesh)
cn_final1<-subset(cn_final, mesh == "Small" )

cn_final2 <-cn_final1%>% group_by(species) %>% summarize (mean_cn_fin = mean(cn_ratio),                                                       sd_cn_fin= sd(cn_ratio))

#CN1$site<-factor(CN1$site, levels=c("CC","W","G","BA"))
#CN1$species <- factor(CN1$species, levels = c("Gingko","Sweetgum","Red Oak","Sycamore"))

### COMPARISON ### 
cn_initial4<-mutate(cn_initial2, time = "initial")
cn_initial5<-dplyr::select(cn_initial4, "species",'time','cn_ratio')


cn_final3<-mutate(cn_final1, time='final')
cn_final4<-dplyr::select(cn_final3, 'species', 'time','cn_ratio')

cn<-rbind(cn_initial5, cn_final4)
table(cn$species,cn$time)

#Group sizes are highly unequal, which is why we choose a Wilcoxon rank sum test (i.e Mann-whitney U) for each species, to compare differences between initial and final C:N. 

cn_gk<-subset(cn, species == "Gingko")
cn_sg<-subset(cn, species == "Sweetgum")
cn_ro<-subset(cn, species == "Red Oak")
cn_sy<-subset(cn, species == "Sycamore")

wilcox.test(cn_ratio ~ time, data=cn_gk) #significant
wilcox.test(cn_ratio ~ time, data=cn_sg)#not
wilcox.test(cn_ratio ~ time, data=cn_ro)#significant 
wilcox.test(cn_ratio ~ time, data=cn_sy)#not 

#Now, determine percent differences between final and (mean) initial for each species and then do anova. 

cn1<-merge(cn_final1, cn_initial3, by="species")
cn2<-mutate(cn1, perc_diff = ((abs(cn_ratio-mean_cn_in))/((cn_ratio+mean_cn_in)/2))*100
)

mod328<-aov(perc_diff ~ species, data=cn2)
anova(mod328)
summary(mod328)
TukeyHSD(x=mod328, 'species', conf.level=0.95)
shapiro.test(resid(mod328))
leveneTest(mod328)
hist(resid(mod328))
plot(mod328)

### Percent C ### 

cn_initial_c<-cn_initial2%>% group_by(species) %>% summarize (mean_percentc = mean(percC),                                                       sd_percentc= sd(percC))

cn_final_c<-cn_final1%>%group_by(species) %>% summarize (mean_percentc = mean(perc_c), sd_percentc=sd(perc_c))

cn_c<-merge(cn_initial_c, cn_final1, by='species')

cn_c1<-mutate(cn_c, perc_diff = ((abs(mean_percentc-perc_c))/((mean_percentc+perc_c)/2))*100
)

mod330<-aov(perc_diff ~ species, data=cn_c1)
anova(mod330) #significant
summary(mod330)
TukeyHSD(x=mod330, 'species', conf.level=0.95)

shapiro.test(resid(mod330))
leveneTest(mod330)
hist(resid(mod330))
plot(mod330)


cn_initial9<-dplyr::select(cn_initial2, "species", "percC")
cn_initial9.1<-mutate(cn_initial9, time = 'initial')

cn_final8<-dplyr::select(cn_final1, "species", "perc_c")
cn_final8.1<-mutate(cn_final8, time = "final")
cn_final8.2<-mutate(cn_final8.1, percC = perc_c)
cn_final8.3 <-cn_final8.2[,-c(2)]

cn4<-rbind(cn_initial9.1, cn_final8.3)

table(cn4$species, cn4$time)

cn_gk_c<-subset(cn4, species == "Gingko")
cn_sg_c<-subset(cn4, species == "Sweetgum")
cn_ro_c<-subset(cn4,  species == 'Red Oak')
cn_sy_c<-subset(cn4, species == "Sycamore")


wilcox.test(percC ~ time, data=cn_gk_c)#significant 
wilcox.test(percC ~ time, data=cn_sg_c)#significant 
wilcox.test(percC ~ time, data=cn_ro_c)#n.s.
wilcox.test(percC ~ time, data=cn_sy_c)#significant 

### Percent N ### 

cn_initial_n<-cn_initial2%>% group_by(species) %>% summarize (mean_percentn = mean(percN),                                                       sd_percentn= sd(percN))

cn_final_n<-cn_final1%>%group_by(species) %>% summarize (mean_percentn = mean(perc_n), sd_percentn=sd(perc_n))

cn6<-merge(cn_initial_n, cn_final1, by='species')

cn6.1<-mutate(cn6, perc_diff = ((abs(mean_percentn-perc_n))/((mean_percentn+perc_n)/2))*100
)

mod331<-aov(perc_diff ~ species, data=cn6.1)
anova(mod331) #significant
summary(mod331)
TukeyHSD(x=mod331, 'species', conf.level=0.95)

shapiro.test(resid(mod331))
leveneTest(mod331)
hist(resid(mod331))
plot(mod331)


cn_initial9.2<-dplyr::select(cn_initial2, "species", "percN")
cn_initial9.3<-mutate(cn_initial9.2, time = 'initial')

cn_final9<-dplyr::select(cn_final1, "species", "perc_n")
cn_final9.1<-mutate(cn_final9, time = "final")
cn_final9.2<-mutate(cn_final9.1, percN = perc_n)
cn_final9.3 <-cn_final9.2[,-c(2)]

cn7<-rbind(cn_initial9.3, cn_final9.3)

table(cn7$species, cn7$time)

cn_gk_n<-subset(cn7, species == "Gingko")
cn_sg_n<-subset(cn7, species == "Sweetgum")
cn_ro_n<-subset(cn7,  species == 'Red Oak')
cn_sy_n<-subset(cn7, species == "Sycamore")


wilcox.test(percN ~ time, data=cn_gk_n)#significant 
wilcox.test(percN ~ time, data=cn_sg_n)#n.s.
wilcox.test(percN ~ time, data=cn_ro_n)#significant
wilcox.test(percN ~ time, data=cn_sy_n)#n.s. 

#### end ####

#### Tea Bags #### 

## Tea bag samples were deployed at the same time as natural leaf litter and collected in 2-month intervals to establish a time series of decay by green and rooibos tea, in accordance with the Tea Bag Index (Keuskamp et al., 2013). 9 individual samples were lost during the incubation and thus do not have data (8 rooibos compartments, 1 green tea compartment). Most samples were lost at the Garden site. 

teabags<-read.csv('teabag_data.csv')
teabags2<-teabags[complete.cases(teabags), ] #removed samples with missing data. 

table(teabags2$Site, teabags2$duration_months, teabags2$tea_type)
(length(which(teabags2$duration_months == "6")))

#calculate decay rate based on days incubated and mass loss. 

teabags3 <-mutate(teabags2, date_deployed = case_when(Site == "BA4" ~ "2024/05/02", Site == "CC1" ~ "2024/05/06", Site == "CG3" ~ "2024/05/07", Site == "W2"~ "2024/05/02"))

teabags4<-mutate(teabags3, date_collected = case_when(duration_months == "2" ~ "2024/07/01", duration_months == "4" ~ "2024/08/29", duration_months == "6" ~ "2024/11/01"))

teabags4$diff_in_days<- difftime(teabags4$date_collected ,teabags4$date_deployed, units = c("days"), tz="EST")

teabags4$diff_in_days<-as.numeric(teabags4$diff_in_days)

# Total mass loss values need to be adjusted for samples collected at the 6-month mark. While all previous samples were weighed without the actual tea bag, these were weighed still contained in the bag. Thus, we need to subtract the average bag weight (accepted values) from the mass loss values. Bags collected at 2 and 4 months do not need a correction. 



teabags4.5<-subset(teabags4, duration_months == "6")
teabags4.5$total_mass_loss<-as.numeric(teabags4.5$total_mass_loss)
teabags4.6 <-dplyr::mutate(teabags4.5, total_mass_loss_cor = case_when(tea_type == "Green" ~ total_mass_loss -0.246, tea_type == "Rooibos" ~ total_mass_loss - 0.245 ))

teabags4.7 <- subset(teabags4, duration_months == c("2"))
teabags4.8<-mutate(teabags4.7, total_mass_loss_cor = total_mass_loss)

teabags4.9 <- subset(teabags4, duration_months == c("4"))
teabags4.95<-mutate(teabags4.9, total_mass_loss_cor = total_mass_loss)

teabags4.96<-rbind(teabags4.95, teabags4.8)
teabags4.97<-rbind(teabags4.96, teabags4.6)
teabags4.97$total_mass_loss_cor<-as.numeric(teabags4.97$total_mass_loss_cor)

teabags4.98<-mutate(teabags4.97, decay_rate_g_day = total_mass_loss_cor / diff_in_days)
teabags5<-mutate(teabags4.98, decay_rate_mg_day = decay_rate_g_day * 1000)

#Summarize some mean values:

teabags6<-teabags5 %>% group_by(tea_type, duration_months, Site) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))


##Subset data to include 6-month collection only, to allow for later comparison with natural leaf litter also collected at 6 months. 

teabags7<-subset(teabags5, duration_months == "6")

#summarize decay by site and tea type: 
teabags7.1<- teabags7 %>% group_by(tea_type, Site) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))

#summarize by tea type only:
teabags7.2 <-teabags7 %>% group_by(tea_type) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))



#With the teabag data, we have a continuous response variable that has normally distributed error, but we have "clusters" or a grouped structure with the time factor. Thus, we want to model that not as a fixed effect but as a random effect. We are not interested in the differences imposed by the time variable, but want to account for it. We would use a linear mixed model for this. If our data is also not normal, then we would use a generalized linear mixed model. 

table(teabags5$Site, teabags5$tea_type, teabags5$duration_months)

#linear mixed model solution:
mod008<-lmer(decay_rate_mg_day ~ tea_type * Site + (1|duration_months), data=teabags5, REML = F)
summary(mod008)
car::Anova(mod008, type=2)
#tea type is the only significant factor. 

#reduced model:
model_reduced <- lmer(decay_rate_mg_day ~ tea_type + (1 | duration_months), data = teabags5, REML = F)
summary(model_reduced)
#Anova type 2 does wald test for hypothesis testing
Anova(model_reduced, type=2)
emmeans(model_reduced,pairwise~tea_type)

#likelihood ratio test to see if site inclusion improves model fit - since we are changing fixed effects and comparisons between models with different fixed effects needs to have REML = F. 
anova(mod008, model_reduced) #including site does not improve model fit. 


#likelihood ratio test on the random effect of time 
mod0081<-lmer(decay_rate_mg_day ~ tea_type + (1|duration_months), data=teabags5)
mod_reduced1<-lm(decay_rate_mg_day ~ tea_type, data=teabags5 )
anova(mod0081, mod_reduced1) #time significantly improves model fit. 


#check assumptions of anova. 
shapiro.test(resid(mod008)) #good 
leveneTest(residuals(mod008) ~ teabags5$Site) #good
leveneTest(residuals(mod008) ~ teabags5$tea_type) 
hist(resid(mod008))



#alternative solution with time as a fixed factor: 
mod009<-lm(decay_rate_mg_day ~ tea_type * Site * factor(duration_months), data=teabags5)
car::Anova(mod009, type=2)
summary(mod009)
#site effects show up marginally in this approach.

#alternative approach doing a glm on 6-month data only.
mod009b <- glm(decay_rate_mg_day ~ tea_type * Site, data = teabags7, family = Gamma(link="log"))
plot(fitted(mod009b), residuals(mod009b, type="deviance"))
anova(mod009b, test = "Chisq")
#still no real site effects. 





# Since we handled time as a random factor previously, we will fit a new model here to help with visualization purposes. We will create an exponential decay line for each tea type over time, disregarding the site factor (since it is non-significant anyways).

#fit separate models by tea type to get distinct r^2, decay lines,  and k-values (for visualization purposes on our figure). 

teabags5_green<-subset(teabags5, tea_type == "Green")

fit_green = glm(decay_rate_mg_day ~ duration_months, family=Gamma(link='log'), data=teabags5_green)
summary(fit_green)#k for green: ~0.28 (-(slope))

#determine predicted values and confidence intervals. 
#green tea

pred_green <- predict(fit_green, type = "link", se.fit = TRUE)
crit <- qnorm(0.975)

teabags5_green$pred    <- exp(pred_green$fit)
teabags5_green$CI_lwr  <- exp(pred_green$fit - crit * pred_green$se.fit)
teabags5_green$CI_upr  <- exp(pred_green$fit + crit * pred_green$se.fit)

teabags5_greenplot<-subset(teabags5_green, !duplicated(teabags5_green$pred))

#rooibos tea.

teabags5_rooibos<-subset(teabags5, tea_type == "Rooibos")

fit_rooibos = glm(decay_rate_mg_day ~ duration_months, family=Gamma(link='log'), data=teabags5_rooibos)
summary(fit_rooibos)
#k for rooibos: ~0.21,

pred_rooibos <- predict(fit_rooibos, type = "link", se.fit = TRUE)

teabags5_rooibos$pred    <- exp(pred_rooibos$fit)
teabags5_rooibos$CI_lwr  <- exp(pred_rooibos$fit - crit * pred_rooibos$se.fit)
teabags5_rooibos$CI_upr  <- exp(pred_rooibos$fit + crit * pred_rooibos$se.fit)

teabags5_rooibosplot<-subset(teabags5_rooibos, !duplicated(teabags5_rooibos$pred))


#find deviance-based pseudo r squared:  
pesudo_r2_rooibos<- 1-fit_rooibos$deviance / fit_rooibos$null.deviance
#~0.74
pseudo_r2_green <- 1-fit_green$deviance / fit_green$null.deviance
#~0.96 

#glm model fit was estimated using deviance-based pseudo r-squared, which showed a strong fit for both models (but especially fo green). 





#### end ####

#### Comparison: Tea Bag and Natural Leaf Litter Decay #### 

#Merge our data frames in an elegant way by making a new column for both data frames so we can use that as the x axis variable. 

tb <-mutate(data8.1, treatment = case_when(C_source == "L" ~ "Labile Carbon Species", C_source == "R" ~ "Recalcitrant Carbon Species"))
tb1<-tb[,-c(1)]

teabags7.3 <-mutate(teabags7.2, treatment = case_when(tea_type =="Green" ~ "Green teabags", tea_type == "Rooibos" ~ "Rooibos teabags"))
teabags7.4<-teabags7.3[,-c(1)]

tb2<-rbind(teabags7.4, tb1) #showing means and standard deviation of each treatment, useful for plotting later
tb2$treatment<-factor(tb2$treatment, levels= c("Green teabags", "Labile Carbon Species", "Rooibos teabags","Recalcitrant Carbon Species"))

#Analysis: combine data8 (natural leaf litter with carbon source column) with teabags 7 (tea bags collected at 6 months)


data8.3 <-subset(data8, select = c("Code", "Site", "decay_rate_mg_day","C_source"))
data8.4<-mutate(data8.3, treatment = case_when(C_source == "L" ~ "Labile Carbon", C_source == "R" ~ "Recalcitrant Carbon"))  
data8.5 <-data8.4[,-c(4)]


teabags7.5<-subset(teabags7, select = c("Code", "Site", "decay_rate_mg_day","tea_type"))
teabags7.6<-mutate(teabags7.5, treatment = case_when(tea_type == "Green" ~ "Green teabags", tea_type == "Rooibos" ~ "Rooibos teabags"))            
teabags7.7<-teabags7.6[,-c(4)]                


tb3<-rbind(teabags7.7, data8.5) 
table(tb3$treatment, tb3$Site)


#For the analysis, we will do a Kruskal-Wallis and Dunn's post-hoc. Group sizes are unequal, heteroscedastic, and residuals are non-normal. 

kruskal.test(data=tb3, decay_rate_mg_day ~ treatment)
DunnTest(decay_rate_mg_day ~ treatment, data=tb3, method='bonferroni')



#### end ####

#### Coefficient of Variation ####
# CV = standard deviation / mean (decay rate)

#for teabags
table(teabags5$Site, teabags5$duration_months, teabags5$tea_type)
#we only have one observation for CG at 4 months, and 2 for G site at 6 months, all others have at least 3 observations. 

teabags8<-teabags5%>% group_by(Site, duration_months, tea_type) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       stdev_decay_rate_mg_day = sd(decay_rate_mg_day))
teabags8.1<-mutate(teabags8, CV = stdev_decay_rate_mg_day/ mean_decay_rate_mg_day)

#for natural leaf litter
data9<-data6%>% group_by(Site,Common_Name) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       stdev_decay_rate_mg_day = sd(decay_rate_mg_day))

data9.1<-mutate(data9, CV = stdev_decay_rate_mg_day/ mean_decay_rate_mg_day)


mod046<-lm(data=data9.1, CV ~ Site)

anova(mod046)#almost significant. 
summary(mod046)
hist(resid(mod046))
shapiro.test(resid(mod046))#good
LeveneTest(mod046) #good 

#### end ####



### FIGURES ###

#### Figure 1: Site * Species Matrix Plot ####

data6$Site<-factor(data6$Site, levels = c("CC1","W2","CG3","BA4"))
data6$Common_Name<-factor(data6$Common_Name, levels = c("Gingko","Sweetgum","Red Oak","Sycamore"))

letters101c<-tibble(
  Common_Name=c("Gingko","Gingko","Gingko","Gingko","Sweetgum","Sweetgum","Sweetgum","Sweetgum","Red Oak","Red Oak","Red Oak","Red Oak","Sycamore","Sycamore","Sycamore","Sycamore"),
  Site = c("CC1","W2","CG3","BA4","CC1","W2","CG3","BA4","CC1","W2","CG3","BA4","CC1","W2","CG3","BA4"),
  x = c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4),
  y =  c(10.7,9.7,7.3,8.8, 9.8, 8.5, 7,8.5,7,15,4.3,5, 6.3,6.7, 5,4), 
  label = c('a','ab','abcd','abcde','ab','abc','abcde','bcde','abcde','abcd','de','de','abcde','abcde','cde','e')
)

ggplot(data=data6, aes(x=Site, y=decay_rate_mg_day, color=Site, fill=Site))+theme_bw()+geom_point( alpha=0.95, show.legend=FALSE)+geom_boxplot(alpha=0.5, outlier.color='black', outlier.shape=16, outlier.size = 2, show.legend=FALSE)+facet_wrap(~factor(Common_Name, levels=c('Gingko','Sweetgum','Red Oak','Sycamore')), ncol=4, scales="free_x")+scale_color_manual(breaks=c("CC1","W2","CG3","BA4"), values=c('seagreen3','deepskyblue2',"cornflowerblue",'skyblue4'))+scale_fill_manual(breaks=c("CC1","W2","CG3","BA4"), values=c('seagreen3','deepskyblue2',"cornflowerblue",'skyblue4'))+ scale_x_discrete(labels=c("CC1" = "CC", "W2" = "W","CG3"="G","BA4"="BA"))+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(),strip.background = element_rect(fill="white"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0), axis.title.y=element_text(size=14), axis.text.x=element_text(size=13),axis.text.y=element_text(size=13)) +geom_text(data=letters101c, aes(x=x, y=y, label=label), size=3, fontface='bold', inherit.aes=FALSE)+labs(y="Decay Rate (mg" ~day^-1*")", x="")

#### end ####

#### Figure 2a: Site only Boxplot ####

letters101b<-tibble(
  Site = c ("CC1","W2","CG3","BA4"), 
  x = c(1,2,3,4),
  y =  c(12,12,9,9), 
  label = c('a','a','b','b')
)

ggplot(data=data6, aes(x=Site, y=decay_rate_mg_day, color=Site, fill=Site))+theme_bw()+geom_point(position='jitter', alpha=0.9, show.legend=FALSE, aes(x=Site, y=decay_rate_mg_day, fill=Site))+geom_boxplot(alpha=0.6, color='black',outlier.color='black', outlier.shape=16, outlier.size = 2, show.legend=FALSE)+scale_color_manual(breaks=c("CC1","W2","CG3","BA4"), values=c('seagreen3','deepskyblue2',"cornflowerblue",'skyblue4'))+scale_fill_manual(breaks=c("CC1","W2","CG3","BA4"), values=c('seagreen3','deepskyblue2',"cornflowerblue",'skyblue4'))+ geom_text(data=letters101b, aes(x=x, y=y, label=label), size=4, fontface='bold', inherit.aes=FALSE)+labs(x="", y="Decay Rate (mg" ~day^-1*")", title= "")+
  theme(axis.title.y=element_text(size=14), axis.text.x=element_text(size=13),axis.text.y=element_text(size=13),panel.grid.major=element_blank(), panel.grid.minor=element_blank() )+ scale_x_discrete(labels=c("CC1" = "CC", "W2" = "W","CG3"="G","BA4"="BA"))


#### end ####

#### Figure 2b: Site / Temp Linear Regression ####

ggplot(site_data4, aes(x=temp_mean, y=mean_decay_rate_mg_day))+theme_bw()+labs(x="Mean Temperature (°C)", y="Mean Decay Rate (mg" ~day^-1*")", title= "")+ 
  theme(axis.text.x=element_text(size=13),axis.title.x = element_text(size=14), axis.text.y=element_text(size=13), axis.title.y=element_text(size=14), legend.text=element_text(size=13), legend.title=element_text(size=14), legend.position = c(0.7,0.9),legend.direction='horizontal', legend.box.background=element_rect(),legend.box.margin=margin(5,5,5,5), panel.grid.major=element_blank(), panel.grid.minor=element_blank())+stat_poly_eq(use_label(c("eq", "R2")), size=5, label.y=0.17)+
  stat_poly_line(color='black')+ geom_point(aes(fill=Site, shape=Site),size=(9), color='black')+geom_errorbar(aes(ymin = mean_decay_rate_mg_day - stdev_decay_rate_mg_day, ymax = mean_decay_rate_mg_day+ stdev_decay_rate_mg_day), width=0.1)+
  scale_shape_manual(name = 'Site', values=c(21,22,24,23),labels=c("CC","W", "G","BA"), breaks=c("CC1","W2", "CG3","BA4"))+scale_fill_manual(breaks=c("CC1","W2","CG3","BA4"),values=c('seagreen3','deepskyblue2',"cornflowerblue",'skyblue4'), labels=c("CC","W", "G","BA"), name="Site")

#### end ####

#### Figure 3a: Species only Boxplot ####

letters101a<-tibble(
  C_source = c("L","L","R","R"),
  Site = c ("Gingko ","Sweetgum","Red Oak","Sycamore"), 
  x = c(1,2,1,2),
  y =  c(12,12,9,9), 
  label = c('a','a','b','b')
)



ggplot(data=data8, aes(x=Common_Name, y=decay_rate_mg_day, color=Common_Name, fill=Common_Name))+theme_bw()+geom_point(position='jitter', alpha=0.9, show.legend=FALSE)+geom_boxplot(alpha=0.6, color='black',outlier.color='black', outlier.shape=16, outlier.size = 2, show.legend=FALSE)+scale_color_manual(breaks=c("Gingko","Sweetgum","Red Oak","Sycamore"), values=c("pink2",'orange2','maroon2','brown2'))+scale_fill_manual(breaks=c("Gingko","Sweetgum","Red Oak","Sycamore"), values=c("pink2",'orange2','maroon2','brown2'))+ geom_text(data=letters101a, aes(x=x, y=y, label=label), size=4, fontface='bold', inherit.aes=FALSE)+labs(x="", y="Decay Rate (mg" ~day^-1*")")+
  theme(axis.text.x=element_text(size=13),axis.title.x = element_text(size=14), axis.text.y=element_text(size=13), axis.title.y=element_text(size=14),strip.background = element_rect(fill="white", size=1, color="black"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0), panel.grid.major=element_blank(),panel.grid.minor=element_blank())+facet_wrap(vars(C_source), scales='free_x', labeller = labeller( C_source=   c("L" = "Labile C", "R" = "Recalcitrant C")))

#### end ####

#### Figure 3b: Decay rate by Site & Carbon Quality ####


data8.2$C_source<-factor(data8.2$C_source, levels = c("L","R"))
data8.2$Site<-factor(data8.2$Site, levels = c("CC1","W2","CG3","BA4"))

legend <- data8.2 %>% distinct(Site, C_source) %>% 
  mutate(
    fill = case_when(
      Site == "CC1" ~ 'seagreen3',
      Site == "W2" ~ 'deepskyblue2',
      Site == "CG3" ~ 'cornflowerblue',
      Site == "BA4"~ 'skyblue4',
    ),
    shape = case_when(
      C_source == "L"~ 21,
      C_source == "R"~ 22
    ),
    label = paste0(Site, C_source), 
  )

legend1<-mutate(legend, label=case_when(
  label == "CC1L"~ "CC + LC",
  label == "CC1R"~ "CC + RC",
  label == "W2L"~ "W + LC",
  label == "W2R"~ "W + RC",
  label == "CG3L"~ "G + LC",
  label == "CG3R"~ "G + RC",
  label == "BA4L"~ "BA + LC",
  label == "BA4R"~ "BA + RC"
  
))



ggplot(data8.2, aes(x=Site, y=mean_decay_rate_mg_day,fill = interaction(as.factor(Site), as.factor(C_source)),
                             shape = interaction(as.factor(Site), as.factor(C_source))))+theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_text(size=13),strip.background = element_rect(fill="white", size=1, color="black"), legend.text=element_text(size=13), legend.title=element_text(size=14), legend.position ='top',legend.direction='horizontal' ,legend.box.background=element_rect(color='black', size = 1),axis.text.y=element_text(size=13),axis.title.y=element_text(size=14))+
  labs(y="Mean Decay Rate (mg" ~day^-1*")",x="")+
  scale_shape_manual(
    name = "Site &
    C Quality", 
    labels = legend1 %>% pull(label),
    values = setNames(legend1$shape, interaction(legend1$Site, legend1$C_source))
  )+scale_fill_manual(
    name = "Site &
    C Quality", 
    labels = legend1 %>% pull(label),
    values = setNames(legend1$fill, interaction(legend1$Site, legend1$C_source))
  ) +
  geom_point(size=5, color='black')+
  scale_x_discrete(labels=c("CC1" = "CC", "W2" = "W","CG3"="G","BA4"="BA"))+
  geom_errorbar(aes(ymin = mean_decay_rate_mg_day - sd_decay_rate_mg_day, ymax = mean_decay_rate_mg_day+ sd_decay_rate_mg_day), width=0.05)

#### end ####

#### Figure 4a: Organic Matter####

om_final7 <-om_final4%>% group_by(species) %>% summarize (mean_om = mean(om_percent),                                                       sd_om= sd(om_percent))
om_final7.1<-mutate(om_final7, time = 'final')
om_initial4<-mutate(om_initial1, time = 'initial')
om8<-rbind(om_final7.1, om_initial4)

om8$time<-factor(om8$time, levels=c("initial",'final'))

letters001<-tibble(
  species=c("Gingko","Sweetgum","Red Oak ",'Sycamore'),
  x=c(1.5,1.5,1.5,1.5),
  y=c(93,96,96,96),
  label=c('***','***','***','***')
)


ggplot(data=om8, aes(y=mean_om, x=time, fill=time))+geom_col(show.legend=F, color='black', alpha=0.8)+facet_wrap(~factor(species, levels=c("Gingko","Sweetgum","Red Oak ",'Sycamore')), ncol=4)+geom_errorbar(aes(ymin = mean_om - sd_om, ymax = mean_om+sd_om), width=0.2, color='black')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=13,color='black'), axis.title.y=element_text(size=14), axis.text.y=element_text(size=13),strip.background = element_rect(fill="white", linewidth=.5, color="black"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0))+labs(y="Mean Percent OM")+geom_text(data=letters001, aes(x=x, y=y, label=label), color='black', inherit.aes=F, fontface='bold', size=7)+scale_fill_manual(values=c("thistle",'steelblue'), breaks=c("initial", 'final'))


#### end ####

#### Figure 4b: Percent C  ####
cn_initial_c1<-mutate(cn_initial_c, time = 'initial')
cn_final_c1<-mutate(cn_final_c, time = 'final')

cn5<-rbind(cn_initial_c1, cn_final_c1)
cn5$time<-factor(cn5$time, levels = c("initial", 'final'))

letters002<-tibble(
  species = c("Gingko","Sweetgum", "Red Oak", "Sycamore"),
  x=c(1.5,1.5,1.5,1.5),
  y=c(50,50,50,50),
  label=c("***",'***','','***')
)
ggplot(data=cn5, aes(y=mean_percentc, x=time, fill=time))+geom_col(show.legend=F, color='black', alpha=0.8)+facet_wrap(~factor(species, levels=c("Gingko","Sweetgum","Red Oak",'Sycamore')), ncol=4)+geom_errorbar(aes(ymin = mean_percentc - sd_percentc, ymax = mean_percentc+sd_percentc), width=0.2, color='black')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=13,color='black'), axis.title.y=element_text(size=14), axis.text.y=element_text(size=13),strip.background = element_rect(fill="white", linewidth=.5, color="black"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0))+labs(y="Mean Percent C")+geom_text(data=letters002, aes(x=x, y=y, label=label), color='black', inherit.aes=F, fontface='bold', size=7)+scale_fill_manual(values=c("thistle",'steelblue'), breaks=c("initial", 'final'))

#### end ####

#### Figure 4c: Percent N  ####
cn_initial_n1<-mutate(cn_initial_n, time = 'initial')
cn_final_n1<-mutate(cn_final_n, time = 'final')

cn8<-rbind(cn_final_n1, cn_initial_n1)
cn8$time<-factor(cn8$time, levels = c("initial", 'final'))


letters003<-tibble(
  species = c("Gingko","Sweetgum", "Red Oak", "Sycamore"),
  x=c(1.5,1.5,1.5,1.5),
  y=c(1.6,1.6,1.6,1.6),
  label=c("***",'','***','')
)
ggplot(data=cn8, aes(y=mean_percentn, x=time, fill=time))+geom_col(show.legend=F, color='black', alpha=0.8)+facet_wrap(~factor(species, levels=c("Gingko","Sweetgum","Red Oak",'Sycamore')), ncol=4)+geom_errorbar(aes(ymin = mean_percentn - sd_percentn, ymax = mean_percentn+sd_percentn), width=0.2, color='black')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=13,color='black'), axis.title.y=element_text(size=14), axis.text.y=element_text(size=13),strip.background = element_rect(fill="white", linewidth=.5, color="black"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0))+labs(y="Mean Percent N")+geom_text(data=letters003, aes(x=x, y=y, label=label), color='black', inherit.aes=F, fontface='bold', size=7)+scale_fill_manual(values=c("thistle",'steelblue'), breaks=c("initial", 'final'))

#### end ####

#### Figure 4d: C:N ####
cn_final5<-mutate(cn_final2, time = 'final')
cn_final6<-mutate(cn_final5, mean_cn = mean_cn_fin, sd_cn = sd_cn_fin)
cn_final7<-dplyr::select(cn_final6, 'species','mean_cn','sd_cn', 'time')

cn_initial6<-mutate(cn_initial3, time = 'initial')
cn_initial7<-mutate(cn_initial6, mean_cn = mean_cn_in, sd_cn = sd_cn_in)
cn_initial8<-dplyr::select(cn_initial7, 'species','mean_cn','sd_cn', 'time')

cn3<-rbind(cn_final7, cn_initial8)
cn3$time<-factor(cn3$time, levels=c("initial",'final'))

letters000<-tibble(
  species=c("Gingko","Sweetgum","Red Oak",'Sycamore'),
  x=c(1.5,1.5,1.5,1.5),
  y=c(75,75,90,75),
  label=c('***','','***','')
)


ggplot(data=cn3, aes(y=mean_cn, x=time, fill=time))+geom_col(show.legend=F, color='black', alpha=0.8)+facet_wrap(~factor(species, levels=c("Gingko","Sweetgum","Red Oak",'Sycamore')), ncol=4)+geom_errorbar(aes(ymin = mean_cn - sd_cn, ymax = mean_cn+sd_cn), width=0.2, color='black')+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=13,color='black'), axis.title.y=element_text(size=14), axis.text.y=element_text(size=13),strip.background = element_rect(fill="white", linewidth=.5, color="black"),  strip.text.x = element_text(size = 13, colour = "black", angle = 0))+labs(y="Mean C:N")+geom_text(data=letters000, aes(x=x, y=y, label=label), color='black', inherit.aes=F, fontface='bold', size=7)+scale_fill_manual(values=c("thistle",'steelblue'), breaks=c("initial", 'final'))

#### end ####

#### Figure 5a: Tea Bag Decay Rate ####
legend2 <- teabags5 %>% distinct(Site, tea_type) %>% 
  mutate(
    fill = case_when(
      tea_type == "Green" ~ 'palegreen3',
      tea_type == "Rooibos" ~ 'orchid4'
    ),
    shape = case_when(
      Site == "CC1" ~ 21,
      Site == "W2" ~ 22,
      Site == "CG3" ~ 24,
      Site == "BA4"~ 23,
    ),
    label = paste0(tea_type, Site) 
  )

legend3<-mutate(legend2, label=case_when(
  label == "GreenBA4"~ "Green: BA",
  label == "GreenCC1"~ "Green: CC",
  label == "GreenCG3"~ "Green: G",
  label == "GreenW2"~ "Green: W",
  label == "RooibosBA4"~ "Rooibos: BA",
  label == "RooibosCC1"~ "Rooibos: CC",
  label == "RooibosCG3"~ "Rooibos: G",
  label == "RooibosW2"~ "Rooibos: W"
  
))



ggplot(teabags5, aes(x=diff_in_days, y=decay_rate_mg_day,fill = interaction(as.factor(Site), as.factor(tea_type)),
                     shape = interaction(as.factor(Site), as.factor(tea_type))))+theme_bw()+theme(legend.box.background =element_rect( color='black'),legend.background=element_blank(), axis.text.x=element_text(size=13),axis.title.x = element_text(size=14), axis.text.y=element_text(size=13), axis.title.y=element_text(size=14), legend.text=element_text(size=13), legend.title=element_text(size=14), legend.position = 'top', legend.direction='horizontal', panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +labs(x="Days", y="Decay Rate (mg" ~day^-1*")")+geom_point(size=3)+
  scale_shape_manual(
    name = "Site &
    Tea Type", 
    labels = legend3 %>% pull(label),
    values = legend3%>% pull(shape)
  )+scale_fill_discrete(
    name = "Site &
    Tea Type", 
    labels = legend3 %>% pull(label),
    type = legend3 %>% pull(fill)
  )+
  geom_line(aes(y = pred, group=tea_type), data = teabags5_greenplot, color = "palegreen3") +
  geom_line(aes(y = pred, group=tea_type), data = teabags5_rooibosplot, color = "orchid4")+
  geom_ribbon(data=teabags5_greenplot,aes(ymin=CI_lwr, ymax=CI_upr, xmin=60, xmax=180, x=diff_in_days, y=decay_rate_mg_day), inherit.aes=F, alpha=0.2, linetype=2)+ annotate("rect", xmin = 170, xmax = 190, ymin = 0, ymax = 8,alpha = 0, color='black')+scale_x_continuous(breaks=c(60,120,180))+
  geom_ribbon(data=teabags5_rooibosplot,aes(ymin=CI_lwr, ymax=CI_upr, xmin=60, xmax=180, x=diff_in_days, y=decay_rate_mg_day), inherit.aes=F, alpha=0.2, linetype=2)+
  geom_text(label="R^2", aes(x=177, y=18.2), size=5, parse=TRUE, color='palegreen3')+
  geom_text(label="= 0.96", aes(x=186, y=18), size=5, color='palegreen3')+
  geom_text(label='k = 0.28,', aes(x=165, y=18), size=5, color='palegreen3')+
  geom_text(label="R^2", aes(x=177, y=16.2), size=5, parse=TRUE, color='orchid4')+
  geom_text(label=" = 0.74", aes(x=186, y=16), size=5, color='orchid4')+
  geom_text(label="k = 0.21,", aes(x=165, y=16), size=5, color='orchid4')
####end####

#### Figure 5b: Comparison - Tea Bag and Natural Leaf Litter Decay ####


tb4<-tb3 %>% group_by(treatment) %>% summarize (mean_decay_rate_mg_day = mean(decay_rate_mg_day),                                                       sd_decay_rate_mg_day= sd(decay_rate_mg_day))


tb4$treatment<-factor(tb4$treatment, levels=c("Green teabags","Labile Carbon", "Rooibos teabags", "Recalcitrant Carbon"))
tb3$treatment<-factor(tb3$treatment, levels=c("Green teabags","Labile Carbon", "Rooibos teabags", "Recalcitrant Carbon"))


letters948a<-tibble(
  treatment = c ("Green teabags", "Labile Carbon Species", "Rooibos teabags", "Recalcitrant Carbon Species"), 
  x = c(1,2,3,4),
  y =  c(8,11, 6.3,7.1), 
  label = c('a','a','b','b')
)          

ggplot(data=tb4, aes(x=treatment, y=mean_decay_rate_mg_day, fill=treatment, color=treatment), show.legend=F)+geom_col(alpha=0.8,show.legend=F)+geom_errorbar(aes(ymin = mean_decay_rate_mg_day - sd_decay_rate_mg_day, ymax = mean_decay_rate_mg_day+ sd_decay_rate_mg_day), width=0.2, color='black')+theme_bw()+geom_point(data=tb3, aes(x=treatment, y=decay_rate_mg_day, fill=treatment), shape=21, size=2, show.legend=F)+scale_color_manual(breaks=c("Green teabags","Labile Carbon","Rooibos teabags","Recalcitrant Carbon"), values=c("black","black",'black','black'))+scale_fill_manual(breaks=c("Green teabags","Labile Carbon","Rooibos teabags","Recalcitrant Carbon"), values=c('palegreen3','aquamarine2','orchid4','pink4'))+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.x=element_text(size=13, color='black'), axis.title.x=element_blank(),axis.title.y=element_text(size=14), axis.text.y=element_text(size=13))+labs(y="Decay Rate (mg" ~day^-1*")")+ geom_text(data=letters948a, aes(x=x, y=y, label=label), size=5, fontface='bold', inherit.aes=FALSE)

#### end ####



