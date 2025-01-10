### Use of exclusive roosting ranges in disc-winged bats, Thyroptera tricolor
#Authors: Silvia Chaves-Ramírez, Maria Sagot, Hellen Solís Hernández, Gloriana Chaverri

#Install and load the following libraries


# Statistical analysis

install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg"))
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(MASS)
library(fitdistrplus)
library(multcomp)
library(ggstatsplot)
library(extrafont)

### DATA 
list.files()

data2 <- read.csv("leaves_inHR.csv"  , 
                  stringsAsFactors = FALSE)

data3 <- read.csv("HR_Gsize.csv"    , 
                  stringsAsFactors = FALSE) 


data4 <- read.csv("TGH_prom.csv"    , 
                  stringsAsFactors = FALSE) 


### Analysis and results####

## Analysis: Daily leaf density available for each group in their roosting ranges

# data2: This object contains information about the size of the roosting ranges for each group of Thyroptera tricolor, 
#as well as the number and density of available shelters (leaves) within each roosting range.

#Normality tests for data 2
head(data2)
class(data2$ID_group) 
data2$ID_group <- factor(data2$ID_group)
class(data2$N_leaves)
data2$N_leaves <-  as.numeric(data2$N_leaves)
class(data2$den_leavesHR)
data2$Gcrow_size <-  as.numeric(data2$Gcrow_size)

shapiro.test(data2$N_leaves) # no hay normalidad
rnorm(data2$N_leaves)
hist(data2$N_leaves)
qqnorm(data2$N_leaves)
shapiro.test(data2$den_leavesHR)# no hay normalidad
hist(data2$den_leavesHR)
shapiro.test(data2$N_leaves)

#To compare distributions
descdist(data2$den_leavesHR, discrete = FALSE)
pois_dist <- fitdist(data2$den_leavesHR, "pois",  method = "mme")
pois_dist
nbin_dist <- fitdist(data2$den_leavesHR, "nbinom", method = "mme")
nbin_dist
 
par(mfrow=c(1,2), mar=c(4,4,2,2))
cdfcomp(list(pois_dist,nbin_dist))
qqcomp(list(pois_dist,nbin_dist))
par(mfrow=c(1,1)) #The leaf density data "den_leavesHR" best fit a negative binomial distribution

# Negative Binomial Generalized Linear Model

#To find out if the density of available shelters (leaves) within the HR of the 11 groups differs
nb1 <- glm.nb(den_leavesHR ~ ID_group, data = data2)
nb1
summary(nb1)
confint(nb1)
anova(nb1,test ="Chi")
plot(nb1,which=c(1,2,4))
plot(residuals(nb1) ~ data2$ID_group); abline(0,0)
par(mfrow=c(1,1))
drop1(nb1, test = "LRT") # para  ver  el efecto Buscarque es #LRT
drop1(nb1, test = "F")

summary(glht(nb1, mcp(ID_group="Tukey")))

#Plot Negative Binomial Generalized Linear Model
#Figure 4. Daily leaf density available for each group in their roosting ranges
ggbetweenstats(data2, ID_group,den_leavesHR,
 plot.type = "boxviolin",
               type = "nonparametric",
               pairwise.comparisons = FALSE,
               pairwise.display = "s",
               p.adjust.method = "holm",
               #effsize.type = "unbiased",
               #bf.prior = 0.707,
               #bf.message = TRUE,
               results.subtitle = TRUE,
               xlab = "Group",
               ylab = "Leaf density (leaves/ha)",  
              centrality.point.args = list(size = 2, color = "#CD0000"),
               centrality.label.args = list(size = 4, nudge_x = 0.4, segment.linetype = 2,
              min.segment.length = 0),  ggtheme = ggplot2::theme_classic()) +        
              ggplot2::scale_color_manual(values = c("#006400","#006400","#006400", "#006400", "#006400","#006400", "#006400","#006400", "#006400", "#006400","#006400"))+
              theme(axis.title.y = element_text( vjust=1.5, size=rel(2.5))) + theme(axis.title.x = element_text( vjust=1.5, size=rel(2.5))) +  theme(axis.text= element_text(size= 16))

#ggsave("Figure 4.png", height = 8, width = 15, dpi = 600)


# data3: This object contains information about average leaf values, group size and overlap between roosting ranges
## Analysis:  Regression between roosting range size and group size

head(data3)
class(data3$HR)
class(data3$G_crowding2)
data3$G_crowding2 <-  as.numeric(data3$G_crowding2)

#Normality tests for data 3
shapiro.test(data3$G_crowding2) 
shapiro.test(data3$HR) #there is normality in all data sets

#Regression between roosting range size and group size
x <- lm(HR ~ G_crowding2, data=data3)
x
summary(x)
#F-statistic:  0.0003425 on 1 and 9 DF,  p-value: p-value: 0.985
#Residual standard error: 0.1114 on 9 degrees of freedom
#Multiple R-squared:  3.805e-05,	Adjusted R-squared:  -0.1111 

#Plot regression between roosting range size and group size

x1 <- ggplot (data3, aes (x=G_crowding2, y=HR)) +
  labs(x = "Group size",y = "Roosting range (ha)")  +  
  geom_point(color = "purple3", size=3) + 
  theme_classic()
X1 <-x1 + theme(axis.title.y = element_text( vjust=1.5, size=rel(2.2))) + 
  theme(axis.title.x = element_text( vjust=1.5, size=rel(2.2))) +
  theme(axis.text= element_text(size= 16))

ggsave("roosting range and group size.png", height = 7, width = 8, dpi = 600)

#Analysis: Average daily density of leaves available for roosting within the roosting ranges according to group size. 
#Average values for both variables

 c <-cor.test(data4$DenH_prom, data4$G_crowding2)
 c
   plot(data4$G_crowding2, data4$DenH_prom)

# Plot Average daily density of leaves available for roosting within the roosting ranges according to group size.
  
 cd <- ggplot (data4, aes (x=G_crowding2, y=DenH_prom )) +
     xlab("Group zise") + 
      ylab("Average daily density of leaves  (leaves/ha)") +
       geom_point (color = "purple3", size=3) +
      theme_classic()
   cd
   
 M  <-   cd + theme(axis.title.y = element_text( vjust=1.5, size=rel(2.2))) + 
   theme(axis.title.x = element_text( vjust=1.5, size=rel(2.2))) +
   theme(axis.text= element_text(size= 16))
   ggsave("Den_tam.png", height = 7, width = 8, dpi = 600)   
    

##Analysis: Average daily number of leaves available according to roosting range size (ha)

sprom <-lm( data4$Nleaves ~  data4$Area) 
sprom
summary(sprom)
o <- ggplot (data4, aes (x=Area, y=Nleaves)) + xlab("Roosting range(ha)") + ylab("Average daily number of leaves")+
    geom_point(color = "purple3", size=3) + 
     theme_classic()

O <-o + theme(axis.title.y = element_text( vjust=1.5, size=rel(2.2))) + 
  theme(axis.title.x = element_text( vjust=1.5, size=rel(2.2))) + 
  theme(axis.text= element_text(size= 16))
O
 ggsave("NhojasyTerriA.png", height = 7, width = 8, dpi = 600)  

  
#Figure 3. A) Size of the roosting range (ha) according to group size. 
#B) Average daily density of leaves available for roosting within the roosting ranges according to group size. 
#C) Average daily number of leaves available according to roosting range size (ha).
  
ggarrange(X1 ,M ,O , ncol = 3,
            labels = c("A)","B)", "C)"),
            font.label = list(size = 14, color = "black")) 

ggsave("Figure 3.", height = 8, width = 19, dpi = 600)   
  

## Analysis: Percent overlap between roosting ranges according to average density of tubular leaves suitable for roosting in Thyroptera tricolor

# data4: This object contains information about percent overlap between roosting ranges and average density of tubular leaves suitable for roosting

data4
head(data4)
class(data4$Overlapping)
class(data4$G_crowding2)
data4$Crowding <-  as.numeric(data4$G_crowding2)

# Linear Model
r  <-lm(Overlapping  ~ DenH_prom, data=data4) 
summary(r)
#Create a new variable to make a quadratic model
data4$DenH_prom2 <- data4$DenH_prom^2
# Fit quadratic regression model
quad <- lm(Overlapping  ~ DenH_prom + DenH_prom2, data=data4)
quad

summary(quad)
# R-squared:  0.5843,	Adjusted R-squared:  0.4803 
#F-statistic: 5.621 on 2 and 8 DF,  p-value: 0.02987


## Plot quadratic regression model
# Figure 5 .Percent overlap between roosting ranges according to average density of tubular leaves suitable for roosting in Thyroptera tricolor

z <- ggplot (data4, aes (x=DenH_prom)) +
  xlab("Average leaf density (leaves/ha)") +
  ylab("overlap between roosting ranges  (%)")+
  geom_point(aes(y= Overlapping), color = "purple3", size= 2.5) +
  stat_smooth(aes(y = Overlapping),method = "lm",col="black", formula = y ~ x + I(x^2), size = 0.60, alpha= 0.2) +
  theme_classic()
z + theme(axis.title.y = element_text( vjust=1.5, size=rel(2.5))) + theme(axis.title.x = element_text( vjust=1.5, size=rel(2.5))) +
    theme(axis.text= element_text(size= 15))

ggsave("Figure5.png", height = 8, width = 11, dpi = 600)

 
