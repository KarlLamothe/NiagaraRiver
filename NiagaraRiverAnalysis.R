################################################################################
################################################################################
# Spatial and temporal diversity trends in the Niagara River fish assemblage
# Karl A. Lamothe, Justin A. G. Hubbard, D. Andrew R. Drake
# R Code prepared by Karl A. Lamothe, PhD - Karl.Lamothe@dfo-mpo.gc.ca
# 2020-10-28 revision; R version 4.0.2
################################################################################
################################################################################
# Load libraries
library(pacman) # For p_load function
p_load(xlsx)    # For importing xlsx documents
p_load(ggplot2) # For plotting
p_load(ggrepel) # For plotting
p_load(psych)   # For correlation calculations
p_load(vegan)   # For multivariate analyses
p_load(FD)      # For functional diversity analysis
p_load(corrplot)# For correlation analysis
p_load(tidyr)   # For reshaping Data
p_load(dplyr)   # For reshaping Data

# Personal ggplot theme
theme_me <- theme_bw() + 
  theme(axis.title=element_text(family="sans", colour="black"),
        axis.text.x=element_text(size=10, family="sans", colour="black"),
        axis.text.y=element_text(size=10, family="sans", colour="black"),
        legend.title=element_text(size=10, family="sans", colour="black"),
        legend.text=element_text(size=10, family="sans", colour="black"),
        panel.border = element_rect(colour = "black", fill=NA),
        axis.ticks = element_line(colour="black"))

# To plot multiple ggplots in one pane
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))}
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))}}
}

################################################################################
# Load data
################################################################################
Fish.traits <- read.xlsx("NiagaraRiverFishDiv.xlsx", header=T, 
                         sheetName= "NiagaraStudyTraits") 
Habitat <- read.xlsx("NiagaraRiverFishDiv.xlsx", header=T, 
                     sheetName= "Habitat") 
Fish.comm <- read.xlsx("NiagaraRiverFishDiv.xlsx", header=T, 
                       sheetName= "FishCommunity") 

################################################################################
# Check data
################################################################################
head(Fish.traits) # Data are ordered by common name
str(Fish.traits)

# Fish community CPUE data frame
Fish.comm.CPUE <- Fish.comm[10:74]/Fish.comm$Effort_sec # CPUE
colnames(Fish.comm.CPUE) <- Fish.traits$CODE # Shortens species names

# Observations summed across 1000 m transects
Fishagg<-aggregate(Fish.comm[8:74], list(Site.No = Fish.comm$Site.No,
                                         Year   = Fish.comm$Year,
                                         Season = Fish.comm$Season), sum)
head(Fishagg)
Fishagg <- Fishagg[-5] # Remove aggregated sampling pass variable

# Order observations by Site
Fishagg <- Fishagg[order(Fishagg$Site.No),]

# Set up data frame for future RDA
Fish.RDA <- Fishagg[5:69]
colnames(Fish.RDA) <- Fish.traits$CODE

# Presence absence data frame
Fish.comm.PA <- Fish.comm.CPUE
Fish.comm.PA[Fish.comm.PA > 0] <- 1

################################################################################
################################################################################
############ Figure 1. Abundance and CPUE of species ###########################
################################################################################
################################################################################
#Total Counts
Fish.Counts <- as.data.frame(cbind(Count=colSums(Fishagg[5:69]), 
                                   Common=as.character(Fish.traits$COMMONNAME)))

# Ensure proper class
str(Fish.Counts)
Fish.Counts$Count <- as.numeric(Fish.Counts$Count)

# Only plot species where more than 100 were caught
Fish.Counts2 <- Fish.Counts[Fish.Counts$Count>100,]

################################################################################
# Plot Total Catch
################################################################################
Fish.Counts2$Common <- factor(Fish.Counts2$Common, levels = Fish.Counts2$Common[order(Fish.Counts2$Count)])

Countplot<-ggplot(Fish.Counts2, aes(x = reorder(Common,Count), y = Count))+
  geom_bar(stat="identity")+ylab("Total caputured")+ xlab("Species")+
  coord_flip()+
  theme_me
Countplot

################################################################################
# CPUE by river section 
################################################################################
# Lower section
Fish.comm.CPUE.lower <- Fish.comm[Fish.comm$Section=="Lower",]
Fish.comm.CPUE.lower <- Fish.comm.CPUE.lower[10:74]/Fish.comm.CPUE.lower$Effort

# Upper section
Fish.comm.CPUE.upper <- Fish.comm[Fish.comm$Section=="Upper",]
Fish.comm.CPUE.upper <- Fish.comm.CPUE.upper[10:74]/Fish.comm.CPUE.upper$Effort

# Combine
Fishes<-as.data.frame(rbind(Fish.comm.CPUE.lower,Fish.comm.CPUE.upper))
Section<-c(rep("Lower",length(Fish.comm.CPUE.lower$Alewife)),
           rep("Upper",length(Fish.comm.CPUE.upper$Alewife)))
Fishes<-cbind(Fishes, Section)

# Make Data frames
Fish.CPUE.mean.section<-as.data.frame(aggregate(Fishes[1:65],
                                                list(Fishes$Section),mean))
Fish.CPUE.mean.section.t<-t(Fish.CPUE.mean.section) #transpose
Fish.CPUE.mean.section.t<-Fish.CPUE.mean.section.t[-1,] # remove first row
Species<-colnames(Fish.CPUE.mean.section[2:66]) #Species
Fish.CPUE.mean.section.t<-cbind(Fish.CPUE.mean.section.t,Species) # combine 
colnames(Fish.CPUE.mean.section.t)<-c("Lower","Upper","Species")
Fish.CPUE.mean.section.t<-as.data.frame(Fish.CPUE.mean.section.t)

# Check class
str(Fish.CPUE.mean.section.t)
Fish.CPUE.mean.section.t$Upper<-as.numeric(Fish.CPUE.mean.section.t$Upper)
Fish.CPUE.mean.section.t$Lower<-as.numeric(Fish.CPUE.mean.section.t$Lower)

# Remove periods and substitute space for species names
Fish.CPUE.mean.section.t$Species<-gsub(".", " ",
                                       Fish.CPUE.mean.section.t$Species, 
                                       fixed=TRUE)
head(Fish.CPUE.mean.section.t)

# Only use the same fishes from Total catch plot
CPUE2<-Fish.CPUE.mean.section.t[match(Fish.Counts2$Common, 
                                      Fish.CPUE.mean.section.t$Species), ]
CPUE3<-as.data.frame(cbind(Species=rep(CPUE2$Species,2),
                           CPUE = c(CPUE2$Upper,CPUE2$Lower),
                           Section = c(rep("Upper",length(CPUE2$Upper)),
                                       rep("Lower",length(CPUE2$Lower))
                           )))
CPUE3$CPUE<-as.numeric(CPUE3$CPUE)

################################################################################
# Plot CPUE by section
################################################################################

CPUEplot<-ggplot(CPUE3, aes(reorder(Species, CPUE), CPUE, fill=Section))+ 
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=c("black","grey"))+
  ylab("CPUE") + xlab("Species") +
  coord_flip()+ 
  theme_me+ 
  theme(
    legend.position = c(.95, .25),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.background = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  )
CPUEplot
#dev.off()

#tiff("Figure1.tiff", width = 6, height = 4.5, units = 'in', res = 1000)
multiplot(Countplot,CPUEplot,cols=2)
#dev.off()

################################################################################
################################################################################
#################### Summarize taxonomic diversity #############################
################################################################################
################################################################################
length(unique(Fish.traits$GENUS))
length(unique(Fish.traits$FAMILY))

# Summarize abundance of species
colnames(Fish.comm)
colSums(Fish.comm[c(10:74)]) #species specific catch

# Introduced, reintroduced, reportedly introduced, probably introduced
Introduced<-colnames(Fish.comm[c(10,12,13,22,25,26,27,37,55,56,58,59,60,71)])
Introduced

# Total number caught
sum(Fish.comm[c(10,12,13,22,25,26,27,37,55,56,58,59,60,71)]) #3568

# Proportion of catch
sum(Fish.comm[c(10,12,13,22,25,26,27,37,55,56,58,59,60,71)])/sum(Fish.comm[c(10:74)])
# 0.086

# Number of observations for each species that was caught
colnames(Fishagg)
Fish.agg.PA <- Fishagg[5:69]
Fish.agg.PA[Fish.agg.PA > 0] <- 1
Fish.agg.PA <- cbind(Fish.agg.PA, Site.No = Fishagg$Site.No)

# Number of sites present
colSums(aggregate(Fish.agg.PA[1:65], list(Site.No = Fish.agg.PA$Site.No), sum))

# Total number of white sucker captured
sum(Fish.comm$White.Sucker)

# Proportion of sites white sucker captured
84/88 #95.5%

# Emerald Shiner
sum(Fish.comm$Emerald.Shiner)
77/88 #87.5

# Yellow Perch
sum(Fish.comm$Yellow.Perch)
79/88 #89.8

# Seasonal species catch
Seasonal<-aggregate(Fish.comm[c(10:74)], by = list(Year = Fish.comm$Year,
                                                   Season = Fish.comm$Season,
                                                   Section = Fish.comm$Section), 
                    FUN = sum)
Seasonal

z<-aggregate(Fish.comm[c(10:74)], by = list(Season=Fish.comm$Section), FUN = sum)
rowSums(z[2:66])

z<-(aggregate(Fish.comm[c(10:74)], by = list(Season=Fish.comm$Year), FUN = sum))
rowSums(z[2:66])

z<-aggregate(Fish.comm[c(10:74)], by = list(Season=Fish.comm$Season), FUN = sum)
rowSums(z[2:66])

rowSums(Seasonal[4:68])

# Number of species per section
Section2<-aggregate(Fish.comm[c(10:74)], by = list(Season=Fish.comm$Section), 
                    FUN = sum)
Section2[Section2>0] <-1
rowSums(Section2[2:66])

#########################################################################
#########################################################################
##################### Prepare data for RDA ##############################
#########################################################################
#########################################################################
colnames(Habitat)
Habitat1<-as.data.frame(cbind(Site.No = Habitat$Site.No,
                              Season  = as.character(Habitat$Season),
                              Year    = as.character(Habitat$Year),
                              Section = as.character(Habitat$Section),
                              Temp    = Habitat$WaterTemp,
                              Cond    = Habitat$Conductivity,
                              DO      = Habitat$DO, 
                              Turb    = Habitat$Turbidity,
                              Depth   = Habitat$AvDepth,
                              WV      = Habitat$AvWaterVelocity,
                              Veg     = Habitat$Submerged,
                              dMouth  = Habitat$dMouth))
str(Habitat1)

# Convert back to numeric variables
Columns<-colnames(Habitat1[5:12])
Habitat1[Columns] <- sapply(Habitat1[Columns],as.numeric)

# Calculate pearson correlations
cor_mat <- cor(Habitat1[,c(5:12)], method='pearson')

###############################################################################
# Figure S2 ###################################################################
###############################################################################
par(xpd=T)
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=1, number.cex=1, addCoefasPercent=T,
               mar = c(1, 1, 4, 1), tl.col="black")

################################################################################
# Table 1 ######################################################################
################################################################################
length(Habitat1$Section[Habitat1$Section=="Lower"])
length(Habitat1$Section[Habitat1$Section=="Upper"])
aggregate(.~Section, Habitat1[4:12], mean)
aggregate(.~Section, Habitat1[4:12], min)
aggregate(.~Section, Habitat1[4:12], max)

# Aggregate per site
head(Habitat1)
Hab.agg <- aggregate(Habitat1[c(1,5:12)], list(Site.No = Habitat1$Site.No,
                                               Year   = Habitat1$Year,
                                               Season = Habitat1$Season), mean)
head(Hab.agg)
Hab.agg <- Hab.agg[-c(4)] # Remove aggregated site variable 

# Order the data by Site No 
Hab.agg$Site.No<-as.numeric(Hab.agg$Site.No)
Hab.agg <- Hab.agg[order(Hab.agg$Site.No),]

# Scaling continuous covariates to improve model fit and interpretation
Habitat.RDA <- as.data.frame(scale(Hab.agg[4:11], center = TRUE))
head(Habitat.RDA)

################################################################################
################################################################################
#################### Redundancy Analysis (RDA) #################################
################################################################################
################################################################################
#Transform fish data to reduce effects of abundant species
FishesTransformed <- decostand(Fish.RDA, method = "hellinger")

#reproducible results (unnecessary for performing analysis)
set.seed(2336) 

# Stepwise model selection
head(Habitat.RDA)
mod0<-rda(FishesTransformed ~ 1, Habitat.RDA)
mod1<-rda(FishesTransformed ~ ., Habitat.RDA)
rda_select.r <-ordistep(mod0, scope = formula(mod1), direction = "both", 
                        Pin = 0.05, Pout = 0.10, perm.max = 9999)

# Run final model
NR.rda <-rda(FishesTransformed ~ Depth + Temp + dMouth + Veg, 
             data = Habitat.RDA)
summary(NR.rda)

################################################################################
############################### Figure 2: Triplot ##############################
################################################################################
# Observation scores (site scores)
scores <- data.frame(Habitat.RDA,NR.rda$CCA$u)

# Species scores
vscores <- data.frame(NR.rda$CCA$v)

# Covariate scores
var.score <- data.frame(NR.rda$CCA$biplot[,1:2]) 
var.score.sc <- var.score * 0.8 # Scaling covariates for plotting purposes
var.score.sc$variables <- c("Depth", "Temp", "dMouth", "Veg")

# Only plotting the most abundant species for visualization purposes
rownames(vscores)
vscores2.sc <- vscores[c(6,9,18,21,25,30,35,46,54,56,63,65),]
Section <- c(rep("Upper", 52), rep("Lower", 36))

# Create plot of model
ggRDA <- ggplot(scores, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(aes(shape=Section, color=Hab.agg$Season)) +
  scale_shape_manual(name = "River Section", values = c(15:17))+
  scale_color_manual(name = "Season", 
                     values=c("black","lightgrey","darkgrey"),
                     labels = c("Fall","Spring","Summer"))+
  geom_segment(data = vscores2.sc, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow=arrow(length=unit(0.2,"cm")), 
               color = "black",inherit.aes = FALSE,lwd=0.25) +
  geom_text(data = vscores2.sc, 
            aes(x = RDA1, y = RDA2, label = rownames(vscores2.sc)), 
            col = 'black', inherit.aes = FALSE, 
            nudge_y = ifelse(vscores2.sc$RDA2 > 0, 0.02, -0.02),
            nudge_x = ifelse(vscores2.sc$RDA1 > 0, 0.02, -0.02),size=3)+
  geom_segment(data = var.score.sc, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length=unit(0.2,"cm")), 
               color = 'black', inherit.aes = FALSE, lwd=0.25) +
  geom_text(data = var.score.sc, 
            aes(x = RDA1, y = RDA2, label = variables), 
            col = 'black', inherit.aes = FALSE, 
            nudge_y = ifelse(var.score.sc$RDA2 > 0, 0.02, -0.02),
            nudge_x = ifelse(var.score.sc$RDA1 > 0, 0.02, -0.02),size=3) +
  labs(x = "RDA Axis 1", y = "RDA Axis 2") + 
  theme_me + theme(legend.justification = c(1,0), 
                   legend.position = c(1,0),
                   legend.background = element_blank(),
                   legend.box = 'horizontal')
ggRDA

#tiff("ggRDA.tiff", width = 6, height = 4, units = 'in', res = 1000)
#ggRDA
#dev.off()

################################################################################
################################################################################
############ Non-metric multidimensional scaling (NMDS) ########################
################################################################################
################################################################################
# Prepare data
colnames(Fish.RDA)

#habitat data
colnames(Habitat.RDA)
Habitat.NMDS<-Habitat.RDA[c(1,5,7,8)]

#comm data
NMDSdist<-decostand(Fish.RDA, method="hellinger")

#Final
NMDSord <- metaMDS(NMDSdist,k=2, try=20, trymax=1000) #final used
NMDSord

# Stressplot
stressplot(NMDSord)

# Shepards test/goodness of fit
goodness(NMDSord) # Produces a results of test statistics for goodness of fit for each point

# fit environmental variables
en = envfit(NMDSord, Habitat.NMDS, permutations = 9999, na.rm = TRUE,
            choices=c(1,2))
en

# data
data.scores = as.data.frame(scores(NMDSord))
head(data.scores)
data.scores$season = Hab.agg$Season
en_coord_cont = as.data.frame(scores(en, "vectors")) 

# for plotting species
Species<-as.data.frame(NMDSord$species)
Spec<-rownames(Species)
Count<-Fish.Counts$Count
Species<-cbind(Species, Spec,Count)
Species2<-Species[Species$Count>25,]
Species3<-Species[Species$Count>600,] #only plotting abundant species

################################################################################
############################ Figure S3 #########################################
################################################################################
ggNMDSENV <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(aes(color=Hab.agg$Season, shape=Section),
             size=3, alpha=0.5) +
  scale_shape_manual(name = "River Section", values = c(15:17))+
  scale_colour_manual(values = c("orange", "violet","green"),
                      labels = c("Fall","Spring","Summer")) +
  geom_segment(data=Species3, aes(x=0,y=0,xend=MDS1, yend=MDS2), 
               arrow=arrow(length=unit(0.2,"cm")),
               size =0.5, colour = "black") + 
  geom_text(data=Species3, aes(x=MDS1, y=MDS2),
            label=Species3$Spec, check_overlap = F, 
            nudge_y = ifelse(Species3$MDS2 > 0, 0.03, -0.03),
            nudge_x = ifelse(Species3$MDS1 > 0, 0.03, -0.03))+
  geom_segment(aes(x = 0, y = 0, 
                   xend = NMDS1, yend = NMDS2), 
               arrow=arrow(length=unit(0.2,"cm")),
               data = en_coord_cont, size =0.5, colour = "black") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
            colour = "black", 
            label = row.names(en_coord_cont),
            nudge_y = ifelse(en_coord_cont$NMDS2 > 0, 0.03, -0.03),
            nudge_x = ifelse(en_coord_cont$NMDS1 > 0, 0.03, -0.03)) + 
  labs(colour = "Season")+
  theme_me+
  theme(legend.position = c(.95, 0.75),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))

ggNMDSENV

################################################################################
################################################################################
################ Functional Diversity Analysis #################################
################################################################################
################################################################################
# Provide correlation of trait data
head(Fish.traits)
Traitscorr<-Fish.traits[-c(1:8)]
colSums(Traitscorr[2:20])/65
str(Traitscorr)

################################################################################
########################## Figure S4 ###########################################
################################################################################
cor_mat <- cor(Traitscorr, method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=1, number.cex=1, addCoefasPercent=T,
               mar = c(1, 1, 4, 1), tl.col="black")

# Reduce three trait categories into respective trait dimensions
# Diet analysis - what the species eat
# Removed correlated variables and variables without variation
Diet.data <- as.data.frame(cbind(Algae = Fish.traits$ALGPHYTO,
                                 Macrophyte = Fish.traits$MACVASCU, 
                                 Fish       = Fish.traits$FSHCRCRB,
                                 Eggs       = Fish.traits$EGGS))

# PCA of diet preference
F.transformed1 <- decostand(Diet.data, method = "hellinger")
PCA1 <- princomp(F.transformed1, cor=FALSE, scores=TRUE)
summary(PCA1)

# Quick plots
par(xpd=F)
plot(PCA1)
biplot(PCA1, xlim = c(-.4,.4), ylim = c(-.4,.3), cex = 0.8)
abline(v=0,h=0,lty="dashed")

# Substrate analysis
Habitat.data <- as.data.frame(cbind(CLAYSILT = Fish.traits$CLAYSILT, 
                                    SAND     = Fish.traits$SAND,
                                    GRAVEL   = Fish.traits$GRAVEL,
                                    COBBLE   = Fish.traits$COBBLE,
                                    BOULDER  = Fish.traits$BOULDER,
                                    BEDROCK  = Fish.traits$BEDROCK,
                                    VEGETAT  = Fish.traits$VEGETAT,
                                    LWD      = Fish.traits$LWD,
                                    DEBRDETR = Fish.traits$DEBRDETR))

# PCA on Substrate data
F.transformed2 <- decostand(Habitat.data, method = "hellinger")
PCA2 <- princomp(F.transformed2, cor = FALSE, scores = TRUE)
summary(PCA2)

# Quick plots
plot(PCA2)
biplot(PCA2, xlim = c(-.4,.4), ylim = c(-.4,.3), cex = 0.8)
abline(v=0,h=0,lty="dashed")

# Reproduction analysis
Reprod<-Fish.traits[c(17,18)]
F.transformed3 <- decostand(Reprod, method="hellinger")
PCA3<-princomp(F.transformed3, cor=FALSE, scores = TRUE)
summary(PCA3)

# Quick plots
plot(PCA3)
biplot(PCA3, xlim = c(-.4,.4), ylim = c(-.4,.3), cex = 0.8)
abline(v=0,h=0,lty="dashed")

# Fuction to assess significance of the principal components.
sign.pc <-function(x, R=9999, s=10, cor=T,...){
  pc.out <- princomp(x, cor=cor,...)  # run PCA
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]  
  pve.perm <- matrix(NA,ncol=s,nrow=R)
  for(i in 1:R){
    x.perm <- apply(x,2,sample)# permutate each column
    pc.perm.out <- princomp(x.perm,cor=cor,...)# run PCA
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s] 
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R # calcalute the p-values
  return(list(pve=pve,pval=pval))
}

# Apply the significance function
sign.pc(F.transformed1, cor = FALSE) #Axis 1 is significant
sign.pc(F.transformed2, cor = FALSE) #Axes 1 + 2 are significant
sign.pc(F.transformed3, cor = FALSE) #Axis 1 is significant

################################################################################
# Create full data frame for final ordination analysis
Fish.traits.reduced<-data.frame(Spp          = Fish.traits$CODE,
                                Habitat1     = PCA2$scores[,1],
                                Habitat2     = PCA2$scores[,2],
                                Diet         = PCA1$scores[,1],
                                Reproduction = PCA3$scores[,1],
                                Size = scale(Fish.traits$AV.TL.CM, 
                                             scale = TRUE, center = TRUE))

rownames(Fish.traits.reduced)<-Fish.traits$CODE

# Look at the variability across data axes
sd(Fish.traits.reduced$Habitat1)
sd(Fish.traits.reduced$Habitat2)
sd(Fish.traits.reduced$Diet)
sd(Fish.traits.reduced$Reproduction)
sd(Fish.traits.reduced$Size)

summary(Fish.traits.reduced$Habitat1)
summary(Fish.traits.reduced$Habitat2)
summary(Fish.traits.reduced$Diet)
summary(Fish.traits.reduced$Reproduction)
summary(Fish.traits.reduced$Size)

################################################################################
# Calculate Functional Diversity Measures
################################################################################
Fishfunction<-dbFD(x = Fish.traits.reduced[2:6], a = Fish.RDA,
                   w.abun = TRUE, stand.x = TRUE, calc.FRic = TRUE, m = "max", 
                   stand.FRic = FALSE, scale.RaoQ = FALSE, 
                   print.pco = TRUE, calc.FGR = FALSE, messages = TRUE)

# Extract diversity measures per fish community
SpeciesRichness      <- Fishfunction$nbsp     # Species Richness
FunctionalDivergence <- Fishfunction$FDiv     # Functional divergences
FunctionalDispersion <- Fishfunction$FDis     # Functional dispersion
Eigen                <- Fishfunction$x.values # Eigenvalues
Axes                 <- Fishfunction$x.axes

Eigen[1]/sum(Eigen) # 29.3%
Eigen[2]/sum(Eigen) # 23.7%
Eigen[3]/sum(Eigen) # 20.0%

################################################################################
################# Figure 3 - Functional Trait PCOA #############################
################################################################################
Ord1<-ggplot(Axes, aes(x = A1, y = A2, label=rownames(Axes))) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(pch=20) +
  geom_text_repel(min.segment.length = Inf, seed = 42, box.padding = 0.2, size=2) +
  xlab("Component 1 (29.3%)")+
  ylab("Component 2 (23.7%)")+
  theme_me
Ord1


Ord2<-ggplot(Axes, aes(x = A2, y = A3,label=rownames(Axes))) +
  geom_hline(yintercept=0,linetype="dashed",col="black")+
  geom_vline(xintercept=0,linetype="dashed",col="black")+
  geom_point(pch=20) +
  geom_text_repel(min.segment.length = Inf, seed = 42, box.padding = 0.2, size=2) +
  xlab("Component 2 (23.7%)")+
  ylab("Component 3 (20.0%)")+
  theme_me
Ord2

#tiff("Figure3.tiff", width = 6.5, height = 4, units = 'in', res = 1000)
multiplot(Ord1,Ord2,cols=2)
#dev.off()

################################################################################
# Create data frame for post-analyses
Gradient.Frame<-data.frame(Site    = Hab.agg$Site.No,
                           Year    = Hab.agg$Year,
                           Season  = Hab.agg$Season,
                           FDis    = FunctionalDispersion,
                           FDiv    = FunctionalDivergence,
                           SpRich  = SpeciesRichness,
                           Section = Section)
head(Gradient.Frame)
str(Gradient.Frame)
{Gradient.Frame$Year <- as.character(Gradient.Frame$Year)
  Gradient.Frame$Year <- as.factor(Gradient.Frame$Year)}

# Aggregate functional diversity measures
aggregate(Gradient.Frame[4], list(Gradient.Frame$Season, Gradient.Frame$Section, Gradient.Frame$Year), mean)
aggregate(Gradient.Frame[5], list(Gradient.Frame$Season, Gradient.Frame$Section, Gradient.Frame$Year), mean)
aggregate(Gradient.Frame[4], list(Gradient.Frame$Season, Gradient.Frame$Section, Gradient.Frame$Year), sd)
aggregate(Gradient.Frame[5], list(Gradient.Frame$Season, Gradient.Frame$Section, Gradient.Frame$Year), sd)

###############################################################################
# Permanova to look at differences in functional metrics
###############################################################################
# Functional dispersion
perm.fdis <- adonis(FunctionalDispersion ~ Season + Year + Section, 
                    data = Gradient.Frame, permutations = 9999, 
                    method = "euclidean")
perm.fdis

# Functional divergence
perm.fdiv <- adonis(FunctionalDivergence ~ Season + Year + Section, 
                    data = Gradient.Frame, permutations = 9999, 
                    method = "euclidean")
perm.fdiv

################################################################################
# Summary statistics
################################################################################

aggregate(Gradient.Frame$FDis, by = list(Gradient.Frame$Section), FUN = mean)
aggregate(Gradient.Frame$FDis, by = list(Gradient.Frame$Section), FUN = sd)

aggregate(Gradient.Frame$FDis, by = list(Gradient.Frame$Year), FUN = mean)
aggregate(Gradient.Frame$FDis, by = list(Gradient.Frame$Year), FUN = sd)

aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Section), FUN = mean)
aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Section), FUN = sd)

aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Season), FUN = mean)
aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Season), FUN = sd)

aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Year), FUN = mean)
aggregate(Gradient.Frame$FDiv, by = list(Gradient.Frame$Year), FUN = sd)

################################################################################
# Figure 4: Plotting Functional diversity difference
################################################################################
Fdis.year.box <- ggplot(data = Gradient.Frame, 
                        aes(x = Year, y = FunctionalDispersion)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional dispersion") + theme_me +
  theme(panel.background = element_rect(fill = "lightgrey"),
        axis.text.x = element_text(angle = 35, hjust = 1))

Gradient.Frame$Season <- factor(Gradient.Frame$Season , levels=c("SPRING","SUMMER","FALL"))

Fdis.season.box <- ggplot(data = Gradient.Frame, 
                          aes(x = Season, y = FunctionalDispersion)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional dispersion") + theme_me +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

Fdis.section.box <- ggplot(data = Gradient.Frame, 
                           aes(x = Section, y = FunctionalDispersion)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional dispersion") + theme_me+
  theme(panel.background = element_rect(fill = "lightgrey"),
        axis.text.x = element_text(angle = 35, hjust = 1))

Fdiv.year.box <- ggplot(data = Gradient.Frame, 
                        aes(x = Year, y = FunctionalDivergence)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional divergence") + theme_me+
  theme(panel.background = element_rect(fill = "lightgrey"),
        axis.text.x = element_text(angle = 35, hjust = 1))

Fdiv.season.box <- ggplot(data = Gradient.Frame, 
                          aes(x = Season, y = FunctionalDivergence)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional divergence") + theme_me+
  theme(panel.background = element_rect(fill = "lightgrey"),
        axis.text.x = element_text(angle = 35, hjust = 1))

Fdiv.section.box <- ggplot(data = Gradient.Frame, 
                           aes(x = Section, y = FunctionalDivergence)) +
  geom_boxplot() + 
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  xlab("") + ylab("Functional divergence") + theme_me+
  theme(panel.background = element_rect(fill = "lightgrey"),
        axis.text.x = element_text(angle = 35, hjust = 1))

###############################################################################
################ Figure 4 Functional Diversity Metrics ########################
###############################################################################
#tiff("Figure4.tiff", width = 6.5, height = 5, units = 'in', res = 1000)
multiplot(Fdis.year.box,Fdiv.year.box,
          Fdis.season.box,Fdiv.season.box,
          Fdis.section.box,Fdiv.section.box,cols=3)
#dev.off()

###############################################################################
################# Figure 5 Functional dispersion x Functional divergence ######
###############################################################################
fdisfdivplot<-ggplot(Gradient.Frame, aes(x=FDiv, y=FDis, shape=Section, 
                                         color=Season, linetype=Section))+
  geom_point() +  
  scale_shape_manual(name = "River Section", values = c(15:17))+
  scale_color_manual(name = "Season", 
                     values=c("black","darkgrey","lightgrey"),
                     labels = c("Fall","Spring","Summer"))+
  geom_smooth(method='lm', se=F)  +
  labs(y="Functional Dispersion",x="Functional Divergence")+
  scale_linetype_discrete(name="Season", 
                          breaks=c("Fall","Spring","Summer"), 
                          labels = c("Fall", "Spring", "Summer"))+
  theme_me + theme(legend.justification = c("left", "top"), 
                   legend.position = c(0,1),
                   legend.background = element_blank(),
                   legend.box = 'vertical',
                   legend.key = element_rect(fill = NA, colour = NA, size = 0.25))
fdisfdivplot

#tiff("Figure5.tiff", width = 5, height = 3.5, units = 'in', res = 1000)
#fdisfdivplot
#dev.off()

################################################################################
# Plot relationship between functional diversity and habitat variables
plot(Gradient.Frame$FDis~Hab.agg$Temp, pch=20, las=1)
plot(Gradient.Frame$FDis~Hab.agg$Cond, pch=20, las=1)
plot(Gradient.Frame$FDis~Hab.agg$DO, pch=20, las=1)
plot(Gradient.Frame$FDis~Hab.agg$Turb, pch=20, las=1)
plot(Gradient.Frame$FDis~Hab.agg$Depth, pch=20, las=1)
plot(Gradient.Frame$FDis~Hab.agg$WV, pch=20, las=1)

plot(Gradient.Frame$FDiv~Hab.agg$Temp, pch=20, las=1)
plot(Gradient.Frame$FDiv~Hab.agg$Cond, pch=20, las=1)
plot(Gradient.Frame$FDiv~Hab.agg$DO, pch=20, las=1)
plot(Gradient.Frame$FDiv~Hab.agg$Turb, pch=20, las=1)
plot(Gradient.Frame$FDiv~Hab.agg$Depth, pch=20, las=1)
plot(Gradient.Frame$FDiv~Hab.agg$WV, pch=20, las=1)

################################################################################
# Upper versus lower unique species
################################################################################
Fishagg<- cbind(Fishagg,
                Section=Gradient.Frame$Section,
                FDis=Gradient.Frame$FDis,
                FDiv=Gradient.Frame$FDiv,
                SpRich=Gradient.Frame$SpRich)
str(Fishagg)
Fishagg$Section<-as.factor(Fishagg$Section)

#Upper
Fish.Upper<-Fishagg[1:52,]
summary(Fish.Upper$Section)
colnames(Fish.Upper)
Fish.Upper<-Fish.Upper[c(1:4,6,18,19,24,26,27,33,49,62,65,71,72,73)]
head(Fish.Upper)
Abund.Upper<-rowSums(Fish.Upper[5:14])

#lower
Fish.Lower<-Fishagg[53:88,]
summary(Fish.Lower$Section)
colnames(Fish.Lower)
Fish.Lower<-Fish.Lower[c(1:4,7,8,20,21,38,41,55,57,71,72,73)]
head(Fish.Lower)
Abund.Lower<-rowSums(Fish.Lower[5:12])

################################################################################
###########Figure 6 Functional dispersion and unique sp abundance ##############
################################################################################
Dis<-ggplot()+
  geom_smooth(aes(x=Fish.Upper$FDis, y=Abund.Upper),method="lm", col="black", lwd=0.5)+
  geom_point(aes(x=Fish.Upper$FDis, y=Abund.Upper))+
  geom_smooth(aes(x=Fish.Lower$FDis, y=Abund.Lower), method="lm", col="grey", lwd=0.5)+
  geom_point(aes(x=Fish.Lower$FDis, y=Abund.Lower), col="grey")+ylim(-10,50)+theme_me+
  scale_x_continuous(expand=c(0,0), limits=c(0,2.5)) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50)) +
  coord_cartesian(xlim=c(0.8,2.5), ylim=c(0,50))+
  ylab("Abundance")+xlab("Functional Dispersion")
Dis

Div<-ggplot()+
  geom_smooth(aes(x=Fish.Upper$FDiv, y=Abund.Upper), method="lm", col="black", lwd=0.5)+
  geom_point(aes(x=Fish.Upper$FDiv, y=Abund.Upper))+
  geom_smooth(aes(x=Fish.Lower$FDiv, y=Abund.Lower), method="lm", col="grey", lwd=0.5)+
  geom_point(aes(x=Fish.Lower$FDiv, y=Abund.Lower), col="grey")+ylim(-10,50)+
  
  theme_me+
  scale_x_continuous(expand=c(0,0), limits=c(0,1.2)) +
  scale_y_continuous(expand=c(0,0), limits=c(-50,50)) +
  coord_cartesian(xlim=c(0.3,1), ylim=c(0,50))+
  ylab("Abundance")+xlab("Functional Divergence")

#tiff("Figure6.tiff", width = 5, height = 3, units = 'in', res = 1000)
multiplot(Dis, Div, cols=2)
#dev.off()

################################################################################
# Calculate distance between Rainbow Smelt and other species
################################################################################
D.RS<-matrix(0)
for(x in 1:(length(Axes$A1))){
  D.RS[x]<-sqrt((Axes$A1[x] - Axes$A1[46])^2 +
                  (Axes$A2[x] - Axes$A2[46])^2 +
                  (Axes$A3[x] - Axes$A3[46])^2)
}

D.RS <- data.frame(cbind(D.RS, as.character(Fish.traits$COMMONNAME)))
D.RS$D.RS <- as.character(D.RS$D.RS)
D.RS$D.RS <- as.numeric(D.RS$D.RS)

#calculate distance between Rainbow Smelt and other species in RDA
D.RS<-matrix(0)
for(x in 1:(length(vscores$RDA1))){
  D.RS[x]<-sqrt((vscores$RDA1[x] - vscores$RDA1[46])^2 +
                  (vscores$RDA2[x] - vscores$RDA2[46])^2)
}

D.RS <- data.frame(cbind(D.RS, as.character(Fish.traits$COMMONNAME)))
D.RS
D.RS$D.RS <- as.character(D.RS$D.RS)
D.RS$D.RS <- as.numeric(D.RS$D.RS)
