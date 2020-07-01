#NB: this is the script used for article by Lutz et al. (june 2020).
#Script is written bij M. C. Lutz-Landesbergen, MSc, miranda.c.lutz@gmail.com
#set working directory to appropiate folder.
#paths "C:\Users\mclut\OneDrive\Documenten\R\ERN_PE_meta\Final_analysis_r1."
#check package version through, check updates: 
packageVersion("pkg")

#Load packages:
library(metafor)
library(readxl)
library(rlang)
library(grid)
library(meta)
library(dmetar)

#load data.

#outline of script, five sections:(1) calculate EF, (2) Overall Model, 
#(3) inspection/heterogeneity, (4) moderation/subgroup analysis, 
#(5) small study bias assessment

# Calculate EF ------------------------------------------------------------
 
#note, meta and metafor are used interchangeably. There are slight difference in outcomes but...
#but these did not influence results.

#Rename data.
EF_ERN_FCz<-ERN_FCz_tot1

###
#(1) calculate EF using escalc.
EF_ERN_FCz<-escalc(measure = "SMD", m1i = Me, sd1i = Se, m2i = Mc, sd2i = Sc, n1i = Ne, n2i = Nc, 
                   vtype = "UB", data = ERN_FCz_tot1, append = T)
# we did not use vtype=LS due to lack of large sample, but when running, smd and variances identical.
EF_ERN_FCz
options(max.print = 99999)
print.escalc(EF_ERN_FCz)
#Yi is SMD; vi is variance in effect size.
#in case you want to use multilevel analysis with other electrodes, you can use this function to inspect the variances.
###

# Overall Model ------------------------------------------------------------

###
# (2) ma in metafor: default estimator method is REML.
ma_ern <-rma(yi, vi, ni, data = EF_ERN_FCz, test = "knha")
#show summary with model fit. See notes print.rma {metafor} for formula.
summary(ma_ern)

#Robust function, checking for first authorship. Results, of 27, 21 clusters -therefore futile.
ma_ernR <- robust(ma_ern, cluster = EF_ERN_FCz$FirstA)
ma_ernR

#find outliers.
find.outliers(ma_ern)
#Study 24, Sokhadze, excl, lowers cEF from .43 to .40, drops hetergeneity with almost 9% making it insign.

#check estimator effect; DerSimonian-Laird estimator and Sidik-Jonkman estimator.
ma_ernDL <-rma(yi, vi, ni, data = EF_ERN_FCz, method = "DL", test = "knha")
summary(ma_ernDL)
ma_ernSJ <-rma(yi, vi, ni, data = EF_ERN_FCz, method = "SJ", test = "knha")
summary(ma_ernSJ)
#assessing the model information criteria: REML and DL suits best. SJ not.

#MA using meta
ern_fcz<-ERN_FCz_tot1
ma_ern_fcz<-metacont(Ne, Me, Se, Nc, Mc, Sc, data = ern_fcz, 
                     studlab = paste(Authors), comb.fixed = FALSE, comb.random = TRUE, 
                     method.tau = "REML", hakn = TRUE, prediction = TRUE, sm = "SMD")
ma_ern_fcz
forest(ma_ern_fcz)

#forestplot of MA ERN, exporting figure2.
par(mar=c(0,0,0,0))

#forestplot met meta.
forest(ma_ern_fcz)
forest(ma_ern_fcz, xlim = c(-1,2),   col.predict = "black",  digits.sd = 2, fontsize = 10)
grid.text("Increased ERN          Reduced ERN", x = unit(0.6, "npc"), y = unit(0.05, "npc"), gp= gpar(fontface="bold"))
grid.lines(x = unit(c(0.03, 0.97), "npc"),
           y = unit(c(0.85), "npc"))

#tiff-en.
tiff(filename = "<path>ERN_forest.tiff", units = "in", width = 10, height = 8, res = 300)
forest(ma_ern_fcz, xlim = c(-1,2),   col.predict = "black", digits.sd = 2, fontsize = 10)
grid.text("Increased ERN          Reduced ERN", x = unit(0.57, "npc"), y = unit(0.038, "npc"), gp= gpar(fontface="bold", fontsize = 10))
grid.lines(x = unit(c(0.07, 0.93), "npc"),
           y = unit(c(0.86), "npc"))
dev.off()
###

# Inspection/heterogeneity ------------------------------------------------


###
# (3) Inspection Heterogeneity. 
#influence analysis, if available through library(dmetar). 
infl_ern<-influence(ma_ern)
infl_ern
plot(infl_ern)
#identify influential studies by looking at 'outliers'.Sokhadze keeps coming, and 
#Michelini is an infuential case due to high sample size

#leave one out analysis tests what will happen to all meta information when one study is excluded.
#Also helps to identify influential cases.
leave1out(ma_ern)

#investigating other plots:labbe(ma_ern)
radial(ma_ern)

baujat(ma_ern)

qqnorm(ma_ern)

#[other plots:
plot.rma.uni(ma_ern)
profile.rma.uni(ma_ern)
#]

#note, we can do analysis with and without study 10 and 24.
###

# moderation/subgroup analysis --------------------------------------------


###
#(4) moderation analysis
#Moderators tested apriori in manuscript are
# diagnosis (variable Diagnosis3), comorbidity and ExpParad.
#using meta-package (due to lack of vectors in columns).
#see meta-script.

#Diagnosis as moderator.
diagnosis3.subgroup<-update.meta(ma_ern_fcz, 
                                 byvar=Diagnosis3, 
                                 comb.random = TRUE, 
                                 comb.fixed = FALSE)
diagnosis3.subgroup
forest(diagnosis3.subgroup)

#Diagnosis4 as moderator, used for revision, not in manuscript.
diagnosis4.subgroup<-update.meta(ma_ern_fcz, 
                                 byvar=Diagnosis4, 
                                 comb.random = TRUE, 
                                 comb.fixed = FALSE)
diagnosis4.subgroup
forest(diagnosis4.subgroup)

#Comorbidity as a moderator. 
comor.subgroup<-update.meta(ma_ern_fcz, 
                            byvar=ComCat, 
                            comb.random = TRUE, 
                            comb.fixed = FALSE)
comor.subgroup
forest(comor.subgroup)

#Feedback in task as a moderator
fbg.subgroup<-update.meta(ma_ern_fcz, 
                          byvar=FeedbackGiven, 
                          comb.random = TRUE, 
                          comb.fixed = FALSE)
fbg.subgroup
forest(fbg.subgroup)

fbt.subgroup<-update.meta(ma_ern_fcz, 
                          byvar=Task_FB, 
                          comb.random = TRUE, 
                          comb.fixed = FALSE)
fbt.subgroup
forest(fbt.subgroup)

#Experimental paradigm as a moderator.
task.subgroup<-update.meta(ma_ern_fcz, 
                           byvar=ExpParad, 
                           comb.random = TRUE, 
                           comb.fixed = FALSE)
task.subgroup
forest(task.subgroup)

#Medication as a moderator.#not included in manuscript.
med1.subgroup<-update.meta(ma_ern_fcz, 
                           byvar = Medication1, 
                           comb.random = TRUE, 
                           comb.fixed = FALSE)
med1.subgroup
forest(med1.subgroup)

med2.subgroup<-update.meta(ma_ern_fcz, 
                           byvar = Medication2, 
                           comb.random = TRUE, 
                           comb.fixed = FALSE)
med2.subgroup
forest(med2.subgroup)

#Task adjustment as a moderator.#not included in manuscript.
ta.subgroup<-update.meta(ma_ern_fcz, 
                         byvar = AdjustedParadigm, 
                         comb.random = TRUE, 
                         comb.fixed = FALSE)
ta.subgroup
forest(ta.subgroup)

#metareg for age.
#note: t test for means age between exp and control group:
#The t-value is 0.4557. The p-value is .325252. The result is not significant at p < .05
MRage<-metareg(ma_ern_fcz, AvAge, method.tau = "REML", intercept = T)
MRage
bubble(MRage, xlab = "AvAge", col.line = "blue", studlab = T)
###

#  small study bias assessment --------------------------------------------


###
#(5) Publication bias assessment.
#using meta package.
#Eggers test: first without trim and fill. 
metabias(ma_ern_fcz, method.bias = "linreg")
# then funnel plot. 
funnel1<-funnel(ma_ern, legend =TRUE)
funnel1

#color enhanced 
funnel2<-funnel(ma_ern, xlab="Hedges' g", 
                contour = c(.95,.975,.99),
                col.contour=c("darkgrey","grey","lightgrey"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("darkgrey","grey","lightgrey"))

#fix warnings!

# Eggers test for the intercept from dmetar.
eggers.test(ma_ern_fcz)
funnel(ma_ern_fcz)
#Now with trim and fill procedure. 
ma_ern_fcz.trimfill<-trimfill(ma_ern_fcz)
ma_ern_fcz.trimfill
#7 studies added, to get to .29. 
funnel(ma_ern_fcz.trimfill, xlab = "Hedges g")
summary(ma_ern_fcz.trimfill)
eggers.test(ma_ern_fcz.trimfill)

#FIX this. Then we test everything again but 
#with adjusted funnel plot to investigate small effects.
reg<- regtest(ma_ern, predictor = "sqrtninv")
reg
#this becomes insignificant, meaning no assmmetry?
funnel3 <-funnel(ma_ern, ylab = "sqrtninv")
funnel3

#Funnel plot of MA ERN, exporting figure3. #check path! 
funnel(ma_ern_fcz.trimfill, xlab="Hedges' g", 
       contour = c(.95,.975,.99),
       col.contour=c("grey20","grey50","grey100"))+
  legend(1.4, 0, c("p < 0.05", "p<0.025", "< 0.01"),bty = "n",
         fill=c("grey20","grey50","grey100"))

library(ggplot2)

#export figure.
par(mar=c(0,0,0,0))
tiff(filename = "C:/Users/mclut/OneDrive/Documenten/R/ERN_PE_meta/Final_analysis_r1/Exporting_figures/ERN_funnel.tiff", 
     units = "in", width = 7, height = 5, res = 300)
funnel(ma_ern_fcz.trimfill, xlab="Hedges' g", ylab = "sqrtninv", xlim = c(-1.2, 1.6),
       contour = c(.95,.975,.99),
       col.contour=c("grey20","grey50","grey100"))+
  coord_cartesian(ylim = c(0.5, 0)) +
legend(1, 0, c("p < 0.05", "p < 0.025", "p < 0.01"),bty = "n",
         fill=c("grey20","grey50","grey100")) 
  
dev.off()

#P-curve analysis, get code from book/package dmetar. 
pcurve(ma_ern_fcz)
#see word document for interpretation.
#calculation of 'true' effect size, which is .316.
n.ma_ern_fcz<-ern_fcz$ni
n.ma_ern_fcz
pcurve(ma_ern_fcz, effect.estimation = TRUE, N = n.ma_ern_fcz, dmin = 0, dmax = 1)

#fail safe-N assessment: does not indicate bias, indicate strength of effect.
fsn(yi, vi, data = EF_ERN_FCz, type = "Rosenthal")
fsn(yi, vi, data = EF_ERN_FCz, type = "Orwin")
fsn(yi, vi, data = EF_ERN_FCz, type = "Rosenberg")
#Rosenthal: 603 studies needed to make the ma not sign  
#Rosenberg: 371 studies need to make null results possible.
#Orwin n =27, meaning we need 27 effect sizes to obtain a combined effect size of .2363 (unweighted av effect size).