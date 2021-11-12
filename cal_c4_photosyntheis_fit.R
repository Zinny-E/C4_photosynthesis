## C4 example
############################################
############################################
############################################
library(doBy); library(stats)
# rm(list=ls()); 
#clear console
#dev.off(); 
#cat("\f") # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

############################################
# load minpack.lm package to fit C4 curves using nlsLM
############################################
library(minpack.lm)
############################################
# load implementation of C4 ACi curve fitter (based on von Caemmerer 2000 and function AciC4 from plantecophys)
############################################
#general inputs

#Km,GammaStar Optionally, provide Michaelis-Menten coefficient for Farquhar model, and Gammastar. If not provided, they are calculated with a built-in function of leaf temperature

#Rd Day respiration rate (mu mol m-2 s-1), optional (if not provided, calculated from Tleaf, Rd0, Q10 and TrefR). Must be a positive value (an error occurs when a negative value is supplied

O2=210 #O2 Mesophyll O2 concentration
FRM=0.5 #FRM Fraction of day respiration that is mesophyll respiration (Rm)
alpha=0 #alpha Fraction of PSII activity in the bundle sheath (-)
Q10=2 #Q10 T-dependence parameter for Michaelis-Menten coefficients
x=0.4 #x Partitioning factor for electron transport
THETA=0.7 #THETA Shape parameter of the non-rectangular hyperbola
low_gammastar <- 1.93e-4 # Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)
Vpr=80 #Vpr PEP regeneration (mu mol m-2 s-1)
gbs= 3e-3 #gbs Bundle sheath conductance (mol m-2 s-1)

# enzyme-limited photosynthesis function
A.enzyme.func=a ~ (-(-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
                             K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko)))) - 
                         sqrt((-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
                                   K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko))))^2 - 4 * a.c * 
                                ((Vcmax - Rd) * ((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) - 
                                   (Vcmax * gbs * low_gammastar * O2 + Rd * gbs * K))))/(2 * a.c)-Rd

# light-limited photosynthesis function
A.light.func= a~ (-(-((x * ((1/(2 * THETA)) * 
                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                             Rm + gbs * ci) + ((1 - x) * ((1/(2 * THETA)) * 
                                                            (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 - Rd) +
                            gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 *((1 - x) * 
                                                                                               ((1/(2 * THETA)) * 
                                                                                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7*Rd/3))) - 
                        sqrt((-((x * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                                   Rm + gbs * ci) + ((1 - x) * ((1/(2 * THETA)) * 
                                                                  (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 - Rd) +
                                  gbs * (7 * low_gammastar * O2/3) + alpha * low_gammastar/0.047 *((1 - x) * 
                                                                                                     ((1/(2 * THETA)) * 
                                                                                                        (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7*Rd/3)))^2 - 
                               4 * a.j * (((x * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/2 - 
                                              Rm + gbs * ci) * ((1 - x) * ((1/(2 * THETA)) * (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * 
                                                                                                                  THETA * Qp2 * Jmax)))/3 - Rd)) -gbs * 
                                            low_gammastar * O2 * ((1 - x) * ((1/(2 * THETA)) * 
                                                                               (Qp2 + Jmax - sqrt((Qp2 + Jmax)^2 - 4 * THETA * Qp2 * Jmax)))/3 + 7 * Rd/3))))/(2 * a.j)-Rd

############################################
# set parameters
############################################
id_c4 <- 'KBS_Zmay_1' # identifier for C4 curve of interest
b_tresp_R <- 0.1781075 # value for b parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
c_tresp_R <- -0.00179152 # value for c parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)

############################################
# read in ACi curve
############################################
source<- read.csv("~/git_repo/C4_photosynthesis/data/A_Ci_rawdata.csv") # set path within local environment
reps<-unique(source$curve_no)
aci_c4<-subset(source, curve_no == reps[[6]])
#head(aci_c4)
############################################
# read in SLCE_data
############################################'
#slce <- read.csv('LCE_data.csv') # set path to local environment
#slce_c4 = subset(slce, aci_id == id_c4) # slce data for individual of interest

############################################
# extract dark respiration value from SLCE_data and add to Aci dataset
############################################
#rd_c4 <- slce_c4$Rd # find Rd value for curve of interest
#rd_c4_tphoto <- rd_c4 / exp((b_tresp_R * (slce_c4$Tleaf_R - slce_c4$Tleaf_photo)) + 
                     #         (c_tresp_R * (slce_c4$Tleaf_R^2 - slce_c4$Tleaf_photo^2))) # adjust Rd to similar temperature as the A/Ci data
#aci_c4$Rd <- rd_c4_tphoto # add dark respiration to ACi data

############################################
# set fit parameters (see documentation for function AciC4 from plantecophys)
############################################
Tleaf=25 # in dataframe, in CHR format # mean of measurement_t = 24.7
Rd= 1 #find that article on this. 
PPFD=mean(aci_c4$measurement_ppfd)
# Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and O (Ko, mmol mol-1) and combined (K)
Kc <- 650*Q10^((Tleaf-25)/10)
Kp <- 80*Q10^((Tleaf-25)/10)
Ko <- 450*Q10^((Tleaf-25)/10)
K <- Kc*(1+O2/Ko)
Rm <-  FRM*Rd # Day leaf respiration, umol m-2 s-1
Qp2 <- PPFD*0.85*(1-0.15)/2 # Non-rectangular hyperbola describing light effect on electron transport rate (J)
a.c <- 1 - (alpha*Kc)/(0.047*Ko)
a.j <- 1 - 7 * low_gammastar * alpha/(3 * 0.047)

############################################
# plot data to estimate transition point & find irregularities
############################################

plot(aci_c4$a~aci_c4$ci)
ci_trans <- 50 # estimated Ci transition point
abline(v = ci_trans)
# note the need to remove Ci point below 0

############################################
# fit the enzyme and light limited portions of the curve separately
############################################
#est.Vcmax<-30

fit_enzyme <- nlsLM(A.enzyme.func, data = subset(aci_c4, ci < ci_trans & ci > 0), 
                    start=list(Vcmax=10,Vpmax=100), control=nls.control(maxiter=500, minFactor=1/10000)) # fit the enzyme limited portion of the curve (Vcmax and Vpmax)
fit_light = nlsLM(A.light.func, data= subset(aci_c4, ci >= ci_trans), 
                  start = list(Jmax = 130), control = nls.control(maxiter = 500, minFactor = 1/10000)) # fit the light-limted portion of the curve (Jmax)


# make a temp df for storage ####
#Obser<-subset(aci_c4, ci < ci_trans & ci > 0)
#Obser
#predict(fit_enzyme)


output<-data.frame(Vcmax=as.numeric(),
                    Vpmax=as.numeric(), 
                   Jmax=as.numeric(),
                   ci_trans=as.numeric())



sum.enzyme<-summary(fit_enzyme)
sum.light<-summary(fit_light)

Vcmax<-sum.enzyme$coefficients[1,1]
Vpmax<-sum.enzyme$coefficients[2,1]
Jmax<-sum.light$coefficients[1,3]

temp.df<-data.frame(Vcmax, Vpmax,Jmax, ci_trans)

#output<-rbind(output, temp.df)
