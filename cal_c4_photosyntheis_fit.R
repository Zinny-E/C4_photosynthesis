## C4 example
############################################
############################################
############################################
library(doBy); library(stats); library(dplyr)
#rm(list=ls()); 
#####clear console
#dev.off(); 
#cat("\f") # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

############################################
# load minpack.lm package to fit C4 curves using nlsLM
############################################
library(minpack.lm)
library(Metrics)

##########general inputs for temp and vpr#####
T.ppdk<-c(10,15,20,25,30,35,40) # From Serrano-Romero and Cousins. 2020 New Phytologist
ppdk<-c(13.2,21.6,32.4,48,70.4,98.8,142) # From Serrano-Romero and Cousins. 2020 New Phytologist
ppdk.df<-data.frame(T.ppdk, ppdk); rm(T.ppdk, ppdk)
ppdk.df$T2<-ppdk.df$T.ppdk^2
ppdk.mod<-lm(log(ppdk) ~ T.ppdk + T2, ppdk.df)
mod.sum<-summary(ppdk.mod)
intercept<-mod.sum$coefficients[1,1]
CTleaf.slope<-mod.sum$coefficients[2,1]
CTleaf2.slope<-mod.sum$coefficients[3,1]
temp.ppdk.df<-data.frame(intercept, CTleaf.slope, CTleaf2.slope)
Tleaf= mean(aci_c4$Tleaf)


############################################
# load implementation of C4 ACi curve fitter (based on von Caemmerer 2000 and function AciC4 from plantecophys)
############################################
#Km,GammaStar Optionally, provide Michaelis-Menten coefficient for Farquhar model, and Gammastar. If not provided, they are calculated with a built-in function of leaf temperature
#Rd Day respiration rate (mu mol m-2 s-1), optional (if not provided, calculated from Tleaf, Rd0, Q10 and TrefR). Must be a positive value (an error occurs when a negative value is supplied
O2=210 #O2 Mesophyll O2 concentration
FRM=0.5 #FRM Fraction of day respiration that is mesophyll respiration (Rm)
alpha=0 #alpha Fraction of PSII activity in the bundle sheath (-)
Q10=2 #Q10 T-dependence parameter for Michaelis-Menten coefficients
x=0.4 #x Partitioning factor for electron transport
THETA=0.7 #THETA Shape parameter of the non-rectangular hyperbola
low_gammastar <- 1.93e-4 # Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)

#Vpr=80 #Vpr PEP regeneration (mu mol m-2 s-1)
gbs= 3e-3 #gbs Bundle sheath conductance (mol m-2 s-1)



######enzyme-limited photosynthesis functin#####
# enzyme-limited photosynthesis function
A.enzyme.func=a ~ (-(-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
                             K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko)))) - 
                         sqrt((-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
                                   K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko))))^2 - 4 * a.c * 
                                ((Vcmax - Rd) * ((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) - 
                                   (Vcmax * gbs * low_gammastar * O2 + Rd * gbs * K))))/(2 * a.c)-Rd




#######enzyme-limited predicted photosynthetic function ########
Pre_enz<- function(ci)  { 
  (-(-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
         K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko)))) - 
     sqrt((-(((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) + (Vcmax - Rd) + gbs * 
               K +((alpha/0.047) * (low_gammastar * Vcmax + Rd * Kc/Ko))))^2 - 4 * a.c * 
            ((Vcmax - Rd) * ((pmin(ci * Vpmax/(ci + Kp), Vpr)) - Rm + gbs * ci) - 
               (Vcmax * gbs * low_gammastar * O2 + Rd * gbs * K))))/(2 * a.c)-Rd

}
  
  

######## light-limited photosynthesis function ##########
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


#########light-predicted photosynthesis function ######
pre_light<-function(ci)  {
  (-(-((x * ((1/(2 * THETA)) * 
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
}





############################################
# set parameters
############################################
id_c4 <- 'KBS_Zmay_1' # identifier for C4 curve of interest
b_tresp_R <- 0.1781075 # value for b parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)
c_tresp_R <- -0.00179152 # value for c parameter in temperature response curve for respiration (From Smith and Dukes (2017), doi = 10.1111/gcb.13735)

############################################
# read in ACi curve
############################################
source<- read.csv("~/git_repo/C4_photosynthesis/data/A_Ci_rawdata.csv")# set path within local environment
reps<-unique(source$curve_no)
source<-subset(source, KEEP == "Y")
source<-source[order(source$ci),]

no<- 75

aci_c4<-subset(source, curve_no == reps[[no]])
#head(aci_c4)

# set fit parameters (see documentation for function AciC4 from plantecophys)
############################################
Tleaf=mean(aci_c4$Tleaf) # in dataframe, in CHR format # mean of measurement_t = 24.7
Rd= 1 #find that article on this. (Wang et al., 2014 ; experimental botany; Three distinct biochemical subtypes of C4 photosynthesis)
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
Vpr<-exp(temp.ppdk.df$intercept+(temp.ppdk.df$CTleaf.slope*Tleaf)+(temp.ppdk.df$CTleaf2.slope*Tleaf**2))
############################################
# plot data to estimate transition point & find irregularities
############################################
ci_trans <- 150  # estimated Ci transition point


#plot(aci_c4$a~aci_c4$ci)
#points(aci_c4$a ~ aci_c4$ci, type = "l")
#abline(v = ci_trans)
# note the need to remove Ci point below 0
############################################
# fit the enzyme and light limited portions of the curve separately
############################################
#est.Vcmax<-30
fit_enzyme <- nlsLM(A.enzyme.func, data = subset(aci_c4, ci < ci_trans & ci > 0), 
                    start=list(Vcmax=10,Vpmax=100), control=nls.control(maxiter=500, minFactor=1/10000)) # fit the enzyme limited portion of the curve (Vcmax and Vpmax)



fit_light = nlsLM(A.light.func, data= subset(aci_c4, ci >= ci_trans), 
                  start = list(Jmax = 130), control = nls.control(maxiter = 500, minFactor = 1/10000)) # fit the light-limted portion of the curve (Jmax)




######Summary of enzyme and light####
sum.enzyme<-summary(fit_enzyme)
sum.light<-summary(fit_light)

Vcmax<-sum.enzyme$coefficients[1,1]
Vpmax<-sum.enzyme$coefficients[2,1]
Jmax<-sum.light$coefficients[1,3]


####predicted and observed enzyme #####
Obser_enz<-subset(aci_c4, ci < ci_trans & ci > 0)
enz_ci<-Obser_enz$ci
A_enz_obs<-Obser_enz$a
A_enz_pre<-Pre_enz(enz_ci)
rmse.enz<-rmse(A_enz_obs,A_enz_pre)
A_obs_enz<- data.frame(A_enz_obs, A_enz_pre,no)



#####predicted and observed light ####
Obser_light<-subset(aci_c4, ci >=ci_trans)
light_ci_<-Obser_light$ci
A_light_obs<-Obser_light$a
A_light_pre <- pre_light(light_ci_)
rmse.light<-rmse(A_light_obs, A_light_pre)
A_obs_light<- data.frame(A_light_obs, A_light_pre,no)



                 

#####Calculate Residual ####
#fit model
#model_enz <- lm(A_enz_pre ~ A_enz_obs, data=A_obs_enz)
#resid(model_enz)

#view model summary
#summary(model_enz)

#qqnorm(resid(model)) # A quantile normal plot - good for checking normality
#qqline(resid(model))

#plot(density(resid(model_enz)))

#Stand_re<- rstandard(model)

#final_data <- cbind(A_obs_enz, Stand_re)
#plot predictor variable vs. standardized residuals
#plot(final_data$A_enz_pre, Stand_re, ylab='Standardized Residuals', xlab='A_enz_pre') + abline(0,0)



#Obser
#enz.pr<-predict(fit_enzyme)
#enz.ob<-subset(aci_c4, ci < ci_trans & ci > 0)
#enz.obs<-enz.ob$a

#enz.ob_pr.df<-data.frame(enz.pr, enz.obs, no)

#rmse(enz.obs, enz.pr)


####### Output dataframe ############
#output<-data.frame(Vcmax=as.numeric(),
#                    Vpmax=as.numeric(),               Jmax=as.numeric(),
#          ci_trans=as.numeric(),
#             no=as.numeric(),
 #           Tleaf=as.numeric(),
 #              Vpr=as.numeric(),
 #            Kc=as.numeric(), 
 #             Ko=as.numeric(), Kp=as.numeric(),
 #            rmse.enz=as.numeric(), rmse.light=as.numeric())


########Enzyme limited dataframe #####
output_enz<-data.frame(A_enz_obs=as.numeric(),
                    A_enz_pre=as.numeric(), no=as.numeric())
 


##########Light limited dataframe #####
output_light<-data.frame(A_light_obs=as.numeric(),
                       A_light_pre=as.numeric(), no=as.numeric())



##### Plot graph ######
par(mfrow = c(2,2))
plot(aci_c4$a~aci_c4$ci) + points(aci_c4$a ~ aci_c4$ci, type = "l") +  abline(v = ci_trans) + title("A/ci")
plot(A_obs_enz$A_enz_pre ~ A_obs_enz$A_enz_obs) + abline(1,1) + title("enzymeLimited")
plot(A_obs_light$A_light_pre ~ A_obs_light$A_light_obs) + abline(1,1) + title("LightLimited")





# make a temp df for storage ####
temp.df<-data.frame(Vcmax, Vpmax,Jmax, ci_trans, no,Tleaf,Vpr,Kc,Ko,Kp,rmse.enz,rmse.light)
temp.df

#enz.ob_pr.df


#output<-rbind(output,temp.df)
#output_enz<-rbind(output_enz, A_obs_enz)
#output_light<-rbind(output_light, A_obs_light)




  

## Write csv file
#aci_output<-output
#write.csv (aci_output, "~/git_repo/C4_photosynthesis/data/output2.csv")




