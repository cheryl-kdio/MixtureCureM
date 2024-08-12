#--------------------------------QALY DATA ANALYSIS-------------------------------------------------------------------------
library(survRM2)
library(survival)
library(survminer)
library(cuRe)
library(dplyr)
source("functions1.R")
source("ggcompetingrisks1.R")


#Import
qaly <- read.csv("GRAAPH/graaphR_ensai.csv", 
                 header = TRUE, 
                 sep = ";", 
                 encoding = "latin1", na.strings = c("", " ", "NA", "NI"))

#Formattage en date
qaly$death.dt <- as.Date(as.character(qaly$deathdt),format="%Y-%m-%d")
qaly$rando.dt <- as.Date(as.character(qaly$randodt),format="%Y-%m-%d")
qaly$max.dt <- as.Date(as.character(qaly$Datemax),format="%Y-%m-%d")
qaly$rech.dt <- as.Date(as.character(qaly$relapsdt),format="%Y-%m-%d")



# Délais calculés
qaly$suivi  <- as.numeric(qaly$max.dt-qaly$rando.dt)/ (365.25/12)
qaly$delrech<- as.numeric(qaly$rech.dt-qaly$rando.dt)/ (365.25/12)
qaly$delpfs<- qaly$suivi
qaly$delpfs[!is.na(qaly$delrech)]<-qaly$delrech[!is.na(qaly$delrech)]
qaly$R1<-factor(qaly$R1,levels=c("Intensive arm (A)","Light arm (B)"))


#RMST
rmst2(qaly$delpfs, qaly$pfs, as.numeric(qaly$R1)-1, covariates = NULL, alpha = 0.05)
#---------------------------------------------------------------------------------------------------------
#-------- Estimation des courbes de survie KM
efs <- survfit(Surv(delpfs,pfs)~R1, data=qaly)

efs_p<-ggsurvplot(efs, data=qaly,
                  risk.table = TRUE,
                  ggtheme = theme_test(),
                  conf.int = F,
                  xlab = "Time from randomization (Months)",
                  ylab = "Pr(RFS)"
)

efs_p$plot+
  coord_cartesian(xlim=c(0,80),ylim=c(0,1),expand = F)+
  scale_color_manual(values=c("#009E73","#56B4E9"))


# COX model
coxph_model<-coxph(Surv(delpfs,pfs)~R1, data=qaly)
summary(coxph_model)

# Log rank
survdiff(Surv(delpfs,pfs)~R1, data=qaly)

#-------- Rapport de risques instantanés
x<-coxph_model<-coxph(Surv(delpfs,pfs)~R1, data=qaly)
summary(coxph_model)

# Test de risques proportionnels
t<-cox.zph(coxph_model)
ggcoxzph(t)


ggsurvplot(survfit(x, newdata = data.frame(R1 = c("Intensive arm (A)","Light arm (B)"))),data=qaly,conf.int = F)

#---------------------------------------------------------------------------------------------------------
#-------- Modèle de cure de mélange (pkg cuRe)
cure_model<-fit.cure.model(Surv(delpfs,pfs)~R1,formula.surv = list(~R1,~1),
                           data=qaly, dist = "weibull",link="logit")
summary(cure_model)

#-------- Estimations
### Probability of being cured 
pi_u<-c(1, 1, 0, 0, 0) %*% unlist(cure_model$coefs)
exp(pi_u)
exp(pi_u)/(1+exp(pi_u))

#OR and HR
exp(unlist(cure_model$coefs)[c(2,4)])

### IC
var_u<-c(1, 1, 0, 0, 0) %*% cure_model$covariance %*% c(1, 1, 0, 0, 0)
exp(pi_u- 1.96*sqrt(var_u))/(1+exp(pi_u- 1.96*sqrt(var_u)))
exp(pi_u+ 1.96*sqrt(var_u))/(1+exp(pi_u+ 1.96*sqrt(var_u)))

#---------------------------------------------------------------------------------------------------------
#-------- Ratio tests by arm

data<- qaly %>% select(delpfs,pfs,R1) %>% rename(tte=delpfs,evt=pfs,group=R1) %>% mutate(group=(as.numeric(as.factor(group))-1))
results<-compute_ratio_test(data, cure_model)
unlist(results)*100

#------------- OVERALL
cm_overall <- fit.cure.model(Surv(delpfs,pfs)~1, data=qaly, dist = "weibull",link="logit")
whichTau <- max(qaly$delpfs)
pi <- exp(unlist(cm_overall$coefs)[1]) / (1 + exp(unlist(cm_overall$coefs)[1]))
scale <- exp(unlist(cm_overall$coefs)[2])
shape <- exp(unlist(cm_overall$coefs)[3])

r_n <- compute_rn(whichTau, shape, scale, pi)

#----------------------- MALLER ZHOU
data<- qaly %>% select(delpfs,pfs,R1) %>% rename(tte=delpfs,evt=pfs,group=R1) %>% mutate(group=(as.numeric(as.factor(group))-1))

cm <- fit.cure.model(Surv(tte,evt)~group, formula.surv=list(~group,~1),data=data,dist = "weibull",link="logit")

alphatest(data[data$group==0,])
alphatest(data[data$group==1,])
alphatest(data)

deviance_test(data,cm,0.95)
cm <- fit.cure.model(Surv(tte,evt)~1, formula.surv=list(~1,~1),data=data,dist = "weibull",link="logit")
deviance_test(data,cm,0.95)


#----------------------------------COVID DATA ANALYSIS-----------------------------------------------------------------------

#Import
covid <- read.csv("COVIDICUS/covidicus1.csv", 
                  header = TRUE, 
                  sep = ";", 
                  encoding = "latin1", na.strings = c("", " ", "NA", "NI"))


covid$rando.dt <- as.Date(as.character(covid$randodt),format="%Y-%m-%d")
covid$max.dt <- as.Date(as.character(covid$Dader),format="%Y-%m-%d")

cens_dc_idx<-which(covid$dc == 1 & covid$sor.dt < covid$max.dt)
covid$dc[cens_dc_idx]<-0
covid$survie  <- as.numeric(covid$max.dt-covid$rando.dt)
covid$survie[which(covid$survie>60)] <- 60

#--------------------------------- Analyse de survie
os<-survfit(Surv(survie,dc)~DXM, data=covid)
os_p<-ggsurvplot(os, ggtheme = theme_test(),
                 xlab = "Durée de suivi (jours)",
                 ylab = "Probabilité de survie")

os_p$plot+
  coord_cartesian(xlim=c(0,60),ylim=c(0,1),expand = F)+
  scale_color_manual(values=c("#009E73","#56B4E9"))


rmst2(covid$survie,covid$dc,as.numeric(as.factor(covid$DXM))-1,tau = 60)

#----------- Hazard ratio
cox<-coxph(Surv(survie,dc)~DXM,data=covid)
ggsurvplot(survfit(cox, newdata = data.frame(DXM = c(unique(covid$DXM)))),data=covid,conf.int = F)

#---- Test hypothèse de PH
t<-cox.zph(coxph(Surv(survie,dc)~DXM,data=covid))
ggcoxzph(t)

#----------- Modèle de cure
cm<-fit.cure.model(Surv(survie,dc)~DXM,data=covid,formula.surv =list(~DXM,~1) , dist = "weibull",link="logit")
summary(cm)

pi<-c(1, 0, 0, 0, 0) %*% unlist(cm$coefs)
var<-c(1, 0, 0, 0, 0) %*% cm$covariance %*% c(1, 0 , 0, 0, 0)

round(exp(pi- 1.96*sqrt(var))/(1+exp(pi- 1.96*sqrt(var))),2)
round(exp(pi+ 1.96*sqrt(var))/(1+exp(pi+ 1.96*sqrt(var))),2)

#----------- Ratio test
data <- covid %>% select(survie,dc,DXM) %>% rename(tte=survie,evt=dc,group=DXM) %>% mutate(group=(as.numeric(as.factor(group))-1))
results<-compute_ratio_test(data, cm)

#------------- OVERALL
cm_overall <- fit.cure.model(Surv(survie,dc)~1, data=covid, dist = "weibull",link="logit")
whichTau <- max(covid$survie)
pi <- exp(unlist(cm_overall$coefs)[1]) / (1 + exp(unlist(cm_overall$coefs)[1]))
scale <- exp(unlist(cm_overall$coefs)[2])
shape <- exp(unlist(cm_overall$coefs)[3])

r_n <- compute_rn(whichTau, shape, scale, pi)

#----------------------- MALLER ZHOU
data<- covid %>% select(survie,dc,DXM) %>% rename(tte=survie,evt=dc,group=DXM) %>% mutate(group=(as.numeric(as.factor(group))-1))

cm <- fit.cure.model(Surv(tte,evt)~group, formula.surv=list(~group,~1),data=data,dist = "weibull",link="logit")

alphatest(data[data$group==0,])
alphatest(data[data$group==1,])
alphatest(data)

deviance_test(data,cm,0.95)
cm <- fit.cure.model(Surv(tte,evt)~1, formula.surv=list(~1,~1),data=data,dist = "weibull",link="logit")
deviance_test(data,cm,0.95)

#-----------------------COMPETING RISK ANALYSIS----------------------------------------------------------------------------------

# Formatage en date
covid$sorIcu.dt<-as.Date(as.character(covid$sorIcudt),format="%Y-%m-%d")
covid$sorHop.dt<-as.Date(as.character(covid$sorHopdt),format="%Y-%m-%d")
covid$sor.dt<-pmax(covid$sorIcu.dt,covid$sorHop.dt,na.rm = T)

#censor people dead outside hospital to alive
cens_dc_idx<-which(covid$dc == 1 & covid$sor.dt < covid$max.dt)
covid$dc[cens_dc_idx]<-0

# Construction des délais
covid$suivi  <- as.numeric(covid$max.dt-covid$rando.dt)
covid$delsor_wd<- as.numeric(covid$sorIcu.dt-covid$rando.dt)*ifelse(covid$dc==1,NA,1)

### Time to events
covid$TTE<-covid$delsor_wd
covid$TTE[is.na(covid$TTE)]<-covid$suivi[is.na(covid$TTE)]


covid$status<-with(covid,ifelse(dc==1,1,ifelse(!is.na(delsor_wd),2,0)))


### Troncature à 60 jours
covid$TTE[which(covid$TTE>60)] <- 60
#covid$status[which(covid$TTE==60)] <- 0


#Vérification
sum(with(covid[which(covid$dc==1),], sor.dt<max.dt),na.rm=T)
#Sur 144 décès, 10 dates de sortie de réanimation/hopital sont inférieurs aux dates de dernier suivi)



#----------------------------------- Fit competingrisk model
library(cmprsk)
cmpfit<-with(covid,cuminc(TTE,status,DXM))


# Extraire la CIF pour l'événement d'intérêt
cif_eventA <- cmpfit$`A 2`
cif_eventB <- cmpfit$`B 2`

# La limite en l'infini est la dernière valeur de la CIF
lmt_infA <- tail(cif_eventA$est, 1)
lmt_infB <- tail(cif_eventB$est, 1)


pvalue1<-paste("Gray's Test: p-value = ",round(cmpfit$Tests[1,2],2))
pvalue2<-paste("Gray's Test: p-value = ",round(cmpfit$Tests[2,2],2))


ggcompetingrisks(cmpfit, xlab = "Time from randomization (Days)", legend = "top",ylab = "Pr( In hospital death)"
                 ,multiple_panels = F, ggtheme=theme_test(),title="") +
  geom_line(size = 0.7) +  
  scale_color_manual(name="", values=c("#FFBA49","#20A39E"), labels=c("In-hospital death", "Discharged alive"))+
  scale_x_continuous(breaks = seq(0, 60, by = 10),limits = c(0,60),expand = c(0,0))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
  scale_linetype_manual(name="", values=c("solid", "dotted"), labels=c("DXM=A", "DXM=B"))+
  theme(legend.text = element_text(size = 10), 
        legend.key.size = unit(1.5, "lines")) +
  annotate("text",x = 50, y = 0.75, label =pvalue1)+
  annotate("text",x = 50, y = 0.30, label = pvalue2)

ggcompetingrisks1(
  cmpfit,                     
  xlab = "Time from randomization (Days)", 
  ylab = "Pr( In hospital death)",
  ylim=c(0, 1), 
  title="",
  legend="top",
  legend.title="Strata", ,
  conf.int = F, 
  type_group="color", 
  event_suppr = 2,
  ggtheme=theme_test()                    
) +
  geom_line(size = 0.7) +  
  scale_color_manual(name="Strata", values=c("#FFBA49","#20A39E"), labels=c("DXM=A", "DXM=B"))+
  scale_x_continuous(breaks = seq(0, 60, by = 10),limits = c(0,60),expand = c(0,0))+
  scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
  theme(legend.text = element_text(size = 10), 
        legend.key.size = unit(1.5, "lines")) 

# Test difference in risks
paste("Gray's Test: p-value = ",round(cmpfit$Tests[1,2],2))


#No significant difference in risk of mortality or discharge alive between the two groups
#------------------------------------- Cause specific hazard
library(riskRegression)
library(prodlim)


CSC(Hist(TTE,status) ~ DXM,data=covid)


#------------------------------------- Subhazard/Fine-gray model
FGR(Hist(TTE,status) ~ DXM,data=covid)

#--------------------------------------------------------------------------------------------------
library(NPHMC)
library(dplyr)

#covidtest<-covid %>%mutate(X=as.numeric(as.factor(covid$DXM))-1) %>% select(survie,dc,X)
graaphtest<-qaly %>%mutate(X=as.numeric(as.factor(qaly$R1))-1) %>% select(delpfs,pfs,X)

# /!\ lambda0=lambda**(-1/gamma)
# /!\ gamma =exp(beta)

NPHMC(power=0.8,alpha = 0.05,accrualtime = 40, followuptime = 84, p=0.5,
      accrualdist = "uniform", hazardratio = 1.06, oddsratio = .46, pi0 = 0.74,
      survdist = "weib", k=exp(0.16),lambda0 = exp(-3.36)**(-1/exp(0.16)))

NPHMC(power=0.80,alpha=0.05,accrualtime=40,followuptime=84,p=0.5,accrualdist="uniform",survdist = "weib",
       data=graaphtest)

