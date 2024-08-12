library(boot)
library(dplyr)
library(cuRe)
source("functions.R")

# Compute r_n by bootstrap
compute_rn_boot <- function(data, indices) {
  data_boot <- data[indices, ]
  cm_boot <- fit.cure.model(Surv(tte,evt)~group, formula.surv=list(~group,~1),data=data_boot,dist = "weibull",link="logit")
  results <- compute_ratio_test(data_boot, cm_boot)
  rn_a <- results[["rn_a"]]
  rn_b <- results[["rn_b"]]
  
  cm_boot <- fit.cure.model(Surv(tte,evt)~1, data=data, dist = "weibull",link="logit")
  
  whichTau <- max(data$tte)
  pi <- exp(unlist(cm_boot$coefs)[1]) / (1 + exp(unlist(cm_boot$coefs)[1]))
  scale <- exp(unlist(cm_boot$coefs)[2])
  shape <- exp(unlist(cm_boot$coefs)[3])
  
  rn_global <- compute_rn(whichTau, shape, scale, pi)
  return(c(rn_a, rn_b, rn_global))
}

#------------------------QALY---------------------------------------------------------------------------------
#Import des données
qaly <- read.csv("graaphR_ensai.csv", 
                 header = TRUE, 
                 sep = ";", 
                 encoding = "latin1", na.strings = c("", " ", "NA", "NI"))

qaly$death.dt <- as.Date(as.character(qaly$deathdt),format="%Y-%m-%d")
qaly$rando.dt <- as.Date(as.character(qaly$randodt),format="%Y-%m-%d")
qaly$max.dt <- as.Date(as.character(qaly$Datemax),format="%Y-%m-%d")
qaly$rech.dt <- as.Date(as.character(qaly$relapsdt),format="%Y-%m-%d")
qaly$suivi  <- as.numeric(qaly$max.dt-qaly$rando.dt)/ (365.25/12)
qaly$delrech<- as.numeric(qaly$rech.dt-qaly$rando.dt)/ (365.25/12)
qaly$delpfs<- qaly$suivi
qaly$delpfs[!is.na(qaly$delrech)]<-qaly$delrech[!is.na(qaly$delrech)]
qaly$R1<-factor(qaly$R1,levels=c("Intensive arm (A)","Light arm (B)"))

data<- qaly %>% select(delpfs,pfs,R1) %>% rename(tte=delpfs,evt=pfs,group=R1) %>% mutate(group=(as.numeric(as.factor(group))-1))
boot_res <- boot(data=data, statistic=compute_rn_boot, R=1000)
saveRDS(boot_res, "boot_res_qaly.rds")

#### Inversion des bras
qaly$R1inv <- ifelse(qaly$R1=="Intensive arm (A)",1,0)
datainv<- qaly %>% select(delpfs,pfs,R1inv) %>% rename(tte=delpfs,evt=pfs,group=R1inv)
boot_res_inv <- boot(data=datainv, statistic=compute_rn_boot, R=1000)
saveRDS(boot_res_inv, "boot_res_qaly_inv.rds")

#------------------------COVID---------------------------------------------------------------------------------
#Import des données
covid <- read.csv("covidicus1.csv", 
                  header = TRUE, 
                  sep = ";", 
                  encoding = "latin1", na.strings = c("", " ", "NA", "NI"))


covid$rando.dt <- as.Date(as.character(covid$randodt),format="%Y-%m-%d")
covid$max.dt <- as.Date(as.character(covid$Dader),format="%Y-%m-%d")
cens_dc_idx<-which(covid$dc == 1 & covid$sor.dt < covid$max.dt)
covid$dc[cens_dc_idx]<-0
covid$survie  <- as.numeric(covid$max.dt-covid$rando.dt)
covid$survie[which(covid$survie>60)] <- 60

data<- covid %>% select(survie,dc,DXM) %>% rename(tte=survie,evt=dc,group=DXM) %>% mutate(group=(as.numeric(as.factor(group))-1))
boot_res <- boot(data=data, statistic=compute_rn_boot, R=1000)
saveRDS(boot_res, "boot_res_covid.rds")

#### Inversion des bras
covid$DXMinv <- ifelse(covid$DXM=="A",1,0)
datainv<- covid %>% select(survie,dc,DXMinv) %>% rename(tte=survie,evt=dc,group=DXMinv)
boot_res_inv <- boot(data=datainv, statistic=compute_rn_boot, R=1000)
saveRDS(boot_res_inv, "boot_res_covid_inv.rds")
