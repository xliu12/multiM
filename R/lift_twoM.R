library(haven)
# library(foreign)
library(tidyverse)
library(mice)
library(mitml)
# library(miceadds)
library(fastDummies)
library(missForest)
options(dplyr.print_max = 100, pillar.print_max = 50, dplyr.width = Inf)

# other examples
# load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/.RData")
setwd("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/LIFT/LIFT Data/PITL")


# Read in data --------------------
PARTICIPATION <- read_sav("cross waves data/PARTICIPATION BY WAVE.SAV")

P1_15_CHILD_ANTISOCIAL <- read_sav("cross waves data/P1-15 CHILD ANTISOCIAL.SAV")
P1_15CBCL_SCALES <- read_sav("cross waves data/P1_15CBCL SCALES.sav")
P1_15CSUB <- read_sav("cross waves data/P1_15CSUB.SAV")

P1_12PARENT_SUBSTANCE <- read_sav("cross waves data/P1_12PARENT SUBSTANCE.SAV")

P1_12FAMILY_STATUS <- read_sav("cross waves data/P1_12FAMILY STATUS.SAV")
TC_OFFENSES_SUMMARY <- read_sav("cross waves data/TC OFFENSES SUMMARY.SAV")

# construct scales measured at each wave
P1SCLE_p1 <- read_sav("SCLE/P1SCLE.D13")
demo_p1 <- read_sav("P1_1-3/P1DEMO.D13")

P3SCLE <- foreign::read.spss("SCLE/P3SCLE.D13", to.data.frame = T)

P4SCLE_p1 <- read_sav("SCLE/P4SCLE.D13")
P5SCLE_p1 <- read_sav("SCLE/P5SCLE.D13")
P6SCLE_p1 <- read_sav("SCLE/p6scle.d13")
P7SCLE_p1 <- read_sav("SCLE/P7SCLE.D13")
P8SCLE_p1 <- read_sav("SCLE/P8SCLE.D13")
P10SCLE_p1 <- read_sav("SCLE/P10SCLE.D13")
P12SCLE_p1 <- read_sav("SCLE/P12SCLE.D13")

# cross wave construct scores
P1_12RISK <- read_sav("TECH/P1_12RISK.sav") # number of antisocial behavior risk factors rated by teachers
P1_12DVP <- read_sav("TECH/P1_12DVP.sav") # deviant peers
P1_3IPC <- read_sav("TECH/P1_3IPC.sav") # playground with peers
P1_3PNOMN <- read_sav("TECH/P1_3PNOMN.sav")  
P1_3SAB <- read_sav("TECH/P1_3SAB.sav")  
P1_3SI <- read_sav("TECH/P1_3SI.sav")  
P1_3TOCA <- read_sav("TECH/P1_3TOCA.sav")  
P1_7TPRSK <- read_sav("TECH/P1_7TPRSK.sav")  
P1_10AD <- read_sav("TECH/P1_10AD.sav")  
# P1_10DEPC <- read_sav("TECH/P1_10DEPC.sav")  

P1_10PCR <- read_sav("TECH/P1_10PCR.sav") # parent-child relationship
P1_10PIN <- read_sav("TECH/P1_10PIN.sav") # parent involement
P1_10PIR <- read_sav("TECH/P1_10PIR.sav")  

P1_12DIS <- read_sav("TECH/P1_12DIS.sav") # parent discipline
P1_12POSR <- read_sav("TECH/P1_12POSR.sav") # parent positive reinforcement
P1_12SUP <- read_sav("TECH/P1_12SUP.sav") # parent supervision and monitor
P1_12IPTL <- read_sav("TECH/P1_12IPTL.sav") # parent-child interactions
P1_12SAS <- read_sav("TECH/P1_12SAS.sav") 
P1_10DSM <- read_sav("TECH/P1_10DSM.sav") # Oppositional Defiant Disorder (OD) & Conduct Disorder (CD)
P1_12CSU <- read_sav("TECH/P1_12CSU.sav") # child substance use
P1_10WD <- read_sav("TECH/P1_10WD.sav") # child withdrawn
P1_10REJ <- read_sav("TECH/P1_10REJ.sav")  
P1_10SSK <- read_sav("TECH/P1_10SSK.sav")  
P1_10WL <- read_sav("TECH/P1_10WL.sav") 
P1_12AA <- read_sav("TECH/P1_12AA.sav") 
P1_12CHW <- read_sav("TECH/P1_12CHW.sav") 
P1_12HAS <- read_sav("TECH/P1_12HAS.sav") 
P1_12LAB <- read_sav("TECH/P1_12LAB.sav") 
P1_12PD <- read_sav("TECH/P1_12PD.sav") 
P1_12PIH <- read_sav("TECH/P1_12PIH.sav") 




# Load ---------------------
#load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/replicate_lift_others2.RData")

# data_in2 <- data_in
# datuse2 <- datuse
# imp_miss2 <- imp_miss
# datimp2 <- datimp

# tt (Treatment) & R (Moderator) ----
trt_dat <- PARTICIPATION %>% select(FAMILY, GRADE, GROUP, SEX) 

# Y: Outcome ---------------------
# child substance use 
y_dat_csub <- P1_12CSU %>% # P1_15CSUB %>%
  mutate(
    Y_substance10 = CTSUA1, 
    Y_tobacco10 = C81AJ4,
    Y_drug10 = C81DRGA1, 
    Y_substance8 = CTSU81, 
    Y_tobacco8 = C818J4, 
    Y_drug8 = C81DRG81
  ) %>%
  select(FAMILY, starts_with("Y_"), starts_with("C_"))




# join the outcome data
y_dat <- y_dat_csub #%>% 
# full_join(y_dat_cbehave, by =  "FAMILY") %>%
#full_join(y_dat_schoolantisoc, by =  "FAMILY") 


# # Or:
# y_dat_2 <- datuse %>% select(
#   starts_with("Y_"), c(FAMILY)
# )
# colnames(y_dat_2)[colMeans(apply(y_dat_2, 2, is.na)) > 0.36]
# y_dat <- y_dat_2[, colMeans(apply(y_dat_2, 2, is.na)) <= 0.36]



# M: Intermediate variables --------------
m_dat_wmtprefer <- P1_7TPRSK %>% dplyr::select(
  FAMILY,  TB4TP41, TB4TP51,TB4TP61,TB4TP71
) %>%
  rename_with(~paste0("M_teaprefer", .), c(!FAMILY))

m_dat_schoolantisoc <- P1_12SAS %>% dplyr::select(
  FAMILY,T42AS42, T42AS52, T42AS62,T42AS82
) %>%
  rename_with(~paste0("M_sas", .), c(!FAMILY))



# add more deviant peer time points
m_dat_dvpeers <- P1_12DVP %>% 
  # mutate(
  #   C_deviantpeers_teacher = TDVPR12
  # ,C_dvpP59139 = P59139, C_dvpPA8113 = PA8113, C_dvpP54DVP1=P54DVP1
  # ) %>%
  mutate(
    M_deviantpeers5 = TDVPR52, # P59539, #
    M_deviantpeers6 = TDVPR62, #P59639, 
    M_deviantpeers8 = TDVPR82,  #P59839, #
    M_deviantpeers4 = TDVPR42, #P59439 #
    M_dvpPOCA3 = P54DVP3, M_dvpPOCA4=P54DVP4, M_dvpPOCA5=P54DVP5, M_dvpPOCA6=P54DVP6, M_dvpPOCA8=P54DVP8, M_dvpPOCA10=P54DVPA
  ) %>% dplyr::select(FAMILY, starts_with("M_"), starts_with("C_"))


m_dat_risk <- P1_12RISK %>% mutate(
  # C_risk = TV8RSK11,
  M_risk8 = TV8RSK81, 
  M_risk4= TV8RSK41, M_risk5= TV8RSK51, M_risk6= TV8RSK61
) %>% dplyr::select(FAMILY, starts_with("M_"), starts_with("C_")) 


# m_dat_aversive <- P1_15_CHILD_ANTISOCIAL %>% dplyr::select(
#   FAMILY, CAAVLR41, CAAVLR51, CAAVLR61, CAAVLR81
# ) %>%
#   rename_with(~paste0("M_aversive", .), c(!FAMILY))


m_dat_aggress <- P1_15_CHILD_ANTISOCIAL %>% dplyr::select(
  FAMILY, P594AG2R, P595AG2R, P596AG2R, P598AG2R,
  T424AG2R, T425AG2R, T426AG2R, T428AG2R
) %>%
  rename_with(~paste0("M_aggress", .), c(!FAMILY))


m_dat_delinq <- P1_15_CHILD_ANTISOCIAL %>% dplyr::select(
  FAMILY, P594DL2R, P595DL2R, P596DL2R, P598DL2R,
  T424DL2R, T425DL2R, T426DL2R, T428DL2R
) %>%
  rename_with(~paste0("M_delinq", .), c(!FAMILY))

# m_dat_cbehave <- P1_15CBCL_SCALES %>% # P1_15_CHILD_ANTISOCIAL %>% # 
#   # mutate( 
#   #   C_aggress = T421AG2R, C_withdrawn = T421WT2R, C_somatic = T421SO2R,
#   #   C_anxiousdepress = T421AN2R, C_attentionprob = T421AP2R, 
#   #   C_delinq = T421DL2R, C_socialproblem = T421SP2R, C_thoughtproblem = T421TP2R
#   #   # C_external = T421EX2R, C_internal = T421IN2R, 
#   # ) %>%
#   mutate( 
#     # total score
#     M_totalscore3 = T423TO2R, M_totalscore4 = T424TO2R, M_totalscore5 = T425TO2R, M_totalscore6 = T426TO2R, M_totalscore8 = T428TO2R, M_totalscore10 = T42ATO2R,
#     # # social problems
#     # M_socialproblem3 = T423SP2R, M_socialproblem4 = T424SP2R, M_socialproblem5 = T425SP2R, M_socialproblem6 = T426SP2R, M_socialproblem8 = T428SP2R, M_socialproblem10 = T42ASP2,
#     # aggressive behaviors
#     M_aggress3 = T423AG2R, M_aggress4 = T424AG2R, M_aggress5 = T425AG2R, M_aggress6 = T426AG2R, M_aggress8 = T428AG2R, M_aggress10 = T42AAG2R
#   ) %>% select(FAMILY, starts_with("M_"), starts_with("C_"))


m_dat_parchild <- P1_10PCR %>% dplyr::select(
  FAMILY,
  P85PC42, P85PC52, P85PC62, P85PC82
) %>%
  rename_with(~paste0("M_PCR", .), c(!FAMILY))


m_dat_pardiscipline <- P1_12DIS %>% mutate(
  M_dispappr = P85APD41,
  M_dispinappr.pint = P85IPD41, M_dispinconsis.pint = P85ICD41
) %>% dplyr::select(FAMILY, starts_with("M_"), starts_with("C_"))

m_dat_playground <- P1_3IPC %>%
  mutate(
    M_aversivephy3 = CAAPR32, M_neutralphy = CANPR32,M_posphy = CAPPR32
  ) %>% dplyr::select(FAMILY, starts_with("M_"))



m_dat_familyinteraction <- P1_12IPTL %>%
  mutate(
    # Aversive Engagement
    M_chaversivefamily3 = CAAVLR31,
    M_chaversivefamily4 = CAAVLR41, M_chaversivefamily5 = CAAVLR51, M_chaversivefamily6 = CAAVLR61, M_chaversivefamily8 = CAAVLR81, M_chaversivefamily10 = CAAVLRA1,
    # Positive Engagement
    M_chldposfamily3 = CAPOLR31,
    M_chldposfamily4 = CAPOLR41,M_chldposfamily5 = CAPOLR51,M_chldposfamily6 = CAPOLR61,
    M_chldposfamily8 = CAPOLR81, M_chldposfamily10 = CAPOLRA1
  ) %>%
  dplyr::select(FAMILY, starts_with("M_"))



m_dat <- list(
  m_dat_dvpeers, 
  m_dat_risk,
  # m_dat_wmtprefer, 
  # m_dat_aversive,
  # m_dat_schoolantisoc, # collinear with aggressive behaviors
  m_dat_aggress,
  m_dat_delinq,
  
  m_dat_playground, 
  m_dat_familyinteraction, 
  m_dat_parchild,
  m_dat_pardiscipline
  # m_dat_cbehave
) %>% 
  reduce(full_join, by="FAMILY")

m_dat %>% summarise(across(.cols = starts_with("M_"), .fns = sd, na.rm=T))



# not standardize m, maintain the original scale
# m_dat <- m_dat %>% mutate(
#   # across(!c(FAMILY), ~as.numeric(scale(., center=T, scale=T)) )
#   across(!c(FAMILY), ~((. - mean(., na.rm=T))/sd(., na.rm=T)) )
# )
m_dat %>% summarise(across(everything(.), n_distinct))
m_dat %>% summarise(across(!any_of(c("FAMILY")), ~sd(., na.rm=T)))
aa<-m_dat %>% select(!c(FAMILY)) %>% cor(., use = "pairwise.complete.obs")
data.frame(aa) %>% summarise(across(everything(), ~range(.[. <1], na.rm=T)))


# # Or:
# m_dat_2 <- datuse %>% select(
#   starts_with("M_"), starts_with("X_"), c(FAMILY)
# )
# colnames(m_dat_2)[colMeans(apply(m_dat_2, 2, is.na)) > 0.36]
# m_dat <- m_dat_2[, colMeans(apply(m_dat_2, 2, is.na)) <= 0.36]


# other C: baseline covariates ------------------------------

# c_dat_2 <- datuse %>% dplyr::select(
#   !starts_with("M_") &!starts_with("X_") & !starts_with("Y_") & !c(GROUP,GRADE,SEX)
# )
# sapply(c_dat_2, attr, "label")
# 
# c_dat_2$C_risk <- select(P1_12RISK, TV8RSK11)
# c_dat_2$C_aggress <- select(P1_15_CHILD_ANTISOCIAL, c(P591AG2R,  T421AG2R))
# c_dat_2$C_delinq <- select(P1_15_CHILD_ANTISOCIAL, c(P591DL2R,  T421DL2R))
# c_dat_2$C_internal <- select(P1_15CBCL_SCALES, c(T421IN2R))
# c_dat_2$C_reinforce.lab <- select(P1_12POSR, c(P3VPOS11))
# c_dat_2$C_reinforce.pint <- select(P1_12POSR, c(M851D6, F851D6))
# c_dat_2$C_deviantpeers_teacher <- select(P1_12DVP, c(TDVPR12))
# c_dat_2$X_deviantpeers <- select(P1_12DVP, c(TDVPR32))
# 
# c_dat_2 <- c_dat_2 %>% dplyr::select(!c(C_external, C_drughelp))
# c_dat_2 <- c_dat_2 %>% dplyr::select(!c(C81ICD11, M3V1C64, C81OD11A)) # missing 1st grade
# 
# c_dat_2 <- do.call(data.frame, c_dat_2)
# 
# c_dat_2 %>% summarise(across(.cols = !c(FAMILY), .fns = sd, na.rm=T))
# 
# 
# colnames(c_dat_2)[colMeans(apply(c_dat_2, 2, is.na)) > 0.25]
# c_dat <- c_dat_2[, colMeans(apply(c_dat_2, 2, is.na)) <= 0.25]



# parents substance use
c_dat_parents_substance <- P1_12PARENT_SUBSTANCE %>% select(
  FAMILY,
  FV9108, # FV9109C,FV9109D,
  MV9108, # MV9109C, MV9109D
) %>%
  rename_with(~paste0("C_parentdrink", .), c(!FAMILY))


c_dat_familystatus <- P1_12FAMILY_STATUS %>% select(
  FAMILY,
  Q11SIB1, P11SES1, INCOME1, FEDUC1, MEDUC1
) %>%
  rename_with(~paste0("C_famst", .), c(!FAMILY))


c_dat_risk <- P1_12RISK %>% select(
  FAMILY,
  TV8RSK11
) %>%
  rename_with(~paste0("C_risk", .), c(!FAMILY))


c_dat_deviantpeer <- P1_12DVP %>% select(
  FAMILY,
  TDVPR12,P54DVP1
) %>%
  rename_with(~paste0("C_deviantpeer", .), c(!FAMILY))


c_dat_schlantisoc <- P1_12SAS %>% select(
  FAMILY,
  T42AS12
) %>%
  rename_with(~paste0("C_schlantisoc", .), c(!FAMILY))


c_dat_cbcl <- P1_15CBCL_SCALES %>% select(
  FAMILY, T421DL2R, T421AG2R, 
  M591DL2R,F591DL2R,M591AG2R,F591AG2R, 
  T421IN2R, M591IN2R
  # M591TO2R,F591TO2R,
  # T421WT2R, T421SO2R,
  # T421AN2R, T421AP2R, 
  # T421SP2R, T421TP2R
) %>%
  rename_with(~paste0("C_cbcl", .), c(!FAMILY))




c_dat_schladjust <- P1_3SAB %>% select(
  FAMILY, TV8SA11, Q1V1PE, QV6SA11
) %>%
  rename_with(~paste0("C_schladjust", .), c(!FAMILY))



c_dat_parchild <- P1_10PCR %>% select(
  FAMILY,
  P85PC12, P86PC11
) %>%
  rename_with(~paste0("C_pcrelation", .), c(!FAMILY))


c_dat_pardiscipline <- P1_12DIS %>% select(
  FAMILY,P85APD11, P85ICD11, P86APD11 #, P86IPDC1
) %>% rename_with(~paste0("C_discipline", .), c(!FAMILY))


c_dat_parsupervision <- P1_12SUP %>% select(
  FAMILY,P86120, MV4MN12, P54MN12, FV4MN12,MV4MN12, TB4119  #C82112, 
) %>%
  rename_with(~paste0("C_supervision", .), c(!FAMILY))



c_dat_pinvol <- P1_10PIN %>% select(
  FAMILY, P85PI11, P27PI11
) %>%
  rename_with(~paste0("C_parentinvole", .), c(!FAMILY))


# c_dat_pinvolhw <- P1_12PIH %>% select(
#   FAMILY, P85CA7 # missing 
# ) %>%
#   rename_with(~paste0("C_pinvolhw", .), c(!FAMILY))


c_dat_attdef <- P1_10AD %>% select(
  FAMILY,  P59AD11, PCAD11,T42AD11
) %>%
  rename_with(~paste0("C_attdef", .), c(!FAMILY))

c_dat_dsm <- P1_10DSM %>% select(
  FAMILY, PDSO12A  
) %>%
  rename_with(~paste0("C_dsm", .), c(!FAMILY))


c_dat_homantisoc <- P1_12HAS %>% select(
  FAMILY, PASCD11, PASOD11  
) %>%
  rename_with(~paste0("C_homantisoc", .), c(!FAMILY))


c_dat_parirritable <- P1_10PIR %>% select(
  FAMILY, P86PIR11,P851G4  
) %>%
  rename_with(~paste0("C_pirritable", .), c(!FAMILY))


c_dat_rejected <- P1_10REJ %>% select(
  FAMILY, PCREJ11,TCREJ12, 
) %>%
  rename_with(~paste0("C_rejected", .), c(!FAMILY))


c_dat_socialskill <- P1_10SSK %>% select(
  FAMILY,  C82SS11
) %>%
  rename_with(~paste0("C_socialskill", .), c(!FAMILY))


c_dat_wmscale <- P1_7TPRSK %>% select(
  FAMILY,  TB4PP12,TB4TP12
) %>%
  rename_with(~paste0("C_wmscale", .), c(!FAMILY))


c_dat_withdrawn <- P1_10WD %>% select(
  FAMILY,  C82103,PCWD11
) %>%
  rename_with(~paste0("C_withdrawn", .), c(!FAMILY))


c_dat_weliked <- P1_10WL %>% select(
  FAMILY,  TB4102,PA8WL11
) %>%
  rename_with(~paste0("C_weliked", .), c(!FAMILY))



c_dat_academic <- P1_12AA %>% select(
  FAMILY, T42AA11, P59AA11
) %>%
  rename_with(~paste0("C_academic", .), c(!FAMILY))


c_dat_hwengage <- P1_12CHW %>% select(
  FAMILY,  P85HE11,T421A9
) %>%
  rename_with(~paste0("C_hwengage", .), c(!FAMILY))


c_dat_iptl <- P1_12IPTL %>% select(
  FAMILY, CAAPLR11,CAAVLR11,
  CAPOLR11,CAPPLR11,
  # CACOLR11,CAADLR11,
  # CADRLR11,CADSLR11,
  # CANELR11,CANPLR11,
  MAAPLR11,MAAVLR11,
  MAPOLR11,MAPPLR11
) %>%
  rename_with(~paste0("C_iptl", .), c(!FAMILY))

c_dat_ipc <- P1_3IPC %>% select(
  FAMILY, CAAPR12,CAAVR12,CAPOR12,CAPPR12
  # ,CACOR12,CADRR12,CADSR12,
  ,CANER12,CANPR12
) %>%
  rename_with(~paste0("C_ipc", .), c(!FAMILY))


c_dat_tocasocial <- P1_3TOCA %>% select(
  FAMILY, TV8SC11,TV8AR11,TV8MA11, TV8CA11,TV8CN11,TV8HY11
) %>%
  rename_with(~paste0("C_tocasocial", .), c(!FAMILY))


c_dat_pardep <- P1_12PD %>% select(
  FAMILY, M57PD11,F57PD11
) %>%
  rename_with(~paste0("C_pd", .), c(!FAMILY))


c_dat_posreinforce <- P1_12POSR %>% select(
  FAMILY,  P86123, P3VPOS11
) %>%
  rename_with(~paste0("C_posreinf", .), c(!FAMILY))






# from scales

attr(P1SCLE_p1$TB4TO11, "label")

# c_dat_scale <- P1SCLE_p1 %>% select(
#   FAMILY,# Walker McConnell - Walker-McConnell Scale of Social Competence and School Adjustment,
#   TB4TO11 # total
#   #TB4PP11, # peer preferred of the Walker McConnell
#   # TB4TP11 # teacher preferred
#   # TB4SA11 # school adjustment
# ) %>%
#   rename_with(~paste0("C_ssk", .), c(!FAMILY))





# join the baseline covariates
clist <- ls() %>% grep("c_dat", ., value = T)
clist <- clist[!clist%in%c("c_dat","c_dat_1")]
# > clist
# [1] "c_dat_academic"         
# [2] "c_dat_attdef"           
# [3] "c_dat_cbcl"             
# [4] "c_dat_deviantpeer"      
# [5] "c_dat_dsm"              
# [6] "c_dat_familystatus"     
# [7] "c_dat_homantisoc"       
# [8] "c_dat_hwengage"         
# [9] "c_dat_ipc"              
# [10] "c_dat_iptl"             
# [11] "c_dat_parchild"         
# [12] "c_dat_pardep"           
# [13] "c_dat_pardiscipline"    
# [14] "c_dat_parents_substance"
# [15] "c_dat_parirritable"     
# [16] "c_dat_parsupervision"   
# [17] "c_dat_pinvol"           
# [18] "c_dat_posreinforce"     
# [19] "c_dat_rejected"         
# [20] "c_dat_risk"             
# [21] "c_dat_schladjust"       
# [22] "c_dat_schlantisoc"      
# [23] "c_dat_socialskill"      
# [24] "c_dat_tocasocial"       
# [25] "c_dat_weliked"          
# [26] "c_dat_withdrawn"        
# [27] "c_dat_wmscale"
cdatlist <- mget(clist)
c_dat <- cdatlist %>% reduce(full_join, by="FAMILY")


c_dat_1 <- c_dat
c_dat_1 %>% summarise(across(.cols = !c(FAMILY), .fns = sd, na.rm=T))


colnames(c_dat_1)[colMeans(apply(c_dat_1, 2, is.na)) > 0.25]
c_dat <- c_dat_1[, colMeans(apply(c_dat_1, 2, is.na)) <= 0.25]



# standardize to facilitate imputation
c_dat <- c_dat %>% mutate(
  # across(!c(FAMILY), ~as.numeric(scale(., center=T, scale=T)) )
  across(!c(FAMILY), ~((. - mean(., na.rm=T))/sd(., na.rm=T)) )
)
# check distinct values & collinearity
c_dat %>% summarise(across(!any_of(c("FAMILY")), ~n_distinct(., na.rm=T)))
aa <-c_dat %>% select(!c(FAMILY)) %>% cor(., use = "pairwise.complete.obs")
data.frame(aa) %>% summarise(across(everything(), ~range(.[. <1], na.rm=T)))

# reduce collinear variables
c_dat <- c_dat %>% select(
  !contains("TV8CN11") & !contains("TV8CA11") & # collinear with C_schladjustTV8SA
    !contains("T42AS12") & # col with T421AG2R
    !contains("P59AD11") # col with PCAD11
)

# Joining Data --------------------


df_use <- list(
  trt_dat, m_dat, y_dat, c_dat
) %>% reduce(full_join, by="FAMILY")
colnames(df_use)


glm(I((Y_tobacco10>=1) ) ~ factor(GROUP) *  I(SEX==2) , data = df_use[df_use$GRADE==2, ], family="binomial") %>% summary()

# ___________ -----

# (NOT USED)  --------------------
# 
# load("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/lift_others2.RData")
# aa<-sapply(datuse, attr, "label")
# label_childantisoc <- sapply(P1_15_CHILD_ANTISOCIAL, attr, "label")
# label_P1SCLE %>% grep("rein", ., ignore.case=T,value = T)
# label_childantisoc %>% grep("deviant", ., ignore.case=T,value = T)
#   
# df_y <- P1_12CSU %>% select(FAMILY, C81AJ4) %>% rename(Y_tobacco10=C81AJ4)
# df_m <- P1_12DVP %>% select(FAMILY, TDVPR82, TDVPR32, TDVPR42, TDVPR52, TDVPR62) %>% rename(M_deviantpeers8=TDVPR82) %>%
#   full_join(., select(P1_12RISK, c(FAMILY, TV8RSK81, TV8RSK31, TV8RSK41,TV8RSK51,TV8RSK61)), by = "FAMILY") %>% rename(M_risk8=TV8RSK81)
# 
# 
# df_c1 <- P1_15_CHILD_ANTISOCIAL %>% select(FAMILY, T421DL2R, T421AG2R, T421EX2R) # delinquent, aggressive, external
# df_c2 <- P1_15CBCL_SCALES %>% select(FAMILY, T421IN2) # internal
# df_c3 <- P1SCLE_p1 %>% select(FAMILY, TV8RSK11, TDVPR12, P85POS11) # risk, deviantpeers, positive-refinforce
# 
# 
# df_c4 <- datuse %>% select(
#   !starts_with("Y_risk") & !starts_with("Y_tobacco10") & !starts_with("M_deviantpeers") & !starts_with("X_") &
#   !c(C_deviantpeers_teacher, C_risk, C_external, C_delinq, C_internal, C_aggress, C_reinforce, C_reinforce.lab, C_reinforce.pint, 
#      C_dispappr, M3V1C64,C81ICD11, C81OD11A, C82112, CV5SS11,
#      R43AG11,R43AX11,R43LK11,R43SD11,R43WD11,R43SS11,R43RJ11,CV5AD11,PASCD11) 
# )
# 
# df_use <- df_c4 %>% 
#   full_join(., df_c3, by = "FAMILY") %>% 
#   full_join(., df_c2, by = "FAMILY") %>%
#   full_join(., df_c1, by = "FAMILY") %>%
#   full_join(., df_m, by = "FAMILY") %>%
#   full_join(., df_y, by = "FAMILY")
#   
# 
# 
# 
# # join the data on baseline covariates C, moderator R, treatment tt, posttreatment variables X and M, outcome Y
# 
# datuse <- trt_dat %>% 
#   full_join(., y_dat, by = "FAMILY") %>% 
#   full_join(., m_dat, by = "FAMILY") %>%
#   #full_join(., x_dat, by = "FAMILY") %>%
#   full_join(., c_dat, by = "FAMILY")
# 
# 
# colMeans(apply(datuse, 2, is.na))
# colMeans(apply(datuse[datuse$GRADE==2, ], 2, is.na))



# ___________ -----

# Impute missing data ------------------

datmice <- df_use %>% dplyr::select(
  FAMILY, GRADE,
  starts_with("C_"), colnames(c_dat),
  GROUP, SEX,
  starts_with("X_"), starts_with("M_"), starts_with("Y_")
)


datmice <- datmice %>% dplyr::select( !FAMILY )
datmice <- do.call(data.frame, datmice)
datmice <- datmice %>% mutate_all(~as.numeric(.))
colnames(datmice)

datmice %>% filter(GRADE==2) %>%
  summarise_all(., ~mean(is.na(.)) ) 

datmice %>% summarise_all(., ~sd(., na.rm=T) ) 
# variables with relatively low missing rates
# acceptmiss <- (colMeans(apply(datmice[datmice$GRADE==2, ], 2, is.na)) < 0.36) * (colMeans(apply(datmice, 2, is.na)) < 0.36)

# datmi <- datmice[, acceptmiss==1] 
datmi <- datmice
colnames(datmi)[colMeans(apply(datmi, 2, is.na)) > 0.36]
datmi %>% group_by(GRADE) %>%
  summarise_all(., ~mean(is.na(.)) ) 

# PARTICIPATION$CINT10
datmi <- datmi %>%
  mutate(
    SEX = as.factor(SEX), GROUP = as.factor(GROUP)
  ) #%>% filter(GRADE==2) %>% select(!GRADE)





# Or:
# datmi <- imp_miss2$data
summary(datmi)


mipred <- quickpred(datmi, mincor = 0.1)
mimthd <- make.method(datmi)


eeeeee



# imp_miss134 <- mice(datmi, m=1, method = mimthd, predictorMatrix = mipred, maxit = 50, seed = 12345)


# # sensitivity analysis for "maxit", check 100
imp_maxit100 <- mice(datmi, m=1, method = mimthd, predictorMatrix = mipred, maxit = 100, seed = 12345)

# removed the duplicated m_aversivefamily from the P1_12IPTL
# imp_maxit120 <- mice(datmi, m=1, method = mimthd, predictorMatrix = mipred, maxit = 120, seed = 12345)


# # missing indicators of covariates
# missind <- imp_miss[["where"]]
# C_miss <- data.frame(missind) %>% select(
#   starts_with("C_")
# ) %>% mutate_all(~as.numeric(.)) 
# Cmiss <- C_miss[, colMeans(C_miss) > 0.2]

# datimp1 <- complete(imp_miss134)
datimp1 <- complete(imp_maxit100)
# datimp1 <- complete(imp_maxit120)

colnames(datimp1)[summarise_all(datimp1, ~mean(is.na(.)) )>0 ]
datimp1 <- datimp1[, summarise_all(datimp1, ~mean(is.na(.)) )==0 ]

# Grade 5 only ----
# datimp <- datimp %>% filter(GRADE == 2) %>% dplyr::select( !GRADE )
# datimp <- datimp[complete.cases(datimp), ]

datimp1 <- datimp1 %>% filter(GRADE == 2) %>% dplyr::select( !GRADE ) 


sum(!complete.cases(datimp1))

datimp1 %>% group_by(GROUP, SEX) %>%
  summarise( mean(Y_tobacco10>=1))

datimp1 %>% group_by(SEX) %>%
  summarise_at(vars(starts_with("M_a"), Y_tobacco10 ),
               list(~mean(., na.rm = T)))

glm(I(Y_tobacco10 >= 1) ~ factor(GROUP) *  I(SEX==2) , data = datimp1, family="binomial") %>% summary()

glm(I(Y_tobacco10 >= 1) ~ I(GROUP==1) *  I(SEX==2) , data = datimp1, family="gaussian") %>% summary()

glm(I(Y_tobacco10 >= 1) ~ I(GROUP==1)*I(M_aggressT424AG2R-M_aggressT426AG2R)  +I(GROUP==1) *  I(SEX==2), data = datimp1, family="gaussian") %>% summary()





# C selecting --------------

Cselected <- datimp1 %>% dplyr::select(
  C_parentdrinkMV9108,
  C_famstP11SES1, C_famstINCOME1, 
  
  C_riskTV8RSK11,
  C_cbclT421DL2R,
  C_cbclT421AG2R, C_cbclM591AG2R,
  C_deviantpeerTDVPR12, C_deviantpeerP54DVP1,
  C_cbclT421IN2R,
  
  C_ipcCAAPR12, C_ipcCAAVR12, C_iptlCAAPLR11, C_iptlCAAVLR11,
  
  C_schladjustTV8SA11,
  C_rejectedTCREJ12,
  C_welikedTB4102, C_wmscaleTB4PP12, C_wmscaleTB4TP12,
  C_socialskillC82SS11,
  
  C_homantisocPASOD11,
  C_pcrelationP85PC12, C_pdM57PD11, C_disciplineP85APD11, 
  C_parentinvoleP85PI11, 
  C_academicT42AA11
)



# Cselected <- datimp %>% dplyr::select( !c(GROUP, SEX) & !starts_with("Y_") & !starts_with("M_") & !starts_with("X_") & c(
#   C_parentsdrink,
#   C_risk, C_delinq, C_schlantisoc,
#   C_internal,C_aggress,
#   C_deviantpeers_teacher,
#   P11SES1,  INCOME1,FEDUC1, MEDUC1, Q11SIB1, # C_parage,
#   C_familystatus.x_2, C_familystatus.x_3,
#   C_dispappr.pint, C_dispinconsis.pint,  #C_dispappr,
#   C_reinforce, C_reinforce.lab,
#   C_overtcovertaggr,
#   # C_parchldrel,
#   P85PC12, # mean mom-dad/child relationship [PINT version 2]
#   TAACH11, # academic achievement TOCA/TCBC
#   CC9AA11, # wiscr
#   MV4MN12, # mom monitor
#   # R43AG11,R43AX11,R43LK11,R43SD11,R43WD11, R43RJ11, # PEI scale
#   TV8AR11, # authority rejection, toca
#   TV8MA11, # maturity
#   TV8CA11, # cognitive achievement
#   # TV8CN11, # concentration (colinear with TV8CA11)
#   TV8HY11, # hyperactivity
#   TB4TO13, # walker mcconnell percentile total
#   PSCIN13, # school invovlment
#   P31PIH1, # PAR involve with homework
#   TCWD11, # teacher measure of child withdrawn
#   TWLLK11, # teacher measure of child well-liked
#   TV8SS11, # social skills
#   CV5SS11, # social skills soimp
#   CV5103, # playground impression, rejected child
#   QV6103, # classroom observation, classroom antisocial
#   T42AD11, # teacher measure of attention deficient
#   CV5AD11, # soimp att def
#   C82CAS11, # ciimp global child antisocial
#   PASCD11, # parent measure of conduct disorder
#   M57PD11, # mom CESDD depression
#   # the playground observations and lab task observations could be partly used in constructing the selected covariates, based on the TECH documentation
#   C_posphysical,C_averphysical,C_positiveEngagement,C_aversiveEngagement,
#   C_aversiveeng,C_positveengage, # C_momaversiveeng,
#   # CACOR12,CADRR12,CADSR12, CANPR12, #
#   CACOR11,CADRR11,CADSR11,
#   CANPR11,CANER12,
#   PAAVLR11
# )  )

colnames(c_dat)[!colnames(c_dat)%in%colnames(Cselected)] # collinear nonselection


# Data in ----------------

# !updated to datimp1
data_in <- data.frame(
  Y = as.numeric(datimp1$Y_tobacco10 >= 1),
  # Mdat = dplyr::select(datimp1, c(
  #   M_risk5, M_deviantpeers5 #, M_risk6, M_deviantpeers6
  #   #M_delinqP595DL2R, M_aggressP595AG2R, M_delinqT425DL2R,M_aggressT425AG2R
  #   )),
  Mdat = data.frame(
    #M_riskslope4to5 = datimp1$M_risk5-datimp1$M_risk4,
    M_risk8 = datimp1$M_risk8 , 
    M_deviantpeers8 = datimp1$M_deviantpeers8 #,
    #M_dvpeerslope4to5=datimp1$M_deviantpeers5-datimp1$M_deviantpeers4
  ),
  tt = as.numeric(datimp1$GROUP == 1), # 1 for intervention condition
  R = as.numeric(datimp1$SEX == 1), # 1 for male , 2 for female
  Cdat = dplyr::select(datimp1, !c(GROUP, SEX) & !starts_with("Y_") & !starts_with("M_") & !starts_with("X_") & all_of(colnames(Cselected))  
  )
  # ,Cdat_miss = Cmiss
)



# --------------------------------------
data_in$id <- 1:nrow(data_in)
data_in$obs_weights_Y <- rep(1, nrow(data_in))
Cnames <- colnames(dplyr::select(
  data_in, c(starts_with("Cdat") )
))
Mnames <- colnames(dplyr::select(data_in, starts_with("Mdat"), starts_with("M_")))


# (no need here for loaded data) add interactions
data_in <- data_in %>% mutate(
  ttR = data_in$tt * data_in$R,
  ttM = data_in$tt*data_in[, Mnames],
  RM = data_in$R*data_in[, Mnames],
  ttRM = ttR*data_in[, Mnames]
)

data_in <- do.call(data.frame, data_in)
ttMnames <- paste0("ttM.", Mnames)
RMnames <- paste0("RM.", Mnames)
ttRMnames <- paste0("ttRM.", Mnames)
colnames(data_in)[(ncol(data_in)+1-3*length(Mnames)): ncol(data_in)] <- c(ttMnames, RMnames, ttRMnames)

# try
# train_data <- valid_data <- data_in





# Estimation -----
source("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/simulation/estimators.R")
source("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/simulation/fitmodels.R")

# r1 for male ------------
# data_in$R <- as.numeric(datimp$SEX == 1) # 1 for male , 2 for female
# !updated to datimp1
data_in$R <- as.numeric(datimp1$SEX == 1)

table(data_in$R)
trt_dat %>% filter(GRADE==2) %>% group_by(SEX) %>% summarise( n() )
# update the interactions
data_in[, "ttR"] <- data_in$tt * data_in$R
data_in[, RMnames] <- data_in$R*data_in[, Mnames]
data_in[, ttRMnames] <- data_in$ttR * data_in[, Mnames]



lin.crossfit_Rmale <- crossfit.onestep(
  cv_folds = 0L,
  data_in = data_in,
  Cnames = Cnames,
  Mnames = Mnames,
  # Xnames = Xnames,
  Fit = "linear",
  Yfamily = "binomial"
)
lin.crossfit_Rmale$estimates
lin.crossfit_Rmale$z_intervals



mlr.crossfit_Rmale <- crossfit.onestep(
  cv_folds = 5,
  data_in = data_in, 
  Cnames = Cnames, 
  Mnames = Mnames,
  Fit = "mlr", 
  Yfamily = "binomial"  
)
# SL_user <- c("SL.glm", "SL.xgboost", "SL.ranger") 
# SL_user.vfit <- c("SL.glm.scaledY", "SL.xgboost", "SL.ranger") 
mlr.crossfit_Rmale$estimates
mlr.crossfit_Rmale$z_intervals



# tmle_Rmale_Mrisk8dvp8
mlr.tmle_Rmale <- tmle.medMO(
  crossfit_onestep = mlr.crossfit_Rmale,
  data_in = data_in,
  Cnames = Cnames,
  # Xnames = Xnames,
  Mnames = Mnames
)

mlr.tmle_Rmale$tml_estimates
mlr.tmle_Rmale$tml_intervals




# Bootstrap ---------------------

one.boot <- function(b) {
  boo <- sample(1:nrow(data_in), nrow(data_in), replace = TRUE)
  boo_data_in <- data_in[boo, ]
  
  # one-step
  boo_crossfit_Rmale <- crossfit.onestep(
    cv_folds = 5,
    # bootstrap data
    data_in = boo_data_in, 
    Cnames = Cnames, 
    Mnames = Mnames,
    Fit = "mlr", 
    Yfamily = "binomial"  
  )
  
  boo_crossfit <- boo_crossfit_Rmale$estimates
  
  # tml
  boo_tmle_Rmale <- tmle.medMO(
    crossfit_onestep = boo_crossfit_Rmale,
    data_in = boo_data_in,
    Cnames = Cnames,
    Mnames = Mnames
  )
  
  boo_tmle <- boo_tmle_Rmale$tml_estimates
  
  boo_out <- data.frame(
    effect = names(boo_crossfit),
    boo_crossfit=boo_crossfit, 
    boo_tmle=boo_tmle)
  
  return(boo_out)
}


set.seed(12345)
boot_res <- lapply(1:1000, one.boot)






save.image("~/Library/CloudStorage/Box-Box/Labs/Mediated/R_mediated/lift1229.RData")

# (haven't yet) r1 for female ------------------
data_in$R <- as.numeric(datimp$SEX == 2) # 1 for male , 2 for female
table(data_in$R)
trt_dat %>% filter(GRADE==2) %>% group_by(SEX) %>% summarise( n() )

# update the interactions
data_in[, "ttR"] <- data_in$tt * data_in$R
data_in[, RMnames] <- data_in$R*data_in[, Mnames]
data_in[, ttRMnames] <- data_in$ttR * data_in[, Mnames]


colnames(data_in)


lin.crossfit_Rfe <- crossfit.onestep(
  cv_folds = 5L,
  data_in = data_in, 
  Cnames = Cnames, 
  Mnames = Mnames,
  # Xnames = Xnames,
  Fit = "linear",
  Yfamily = "binomial"  
)
lin.crossfit_Rfe$estimates
lin.crossfit_Rfe$z_intervals
# crossfit_Rfe_Mdvp8



mlr.crossfit_Rfe <- crossfit.onestep(
  cv_folds = 5L,
  data_in = data_in, 
  Cnames = Cnames, 
  Mnames = Mnames,
  # Xnames = Xnames,
  Fit = "mlr",
  Yfamily = "binomial"  
)
mlr.crossfit_Rfe$estimates
mlr.crossfit_Rfe$z_intervals



crossfit_Rfe_Mrisk8dvp8 <- crossfit.onestep(
  cv_folds = 5L,
  data_in = data_in, 
  Cnames = Cnames, 
  Mnames = Mnames,
  # Xnames = Xnames,
  Yfamily = "binomial" # "gaussian" #  "zinflasso" # "zeroinfl" # "zeroinfl", ##c(familyzero = "binomial", familypositive = "gaussian") 
)

crossfit_out <- crossfit_Rfe_Mrisk8dvp8 # crossfit_Rfe_Mdvp8 
crossfit_out$estimates
crossfit_out$z_intervals

tmle_Rfe_Mdvp8

tmle_Rfe_Mrisk8dvp8 <- tmle.medMO(
  crossfit_onestep = crossfit_Rfe_Mrisk8dvp8,
  data_in = data_in,
  Cnames = Cnames,
  # Xnames = Xnames,
  Mnames = Mnames
)

tmle_out <- tmle_Rfe_Mrisk8dvp8 # tmle_Rfe_Mdvp8 # tmle_Rfe_Mdvp456
tmle_out$tml_estimates
tmle_out$tml_intervals



