setwd('C:/Users/pakma/Projects/1. ADHD GRF/For_share/Data (ABCD 4.0)')


############################## @ Libraries ##############################
library(dplyr)
library(grf)
library(stringr)
library(caret)
library(tibble)
library(fastDummies)
library(tidyr)
library(impute)
library(gtools)
library(reshape2)
library(ggplot2)
library(ggstatsplot)
library(car)
library(readxl)


############################## @ Treatment, outcome, covariate data ##############################
pdem02_vars <- c('subjectkey','interview_age','sex','demo_prnt_marital_v2',
                 'demo_prnt_ed_v2','demo_comb_income_v2')
pdem02 <- read.table('pdem02.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  select(all_of(pdem02_vars))

acspsw03_vars <- c('subjectkey','race_ethnicity')
acspsw03 <- read.table('acspsw03.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  select(all_of(acspsw03_vars))

dhx01_vars <- c('subjectkey','devhx_10c3_p','devhx_10a3_p','devhx_10b3_p','devhx_10d3_p','devhx_10e3_p','devhx_10f3_p',
                'devhx_10g3_p','devhx_10h3_p','devhx_10i3_p','devhx_10j3_p','devhx_10k3_p','devhx_10l3_p','devhx_10m3_p',
                'devhx_8_alcohol','devhx_9_alcohol','devhx_8_tobacco','devhx_9_tobacco','devhx_8_marijuana','devhx_9_marijuana',
                'devhx_8_coc_crack','devhx_9_coc_crack','devhx_12a_p','birth_weight_lbs','birth_weight_oz','devhx_14a3_p',
                'devhx_14b3_p','devhx_14c3_p','devhx_14d3_p','devhx_14e3_p','devhx_14f3_p','devhx_14g3_p','devhx_14h3_p',
                'devhx_15','devhx_16_p','devhx_17_p','devhx_18_p','devhx_3_p','devhx_4_p','devhx_6_p')
dhx01 <- read_excel('dhx01.xlsx') %>%
  filter(visit == 'baseline_year_1_arm_1') %>% 
  select(all_of(dhx01_vars))

abcd_mx01_vars <- c('subjectkey','medhx_2b','medhx_6i')
abcd_mx01 <- read.table('abcd_mx01.txt', header = T, sep = "\t", dec = ".") %>%
  filter(eventname == 'baseline_year_1_arm_1') %>% 
  select(all_of(abcd_mx01_vars))

abcd_lt01_vars <- c('subjectkey', 'site_id_l')
abcd_lt01 <- read.table('abcd_lt01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_lt01_vars))

abcd_ssphp01_vars <- c('subjectkey', 'sds_p_ss_total')
abcd_ssphp01 <- read.table('abcd_ssphp01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_ssphp01_vars))

abcd_sscey01_vars <- c('subjectkey', 'crpbi_y_ss_parent','fes_y_ss_fc')
abcd_sscey01 <- read.table('abcd_sscey01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_sscey01_vars))

abcd_fhxssp01_vars <- c('subjectkey', 'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p')
abcd_fhxssp01 <- read.table('abcd_fhxssp01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_fhxssp01_vars))

abcd_asrs01_vars <- c('subjectkey','asr_scr_adhd_t','asr_scr_anxdisord_t','asr_scr_antisocial_t')
abcd_asrs01 <- read.table('abcd_asrs01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_asrs01_vars))

abcd_rhds01_vars <- c('subjectkey','reshist_addr1_pm252016aa','reshist_addr1_adi_wsum')
abcd_rhds01 <- read.table('abcd_rhds01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_rhds01_vars))

abcd_ant01_vars <- c('subjectkey','anthro_1_height_in','anthroweight1lb')
abcd_ant01 <- read.table('abcd_ant01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_ant01_vars))

abcd_cbcls01_vars <- c('subjectkey','cbcl_scr_dsm5_adhd_t','cbcl_scr_syn_anxdep_t','cbcl_scr_syn_withdep_t','cbcl_scr_syn_somatic_t',
                       'cbcl_scr_syn_social_t','cbcl_scr_syn_thought_t','cbcl_scr_syn_rulebreak_t','cbcl_scr_syn_aggressive_t')
abcd_cbcls01 <- read.table('abcd_cbcls01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_cbcls01_vars))

abcd_cbcls01_long_vars <- c('subjectkey','cbcl_scr_dsm5_adhd_t','eventname')
abcd_cbcls01_long <- read.table('abcd_cbcls01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname != 'baseline_year_1_arm_1')  %>%
  select(all_of(abcd_cbcls01_long_vars))

PRS <- read.csv('ABCD_EUR+Multibased_PRScsx_PC1-10resid_raw.csv', header = T, sep = ",", dec = ".") %>%
  select(-all_of('X'))

#PRS_holdout <- read.csv('ABCD_val+test_PRScsx_PC1-10resid_rawALL.csv', header = T, sep = ",", dec = ".") %>%
#  select(-all_of('X')) %>%
#  filter(set == 'val') %>%
#  select(!all_of(c('X.1','set')))



############################## @ Preprocessing (treatment, outcome, covariate Data) ##############################
#pdem02
not_numeric <- c('subjectkey','sex','demo_prnt_marital_v2')
numeric_cols <- names(pdem02)[!(names(pdem02) %in% not_numeric)]
for(i in 2:ncol(pdem02)){
  if(names(pdem02)[i] %in% numeric_cols){
    pdem02[,i] <- as.numeric(pdem02[,i])
  }
}

pdem02$sex[pdem02$sex == 'M'] <- 0  
pdem02$sex[pdem02$sex == 'F'] <- 1 
pdem02$sex <- as.numeric(pdem02$sex)
pdem02$demo_prnt_marital_v2[pdem02$demo_prnt_marital_v2 == '777'] <- NA
pdem02$demo_prnt_marital_v2[pdem02$demo_prnt_marital_v2 == ''] <- NA
pdem02$demo_prnt_ed_v2[pdem02$demo_prnt_ed_v2 == 777] <- NA
pdem02$demo_comb_income_v2[pdem02$demo_comb_income_v2 == 777] <- NA
pdem02$demo_comb_income_v2[pdem02$demo_comb_income_v2 == 999] <- NA

pdem02_1 <- pdem02[, names(pdem02) %in% c('subjectkey','interview_age','sex')]   #NA가 없는 것 
pdem02_2 <- pdem02[, names(pdem02) %in% c('subjectkey','demo_prnt_marital_v2')]  #NA를 제거해서 dummy variable을 만들 것 
pdem02_3 <- pdem02[, names(pdem02) %in% c('subjectkey','demo_prnt_ed_v2','demo_comb_income_v2')]  #NA가 있지만 imputation 할 것 

pdem02_2 <- pdem02_2 %>% drop_na()

cat_vars <- c('demo_prnt_marital_v2')
pdem02_2 <- dummy_cols(pdem02_2, select_columns = cat_vars) %>%
  subset(select=-c(demo_prnt_marital_v2))



#acspsw03
acspsw03$race_ethnicity[acspsw03$race_ethnicity == ''] <- NA
acspsw03 <- acspsw03 %>% drop_na()
cat_vars <- c('race_ethnicity')
acspsw03 <- dummy_cols(acspsw03, select_columns = cat_vars) %>%
  subset(select=-c(race_ethnicity))



#dhx01
dhx01 <- as.data.frame(dhx01)
not_numeric <- c('subjectkey')
numeric_cols <- names(dhx01)[!(names(dhx01) %in% not_numeric)]
for(i in 2:ncol(dhx01)){
  if(names(dhx01)[i] %in% numeric_cols){
    dhx01[,i] <- as.numeric(dhx01[,i])
  }
}

dhx01$devhx_10c3_p[dhx01$devhx_10c3_p == 999] <- NA
names(dhx01)[names(dhx01) == 'devhx_10c3_p'] <- 'preeclampsia'

ind1 <- which(names(dhx01) == 'devhx_10a3_p')
ind2 <- which(names(dhx01) == 'devhx_10m3_p')
for(i in ind1:ind2){
  dhx01[,i][dhx01[,i] == 999] <- NA
}
for(i in 1:nrow(dhx01)){
  dhx01$preg_comp[i] <- sum(dhx01[i, ind1:ind2])
}
dhx01 <- dhx01[, -c(ind1:ind2)]

dhx01$devhx_8_alcohol[dhx01$devhx_8_alcohol == 999] <- NA
dhx01$devhx_9_alcohol[dhx01$devhx_9_alcohol == 999] <- NA
for(i in 1:nrow(dhx01)){
  dhx01$alcohol[i] <- sum(dhx01$devhx_8_alcohol[i], dhx01$devhx_9_alcohol[i])
}
dhx01$alcohol[dhx01$alcohol == 2] <- 1

dhx01$devhx_8_tobacco[dhx01$devhx_8_tobacco == 999] <- NA
dhx01$devhx_9_tobacco[dhx01$devhx_9_tobacco == 999] <- NA
for(i in 1:nrow(dhx01)){
  dhx01$tobacco[i] <- sum(dhx01$devhx_8_tobacco[i], dhx01$devhx_9_tobacco[i])
}
dhx01$tobacco[dhx01$tobacco == 2] <- 1

dhx01$devhx_8_marijuana[dhx01$devhx_8_marijuana == 999] <- NA
dhx01$devhx_9_marijuana[dhx01$devhx_9_marijuana == 999] <- NA
for(i in 1:nrow(dhx01)){
  dhx01$marijuana[i] <- sum(dhx01$devhx_8_marijuana[i], dhx01$devhx_9_marijuana[i])
}
dhx01$marijuana[dhx01$marijuana == 2] <- 1

dhx01$devhx_8_coc_crack[dhx01$devhx_8_coc_crack == 999] <- NA
dhx01$devhx_9_coc_crack[dhx01$devhx_9_coc_crack == 999] <- NA
for(i in 1:nrow(dhx01)){
  dhx01$coc_crack[i] <- sum(dhx01$devhx_8_coc_crack[i], dhx01$devhx_9_coc_crack[i])
}
dhx01$coc_crack[dhx01$coc_crack == 2] <- 1

ind1 <- which(names(dhx01) == 'devhx_8_alcohol')
ind2 <- which(names(dhx01) == 'devhx_9_coc_crack')
dhx01 <- dhx01[, -c(ind1:ind2)]

dhx01$devhx_12a_p[dhx01$devhx_12a_p == 999] <- NA
names(dhx01)[names(dhx01) == 'devhx_12a_p'] <- 'preterm_birth'

ind <- which(names(dhx01) == 'birth_weight_oz')
dhx01 <- dhx01[, -ind]

ind1 <- which(names(dhx01) == 'devhx_14a3_p')
ind2 <- which(names(dhx01) == 'devhx_17_p')
for(i in ind1:ind2){
  dhx01[,i][which(dhx01[,i] == 999)] <- NA
}
dhx01$devhx_16_p[dhx01$devhx_16_p == 9990] <- NA
dhx01$devhx_15[dhx01$devhx_15 > 0] <- 1
dhx01$devhx_16_p[dhx01$devhx_16_p > 0] <- 1
dhx01$devhx_17_p[dhx01$devhx_17_p > 0] <- 1
ind3 <- which(names(dhx01) == 'devhx_16_p')
for(i in 1:nrow(dhx01)){
  dhx01$neo_illness[i] <- sum(dhx01[i,ind1:ind3])
}
dhx01 <- dhx01[, -c(ind1:ind2)]

names(dhx01)[names(dhx01) == 'devhx_18_p'] <- 'breastfeeding'
names(dhx01)[names(dhx01) == 'devhx_3_p'] <- 'maternal_age'
names(dhx01)[names(dhx01) == 'devhx_4_p'] <- 'paternal_age'
dhx01$paternal_age[dhx01$paternal_age %in% c(332, 389)] <- NA

dhx01$devhx_6_p[dhx01$devhx_6_p == 999] <- NA
names(dhx01)[names(dhx01) == 'devhx_6_p'] <- 'planned_preg'



#abcd_mx01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_mx01)[!(names(abcd_mx01) %in% not_numeric)]
for(i in 2:ncol(abcd_mx01)){
  if(names(abcd_mx01)[i] %in% numeric_cols){
    abcd_mx01[,i] <- as.numeric(abcd_mx01[,i])
  }
}
names(abcd_mx01)[names(abcd_mx01) == 'medhx_2b'] <- 'allergies'
names(abcd_mx01)[names(abcd_mx01) == 'medhx_6i'] <- 'head_injury'



#abcd_ssphp01
abcd_ssphp01$sds_p_ss_total <- as.numeric(abcd_ssphp01$sds_p_ss_total)
names(abcd_ssphp01)[names(abcd_ssphp01) == 'sds_p_ss_total'] <- 'sleep_prob'



#abcd_sscey01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_sscey01)[!(names(abcd_sscey01) %in% not_numeric)]
for(i in 2:ncol(abcd_sscey01)){
  if(names(abcd_sscey01)[i] %in% numeric_cols){
    abcd_sscey01[,i] <- as.numeric(abcd_sscey01[,i])
  }
}



#abcd_fhxssp01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_fhxssp01)[!(names(abcd_fhxssp01) %in% not_numeric)]
for(i in 2:ncol(abcd_fhxssp01)){
  if(names(abcd_fhxssp01)[i] %in% numeric_cols){
    abcd_fhxssp01[,i] <- as.numeric(abcd_fhxssp01[,i])
  }
}



#abcd_asrs01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_asrs01)[!(names(abcd_asrs01) %in% not_numeric)]
for(i in 2:ncol(abcd_asrs01)){
  if(names(abcd_asrs01)[i] %in% numeric_cols){
    abcd_asrs01[,i] <- as.numeric(abcd_asrs01[,i])
  }
}
treatment <- data.frame(subjectkey = abcd_asrs01$subjectkey, parent_adhd = abcd_asrs01$asr_scr_adhd_t)
ind <- which(names(abcd_asrs01) == 'asr_scr_adhd_t')
abcd_asrs01 <- abcd_asrs01[, -ind]
abcd_asrs01$asr_scr_anxdisord_t[abcd_asrs01$asr_scr_anxdisord_t < 65] <- 0
abcd_asrs01$asr_scr_anxdisord_t[abcd_asrs01$asr_scr_anxdisord_t >= 65] <- 1
abcd_asrs01$asr_scr_antisocial_t[abcd_asrs01$asr_scr_antisocial_t < 65] <- 0
abcd_asrs01$asr_scr_antisocial_t[abcd_asrs01$asr_scr_antisocial_t >= 65] <- 1
treatment$parent_adhd[which(treatment$parent_adhd < 65)] <- 0
treatment$parent_adhd[which(treatment$parent_adhd >= 65)] <- 1



#abcd_rhds01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_rhds01)[!(names(abcd_rhds01) %in% not_numeric)]
for(i in 2:ncol(abcd_rhds01)){
  if(names(abcd_rhds01)[i] %in% numeric_cols){
    abcd_rhds01[,i] <- as.numeric(abcd_rhds01[,i])
  }
}
names(abcd_rhds01)[names(abcd_rhds01) == 'reshist_addr1_adi_wsum'] <- 'adi'
names(abcd_rhds01)[names(abcd_rhds01) == 'reshist_addr1_pm252016aa'] <- 'pollution'



#abcd_ant01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_ant01)[!(names(abcd_ant01) %in% not_numeric)]
for(i in 2:ncol(abcd_ant01)){
  if(names(abcd_ant01)[i] %in% numeric_cols){
    abcd_ant01[,i] <- as.numeric(abcd_ant01[,i])
  }
}
abcd_ant01$anthroweight1lb[abcd_ant01$anthroweight1lb < 20] <- NA
abcd_ant01$anthro_1_height_in[abcd_ant01$anthro_1_height_in < 20] <- NA
for(i in 1:nrow(abcd_ant01)){
  abcd_ant01$bmi[i] <- abcd_ant01$anthroweight1lb[i] / (abcd_ant01$anthro_1_height_in[i])^2 * 703
}
ind1 <- which(names(abcd_ant01) == 'anthroweight1lb')
ind2 <- which(names(abcd_ant01) == 'anthro_1_height_in')
abcd_ant01 <- abcd_ant01[, -c(ind1, ind2)]



#abcd_cbcls01
not_numeric <- c('subjectkey')
numeric_cols <- names(abcd_cbcls01)[!(names(abcd_cbcls01) %in% not_numeric)]
for(i in 2:ncol(abcd_cbcls01)){
  if(names(abcd_cbcls01)[i] %in% numeric_cols){
    abcd_cbcls01[,i] <- as.numeric(abcd_cbcls01[,i])
  }
}
ind <- which(names(abcd_cbcls01) == 'cbcl_scr_dsm5_adhd_t')
outcome <- data.frame(subjectkey = abcd_cbcls01$subjectkey, child_adhd = abcd_cbcls01$cbcl_scr_dsm5_adhd_t)
abcd_cbcls01 <- abcd_cbcls01[, -ind]



#abcd_cbcls01_long
abcd_cbcls01_long$cbcl_scr_dsm5_adhd_t <- as.numeric(abcd_cbcls01_long$cbcl_scr_dsm5_adhd_t)
names(abcd_cbcls01_long)[names(abcd_cbcls01_long) == 'cbcl_scr_dsm5_adhd_t'] <- 'child_adhd'
outcome_1y <- abcd_cbcls01_long %>%
  filter(eventname == '1_year_follow_up_y_arm_1') %>%
  select(!all_of('eventname'))
outcome_2y <- abcd_cbcls01_long %>%
  filter(eventname == '2_year_follow_up_y_arm_1') %>%
  select(!all_of('eventname'))
outcome_3y <- abcd_cbcls01_long %>%
  filter(eventname == '3_year_follow_up_y_arm_1') %>%
  select(!all_of('eventname'))

incl_long <- full_join(outcome_1y, outcome_2y, by='subjectkey') %>%
  full_join(outcome_3y, by='subjectkey') %>%
  drop_na()
incl_long <- data.frame(subjectkey=incl_long$subjectkey)


#abcd_lt01
abcd_lt01$site_id_l <- as.factor(abcd_lt01$site_id_l)



#PRS
ind <- which(names(PRS) %in% c('sex'))
PRS <- PRS[, -ind]
PRS <- dummy_cols(PRS, select_columns = 'ethnic_g') %>%
  subset(select=-c(ethnic_g, ethnic_g_NA))



#PRS_holdout
#ind <- which(names(PRS_holdout) %in% c('sex','ethnic_g'))
#PRS_holdout <- PRS_holdout[, -ind]



############################## @ Brain data ##############################
inclusion_crit_vars <- c('subjectkey','imgincl_t1w_include','imgincl_t2w_include')
numeric_vars <- c('imgincl_t1w_include','imgincl_t2w_include')

inclusion_crit <- read.table('abcd_imgincl01.txt', header = T, sep = "\t", dec = ".")  %>%
  filter(eventname == 'baseline_year_1_arm_1') %>%
  select(all_of(inclusion_crit_vars)) %>%
  mutate_at(numeric_vars, as.numeric) %>%
  filter(imgincl_t1w_include == 1 & imgincl_t2w_include == 1)
smri_incl_subjects <- data.frame(subjectkey = inclusion_crit[,1])

smri <- read.table('abcd_smrip10201.txt', header = T) %>%
  filter(eventname == 'baseline_year_1_arm_1')
smri_incl <- left_join(smri_incl_subjects, smri, by='subjectkey') 

smri_vol <- data.frame(subjectkey=smri_incl$subjectkey)
smri_vol <- cbind(smri_vol, smri_incl[, which(str_detect(names(smri_incl), 'vol'))])

# exclude non-ROI metrics
not_needed <- c('smri_vol_cdk_totallh','smri_vol_cdk_totalrh','smri_vol_cdk_total',
                'smri_vol_scs_ltventriclelh','smri_vol_scs_inflatventlh','smri_vol_scs_3rdventricle',
                'smri_vol_scs_4thventricle','smri_vol_scs_csf','smri_vol_scs_lesionlh',
                'smri_vol_scs_ltventriclerh','smri_vol_scs_inflatventrh','smri_vol_scs_lesionrh',
                'smri_vol_scs_wmhint','smri_vol_scs_wmhintlh','smri_vol_scs_wmhintrh',
                'smri_vol_scs_wholeb','smri_vol_scs_latventricles','smri_vol_scs_allventricles')
ind <- which(names(smri_vol) %in% not_needed)
smri_vol <- smri_vol[,-ind]
ind <- which(str_detect(names(smri_vol), 'smri_vol_cf'))
smri_vol <- smri_vol[,-ind]


# convert to numeric
for(i in 2:ncol(smri_vol)){   
  smri_vol[,i] <- as.numeric(smri_vol[,i])
}

# adjust for total intracranial volume
for(i in 2:ncol(smri_vol)){   
  smri_vol[,i] <- smri_vol[,i] / smri_vol$smri_vol_scs_intracranialv
}

# remove near-zero variance features
filtered_col <- nearZeroVar(smri_vol, saveMetrics = TRUE) %>%
  rownames_to_column() %>%
  filter(nzv) 
not_needed <- filtered_col$rowname
smri_filt <- smri_vol %>% select(-all_of(not_needed))




############################## @ Final data ##############################
full_data <- full_join(treatment, outcome, by='subjectkey') %>%
  full_join(abcd_lt01,by='subjectkey') %>%
  full_join(abcd_ant01, by='subjectkey') %>%
  full_join(abcd_asrs01, by='subjectkey') %>%
  full_join(abcd_cbcls01,by='subjectkey') %>%
  full_join(abcd_fhxssp01,by='subjectkey') %>%
  full_join(abcd_mx01,by='subjectkey') %>%
  full_join(abcd_rhds01,by='subjectkey') %>%
  full_join(abcd_sscey01,by='subjectkey') %>%
  full_join(abcd_ssphp01, by='subjectkey') %>%
  #full_join(acspsw03,by='subjectkey') %>%
  full_join(dhx01,by='subjectkey') %>%
  full_join(pdem02_1,by='subjectkey') %>%
  full_join(pdem02_2,by='subjectkey') %>%
  full_join(pdem02_3,by='subjectkey') %>%
  full_join(PRS,by='subjectkey') %>%
  full_join(smri_filt,by='subjectkey') 


env_cov <- c('subjectkey','parent_adhd','child_adhd','site_id_l',
             'interview_age','sex','ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR','demo_comb_income_v2',
             'preeclampsia','preg_comp','alcohol','tobacco','marijuana','coc_crack',
             'preterm_birth','birth_weight_lbs','neo_illness','breastfeeding',
             'allergies','sleep_prob','head_injury',
             'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6',
             'crpbi_y_ss_parent','fes_y_ss_fc','famhx_ss_momdad_dprs_p','asr_scr_anxdisord_t','famhx_ss_momdad_dg_p','asr_scr_antisocial_t',
             'maternal_age','paternal_age','demo_prnt_ed_v2','pollution')

ind <- names(full_data)[names(full_data) %in% env_cov]
data_env <- full_data[, ind]

env_gen_cov <- unique(c(env_cov, names(PRS)))
ind <- names(full_data)[names(full_data) %in% env_gen_cov]
data_env_gen <- full_data[, ind]

ind <- which(names(full_data) %in% names(PRS)[-c(1, 32:35)])
data_env_br <- full_data[, -ind]

data_env_gen_br <- full_data



############################## @ Missing data imputation + downsampling ##############################
data <- data_env_gen_br
for(i in 2:ncol(data)){
  print(names(data)[i])
  print(sum(is.na(data[, i])))
}

ind1 <- which(is.na(data$parent_adhd))
ind2 <- which(is.na(data$child_adhd))
ind3 <- which(is.na(data$site_id_l))
ind4 <- which(is.na(data$interview_age))
data <- data[-unique(c(ind1, ind2, ind3, ind4)) ,]

data_tr <- data[which(data$parent_adhd == 1) ,]
num_m <- sum(data_tr$sex == 0)
num_f <- sum(data_tr$sex == 1)
ind <- which(names(data_tr) %in% c('subjectkey','site_id_l'))
data_tr_imp <- impute.knn(as.matrix(data_tr[,-ind]) ,k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
data_tr_imp <- as.data.frame(data_tr_imp$data)
data_tr_imp <- cbind(subjectkey = data_tr$subjectkey, site_id_l = data_tr$site_id_l, data_tr_imp)

ind <- which(data_tr_imp$ethnic_g_AFR  <= 0.5 & data_tr_imp$ethnic_g_AMR <= 0.5 & data_tr_imp$ethnic_g_EAS <= 0.5 & data_tr_imp$ethnic_g_EUR <= 0.5)
for(i in ind){
  print(i)
  print(data_tr_imp$ethnic_g_AFR[i])
  print(data_tr_imp$ethnic_g_AMR[i])
  print(data_tr_imp$ethnic_g_EAS[i])
  print(data_tr_imp$ethnic_g_EUR[i])
}

ind1 <- which(names(data_tr_imp) %in% c('ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR')) 
for(i in ind){
  max_ethn <- max(c(data_tr_imp$ethnic_g_AFR[i], data_tr_imp$ethnic_g_AMR[i], data_tr_imp$ethnic_g_EAS[i], data_tr_imp$ethnic_g_EUR[i]))
  max_ind <- which(data_tr_imp[i,ind1] == max_ethn)
  max_ind <- sample(max_ind, 1)
  data_tr_imp[i,ind1[max_ind]] <- 1
}

ind <- which(data_tr_imp$demo_prnt_marital_v2_1  <= 0.5 & data_tr_imp$demo_prnt_marital_v2_2 <= 0.5 & data_tr_imp$demo_prnt_marital_v2_3 <= 0.5 & 
               data_tr_imp$demo_prnt_marital_v2_4  <= 0.5 & data_tr_imp$demo_prnt_marital_v2_5 <= 0.5 & data_tr_imp$demo_prnt_marital_v2_6 <= 0.5)
for(i in ind){
  print(i)
  print(data_tr_imp$demo_prnt_marital_v2_1[i])
  print(data_tr_imp$demo_prnt_marital_v2_2[i])
  print(data_tr_imp$demo_prnt_marital_v2_3[i])
  print(data_tr_imp$demo_prnt_marital_v2_4[i])
  print(data_tr_imp$demo_prnt_marital_v2_5[i])
  print(data_tr_imp$demo_prnt_marital_v2_6[i])
}

ind1 <- which(names(data_tr_imp) %in% c('demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3',
                                        'demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')) 
for(i in ind){
  max_status <- max(c(data_tr_imp$demo_prnt_marital_v2_1[i], data_tr_imp$demo_prnt_marital_v2_2[i], data_tr_imp$demo_prnt_marital_v2_3[i], 
                    data_tr_imp$demo_prnt_marital_v2_4[i], data_tr_imp$demo_prnt_marital_v2_5[i], data_tr_imp$demo_prnt_marital_v2_6[i]))
  max_ind <- which(data_tr_imp[i,ind1] == max_status)
  max_ind <- sample(max_ind, 1)
  data_tr_imp[i,ind1[max_ind]] <- 1
}
bin_vars <- c('asr_scr_anxdisord_t','asr_scr_antisocial_t','famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p',
              'allergies','head_injury',
              'preeclampsia','preterm_birth','planned_preg',
              'alcohol','tobacco','marijuana','coc_crack','sex',
              'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6',
              'ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR')

ind <- which(names(data_tr_imp) %in% bin_vars)
for(i in ind){
  data_tr_imp[,i] <- round(data_tr_imp[,i], 0)
}
for(i in 3:46){
  print(names(data_tr_imp)[i])
  print(table(data_tr_imp[,i]))
}


data_ct <- data[which(data$parent_adhd == 0) ,]
data_ct <- left_join(incl_long, data_ct, by='subjectkey')
data_ct <- data_ct %>% drop_na()
data_ct_m <- data_ct[which(data_ct$sex == 0),]
data_ct_f <- data_ct[which(data_ct$sex == 1),]
set.seed(1)
data_ct_m <- data_ct_m[sample(nrow(data_ct_m),num_m), ]
data_ct_f <- data_ct_f[sample(nrow(data_ct_f),num_f), ]
data_ct <- rbind(data_ct_m, data_ct_f)
which(names(data_ct) == 'site_id_l')
data_ct <- data_ct[, c(1,4,2,3,5:ncol(data_ct))]

names(data_ct) == names(data_tr_imp)
#names(data_ct) == names(data_tr)

#data.imp <- rbind(data_tr_imp, data_ct)
#data.orig <- rbind(data_tr, data_ct)

fin_data_full <- rbind(data_tr_imp, data_ct)

fin_data_env <- fin_data_full[,names(fin_data_full) %in% names(data_env)]
fin_data_env_gen <- fin_data_full[,names(fin_data_full) %in% names(data_env_gen)]

fin_data_env_gen_br <- fin_data_full[fin_data_full$subjectkey %in% smri_incl_subjects$subjectkey,]
data_ct <- fin_data_env_gen_br[which(fin_data_env_gen_br$parent_adhd == 0),]
data_tr <- fin_data_env_gen_br[which(fin_data_env_gen_br$parent_adhd == 1),]
set.seed(1)
data_ct <- data_ct[sample(nrow(data_ct),nrow(data_tr)), ]
fin_data_env_gen_br <- rbind(data_ct, data_tr)

fin_data_env_br <- fin_data_env_gen_br[,names(fin_data_env_gen_br) %in% names(data_env_br)]



############################## @ GRF analysis ##############################
data <- fin_data_env_gen_br
not_numeric <- c('subjectkey','parent_adhd','site_id_l',
  'asr_scr_anxdisord_t','asr_scr_antisocial_t','famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p',
  'allergies','head_injury',
  'preeclampsia','preterm_birth','planned_preg',
  'alcohol','tobacco','marijuana','coc_crack','sex',
  'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6',
  'ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR')
numeric_cols <- names(data)[!(names(data) %in% not_numeric)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% numeric_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}



start_time <- Sys.time()

results <- list()
Y <- data$child_adhd
W <- data$parent_adhd
X <- subset(data, select=-c(subjectkey, child_adhd, parent_adhd, site_id_l))
site_id <- data$site_id_l


# grow five forests using the chosen sets of seeds
num_seeds <- 5
#seeds <- results_env_br_pca$seeds
seeds <- sample(1:10000, size=num_seeds, replace=FALSE)
print(seeds)
results[[paste0('seeds')]] <- seeds
forest_list <- list()

for(seed_ind in 1:length(seeds)){
  
  print(paste('seed_ind', seed_ind, sep=' = '))
  
  forest <- causal_forest(X, Y, W, Y.hat=NULL, W.hat=NULL, 
                          compute.oob.predictions = FALSE, 
                          clusters = site_id,
                          num.trees = 5000, 
                          tune.parameters = 'all',
                          num.threads = 16,
                          seed = seeds[seed_ind])
  forest_list[[paste0('seed', seed_ind)]] <- forest
}


forest_merged <- merge_forests(forest_list, compute.oob.predictions = TRUE)


# ATE estimation
ate <- average_treatment_effect(forest_merged, method='AIPW', target.sample='all')
ate_t <- ate[1] / ate[2]
ate_p <- 1.96 * (1 - pnorm(abs(ate_t)))
ate_ci_low <- ate[1] - 1 * qnorm(0.975) * ate[2] 
ate_ci_high <- ate[1] + 1 * qnorm(0.975) * ate[2] 
ate_results <- data.frame(ate[1], ate[2], ate_t, ate_p, ate_ci_low, ate_ci_high)
rownames(ate_results) <- 'ATE'
colnames(ate_results)<-c("estimate", "std.error", "t-score", "P-value", "[95% ", " CI]")
print(ate_results)
results[[paste0('ate')]] <- ate_results


# test calibration
test_cal <- test_calibration(forest_merged)
print(test_cal)
results[[paste0('test_cal')]] <- test_cal


# feature importance
var_imp <- c(variable_importance(forest_merged))
names(var_imp) <- colnames(X)
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
sorted_var_imp <- data.frame(importance=sorted_var_imp)
print(sorted_var_imp)
results[[paste0('var_imp')]] <- sorted_var_imp


# best linear projection
blp <- best_linear_projection(forest_merged, X)
#print(blp)
results[[paste0('blp')]] <- blp


# save results 
results_env_gen_br <- results

rm(forest, forest_list, forest_merged)



############################## @ Model optimization ##############################
# compute variable importance significance threshold
start_time <- Sys.time()

var_imp_sign <- list()
num_iter <- 500

for(iter in 1:num_iter){
  print(paste('iter', iter, sep=' = '))
  
  W_perm <- permute(W)
  
  # grow a causal forest
  forest_list <- list()
  
  forest <- causal_forest(X, Y, W_perm, Y.hat=NULL, W.hat=NULL, 
                          compute.oob.predictions = TRUE, 
                          clusters = site_id,
                          num.trees = 5000, 
                          tune.parameters = 'all',
                          num.threads = 16,
                          seed = sample(1:10000, size=1))
  
  # feature importance
  var_imp_iter <- c(variable_importance(forest))
  var_imp_sign[[paste0('iter=', iter)]] <- var_imp_iter
  
}

var_imp_vec <- unlist(var_imp_sign)
threshold <- quantile(var_imp_vec, probs=0.95)

end_time <- Sys.time()
end_time - start_time


# save variable importance significance threshold
results[[paste0('threshold')]] <- threshold


rm(forest, forest_merged, forest_list)


# save important variables
imp_names <- rownames(sorted_var_imp)[sorted_var_imp$importance > threshold]
#blp <- as.data.frame(results$blp[,])
#imp_names <- append(imp_names, rownames(blp)[which(blp$`Pr(>|t|)` < 0.05)])
#imp_names <- unique(imp_names)
results[[paste0('imp_names')]] <- imp_names

results_env_br_pca <- results



############################## @ Longitudinal outcomes ##############################
data <- fin_data_env_gen

ind <- which(names(data) == 'child_adhd')
fin_data_env_gen_1y <- data[,-ind]
fin_data_env_gen_1y <- left_join(fin_data_env_gen_1y, outcome_1y, by='subjectkey')
fin_data_env_gen_1y <- fin_data_env_gen_1y %>% drop_na()

fin_data_env_gen_2y <- data[,-ind]
fin_data_env_gen_2y <- left_join(fin_data_env_gen_2y, outcome_2y, by='subjectkey')
fin_data_env_gen_2y <- fin_data_env_gen_2y %>% drop_na()

fin_data_env_gen_3y <- data[,-ind]
fin_data_env_gen_3y <- left_join(fin_data_env_gen_3y, outcome_3y, by='subjectkey')
fin_data_env_gen_3y <- fin_data_env_gen_3y %>% drop_na()



"
data_env_gen_1y <- data_1y[which(data_1y$subjectkey %in% fin_data_env_gen$subjectkey), 
                           which(names(data_1y) %in% names(data_env_gen))]
data_env_gen_2y <- data_2y[which(data_2y$subjectkey %in% fin_data_env_gen$subjectkey), 
                           which(names(data_2y) %in% names(data_env_gen))]
data_env_gen_3y <- data_3y[which(data_3y$subjectkey %in% fin_data_env_gen$subjectkey), 
                           which(names(data_3y) %in% names(data_env_gen))]

data <- data_env_gen_3y

ind1 <- which(is.na(data$parent_adhd))
ind2 <- which(is.na(data$child_adhd))
ind3 <- which(is.na(data$site_id_l))
ind4 <- which(is.na(data$interview_age))
data <- data[-unique(c(ind1, ind2, ind3, ind4)) ,]

data_tr <- data[which(data$parent_adhd == 1) ,]
num_m <- sum(data_tr$sex == 0)
num_f <- sum(data_tr$sex == 1)
ind <- which(names(data_tr) %in% c('subjectkey','site_id_l'))
data_tr_imp <- impute.knn(as.matrix(data_tr[,-ind]) ,k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
data_tr_imp <- as.data.frame(data_tr_imp$data)
data_tr_imp <- cbind(subjectkey = data_tr$subjectkey, site_id_l = data_tr$site_id_l, data_tr_imp)
bin_vars <- c('famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
              'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
              'preeclampsia','preterm_birth','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
              'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
ind <- which(names(data_tr_imp) %in% bin_vars)
for(i in ind){
  data_tr_imp[,i] <- round(data_tr_imp[,i], 0)
}

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = 'lt.blue')
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = 'lt.pink')
hist(data_tr_imp$neo_illness, col=c1)
hist(data_tr$neo_illness, col=c2, add=T)

data_ct <- data[which(data$parent_adhd == 0) ,]
data_ct <- data_ct %>% drop_na()
data_ct_m <- data_ct[which(data_ct$sex == 0),]
data_ct_f <- data_ct[which(data_ct$sex == 1),]
set.seed(1)
data_ct_m <- data_ct_m[sample(nrow(data_ct_m),num_m), ]
data_ct_f <- data_ct_f[sample(nrow(data_ct_f),num_f), ]
data_ct <- rbind(data_ct_m, data_ct_f)
which(names(data_ct) == 'site_id_l')
data_ct <- data_ct[, c(1,3,2,4:ncol(data_ct))]

names(data_ct) == names(data_tr_imp)

fin_data_env_gen_3y <- rbind(data_tr_imp, data_ct)

ind <- which(names(fin_data_env_gen_eur) == 'child_adhd')
fin_data_env_gen_eur_1y <- fin_data_env_gen_eur[,-ind]
fin_data_env_gen_eur_2y <- fin_data_env_gen_eur[,-ind]
fin_data_env_gen_eur_3y <- fin_data_env_gen_eur[,-ind]

fin_data_env_gen_eur_1y <- full_join(outcome_1y, fin_data_env_gen_eur_1y, by='subjectkey') %>% drop_na()
fin_data_env_gen_eur_2y <- full_join(outcome_2y, fin_data_env_gen_eur_2y, by='subjectkey') %>% drop_na()
fin_data_env_gen_eur_3y <- full_join(outcome_3y, fin_data_env_gen_eur_3y, by='subjectkey') %>% drop_na()
"


############################## @ Holdout dataset ##############################
ind <- which(!names(PRS_holdout) %in% names(PRS))
PRS_holdout <- PRS_holdout[, -ind]
data_env_gen_hold <- full_data[which(full_data$subjectkey %in% PRS_holdout$subjectkey),
                               which(names(full_data) %in% names(data_env_gen))]
ind <- which(names(data_env_gen_hold) %in% names(PRS_holdout)[-1])
data_env_gen_hold <- data_env_gen_hold[, -ind]
data_env_gen_hold <- left_join(data_env_gen_hold, PRS_holdout, by='subjectkey')


data <- data_env_gen_hold

ind1 <- which(is.na(data$parent_adhd))
ind2 <- which(is.na(data$child_adhd))
ind3 <- which(is.na(data$site_id_l))
ind4 <- which(is.na(data$interview_age))
data <- data[-unique(c(ind1, ind2, ind3, ind4)) ,]

data_tr <- data[which(data$parent_adhd == 1) ,]
num_m <- sum(data_tr$sex == 0)
num_f <- sum(data_tr$sex == 1)
ind <- which(names(data_tr) %in% c('subjectkey','site_id_l'))
data_tr_imp <- impute.knn(as.matrix(data_tr[,-ind]) ,k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
data_tr_imp <- as.data.frame(data_tr_imp$data)
data_tr_imp <- cbind(subjectkey = data_tr$subjectkey, site_id_l = data_tr$site_id_l, data_tr_imp)
bin_vars <- c('famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
              'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
              'preeclampsia','preterm_birth','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
              'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
ind <- which(names(data_tr_imp) %in% bin_vars)
for(i in ind){
  data_tr_imp[,i] <- round(data_tr_imp[,i], 0)
}

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(data_tr_imp$neo_illness, col=c1)
hist(data_tr$neo_illness, col=c2, add=T)

data_ct <- data[which(data$parent_adhd == 0) ,]
data_ct <- data_ct %>% drop_na()
"
data_ct_m <- data_ct[which(data_ct$sex == 0),]
data_ct_f <- data_ct[which(data_ct$sex == 1),]
set.seed(1)
data_ct_m <- data_ct_m[sample(nrow(data_ct_m),num_m), ]
data_ct_f <- data_ct_f[sample(nrow(data_ct_f),num_f), ]
data_ct <- rbind(data_ct_m, data_ct_f)
"
which(names(data_ct) == 'site_id_l')
data_ct <- data_ct[, c(1,4,2,3,5:ncol(data_ct))]

names(data_ct) == names(data_tr_imp)

fin_data_env_gen_hold <- rbind(data_tr_imp, data_ct)

data <- fin_data_env_gen_hold
for(i in 2:ncol(data)){
  print(names(data)[i])
  print(sum(is.na(data[,i])))
}



############################## @ Downsampling + No imputation ##############################
data <- data_env_gen_br

data <- left_join(smri_incl_subjects, data, by='subjectkey')

ind1 <- which(is.na(data$parent_adhd))
ind2 <- which(is.na(data$child_adhd))
ind3 <- which(is.na(data$site_id_l))
ind4 <- which(is.na(data$interview_age))
data <- data[-unique(c(ind1, ind2, ind3, ind4)) ,]

data_tr <- data[which(data$parent_adhd == 1) ,]
data_tr <- data_tr %>% drop_na()
num_m <- sum(data_tr$sex == 0)
num_f <- sum(data_tr$sex == 1)

data_ct <- data[which(data$parent_adhd == 0) ,]
data_ct <- data_ct %>% drop_na()
data_ct_m <- data_ct[which(data_ct$sex == 0),]
data_ct_f <- data_ct[which(data_ct$sex == 1),]
set.seed(1)
data_ct_m <- data_ct_m[sample(nrow(data_ct_m),num_m), ]
data_ct_f <- data_ct_f[sample(nrow(data_ct_f),num_f), ]
data_ct <- rbind(data_ct_m, data_ct_f)

names(data_ct) == names(data_tr)

fin_data_noimp_env_gen_br <- rbind(data_tr, data_ct)



############################## @ No downsampling + Imputation ##############################
data <- fin_data_env_gen_br_eur_nodown
for(i in 2:ncol(data)){
  print(names(data)[i])
  print(sum(is.na(data[, i])))
}

data <- inner_join(smri_incl_subjects, data, by='subjectkey')

ind1 <- which(is.na(data$parent_adhd))
ind2 <- which(is.na(data$child_adhd))
ind3 <- which(is.na(data$site_id_l))
ind4 <- which(is.na(data$interview_age))
data <- data[-unique(c(ind1, ind2, ind3, ind4)) ,]

data_tr <- data[which(data$parent_adhd == 1) ,]
ind <- which(names(data_tr) %in% c('subjectkey','site_id_l'))
data_tr_imp <- impute.knn(as.matrix(data_tr[,-ind]) ,k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed=362436069)
data_tr_imp <- as.data.frame(data_tr_imp$data)
data_tr_imp <- cbind(subjectkey = data_tr$subjectkey, site_id_l = data_tr$site_id_l, data_tr_imp)
bin_vars <- c('famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
              'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
              'preeclampsia','preterm_birth','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
              'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2','demo_prnt_marital_v2_3','demo_prnt_marital_v2_4','demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
ind <- which(names(data_tr_imp) %in% bin_vars)
for(i in ind){
  data_tr_imp[,i] <- round(data_tr_imp[,i], 0)
}

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(data_tr_imp$neo_illness, col=c1)
hist(data_tr$neo_illness, col=c2, add=T)

data_ct <- data[which(data$parent_adhd == 0) ,]
data_ct <- data_ct %>% drop_na()
which(names(data_ct) == 'site_id_l')
data_ct <- data_ct[, c(1,4,2,3,5:ncol(data_ct))]

names(data_ct) == names(data_tr_imp)

fin_data_env_gen_br_eur_nodown <- rbind(data_tr_imp, data_ct)


ind <- which(names(fin_data_nodown_env) == 'child_adhd')
fin_data_nodown_env_1y <- fin_data_nodown_env[, -ind]
fin_data_nodown_env_2y <- fin_data_nodown_env[, -ind]
fin_data_nodown_env_3y <- fin_data_nodown_env[, -ind]

fin_data_nodown_env_1y <- left_join(outcome_1y, fin_data_nodown_env_1y, by='subjectkey') %>% drop_na()
fin_data_nodown_env_2y <- left_join(outcome_2y, fin_data_nodown_env_2y, by='subjectkey') %>% drop_na()
fin_data_nodown_env_3y <- left_join(outcome_3y, fin_data_nodown_env_3y, by='subjectkey') %>% drop_na()


ind <- which(names(fin_data_nodown_env_gen) == 'child_adhd')
fin_data_nodown_env_gen_1y <- fin_data_nodown_env_gen[, -ind]
fin_data_nodown_env_gen_2y <- fin_data_nodown_env_gen[, -ind]
fin_data_nodown_env_gen_3y <- fin_data_nodown_env_gen[, -ind]

fin_data_nodown_env_gen_1y <- left_join(outcome_1y, fin_data_nodown_env_gen_1y, by='subjectkey') %>% drop_na()
fin_data_nodown_env_gen_2y <- left_join(outcome_2y, fin_data_nodown_env_gen_2y, by='subjectkey') %>% drop_na()
fin_data_nodown_env_gen_3y <- left_join(outcome_3y, fin_data_nodown_env_gen_3y, by='subjectkey') %>% drop_na()





############################## @ No downsampling + No imputation ##############################
fin_data_full_env <- data_env %>% drop_na()
fin_data_full_env_gen <- data_env_gen %>% drop_na()

fin_data_full_env_br <- left_join(smri_incl_subjects, data_env_br, by='subjectkey') %>% drop_na()
fin_data_full_env_gen_br <- left_join(smri_incl_subjects, data_env_gen_br, by='subjectkey') %>% drop_na()

ind <- which(names(fin_data_full_env_br) == 'child_adhd')
fin_data_full_env_br_1y <- fin_data_full_env_br[, -ind]
fin_data_full_env_br_2y <- fin_data_full_env_br[, -ind]
fin_data_full_env_br_3y <- fin_data_full_env_br[, -ind]

fin_data_full_env_br_1y <- left_join(outcome_1y, fin_data_full_env_br_1y, by='subjectkey') %>% drop_na()
fin_data_full_env_br_2y <- left_join(outcome_2y, fin_data_full_env_br_2y, by='subjectkey') %>% drop_na()
fin_data_full_env_br_3y <- left_join(outcome_3y, fin_data_full_env_br_3y, by='subjectkey') %>% drop_na()



############################## @ Outcome difference distribution ##############################
baseline_1y <- full_join(outcome_1y, outcome, by='subjectkey') %>% drop_na()
names(baseline_1y) <- c('subjectkey','child_adhd_1y','child_adhd')
baseline_1y$diff <- baseline_1y$child_adhd - baseline_1y$child_adhd_1y
hist(baseline_1y$diff)
summary(baseline_1y$diff)

baseline_2y <- full_join(outcome_2y, outcome, by='subjectkey') %>% drop_na()
names(baseline_2y) <- c('subjectkey','child_adhd_1y','child_adhd')
baseline_2y$diff <- baseline_2y$child_adhd - baseline_2y$child_adhd_1y
hist(baseline_2y$diff)
summary(baseline_2y$diff)

baseline_3y <- full_join(outcome_3y, outcome, by='subjectkey') %>% drop_na()
names(baseline_3y) <- c('subjectkey','child_adhd_1y','child_adhd')
baseline_3y$diff <- baseline_3y$child_adhd - baseline_3y$child_adhd_1y
hist(baseline_3y$diff)
summary(baseline_3y$diff)


############################## @ Partial dependence ##############################
is.sorted <- function(num.vec){
  if(length(num.vec) > 1 & all(sort(num.vec, decreasing=F) == num.vec)){
    return('asc')
  } else if(length(num.vec) > 1 & all(sort(num.vec, decreasing=T) == num.vec)){
    return('des')
  } else {
    return('none')
  }
}

features <- results_env_gen_1y$imp_names
univ.sim <- data.frame(features)
names(univ.sim) <- "feature"
univ.sim[, "tau.hat.Q0"] <- 0
univ.sim[, "tau.hat.Q1"] <- 0
univ.sim[, "tau.hat.Q2"] <- 0
univ.sim[, "tau.hat.Q3"] <- 0
univ.sim[, "tau.hat.Q4"] <- 0
univ.sim[, "trend"] <- 0

data <- fin_data_env_gen_1y
not_numeric <- c('subjectkey','parent_adhd','site_id_l',
                 'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
                 'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
                 'preeclampsia','preterm_birth','planned_preg','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
                 'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2',
                 'demo_prnt_marital_v2_3','demo_prnt_marital_v2_4',
                 'demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
numeric_cols <- names(data)[!(names(data) %in% not_numeric)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% numeric_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}


for(i in 1:nrow(univ.sim)){
  selected.covariate <- univ.sim$feature[i]
  X <- subset(data, select=-c(subjectkey, child_adhd, parent_adhd, site_id_l))
  covariates <- names(X)
  other.covariates <- covariates[which(covariates != selected.covariate)]
  
  # Fitting a forest
  # For a continuous variable, select appropriate percentiles
  percentiles <- c(0, .25, .5, .75, 1)
  grid.size <- length(percentiles)
  covariate.grid <- quantile(data[,selected.covariate], probs=percentiles)
  
  # Take median of other covariates 
  medians <- apply(data[, other.covariates, F], 2, median)
  
  # Construct a dataset
  data.grid <- data.frame(sapply(medians, function(x) rep(x, grid.size)), covariate.grid)
  colnames(data.grid) <- c(other.covariates, selected.covariate)
  
  # Expand the data
  fmla <- formula(paste0("~ 0 + ", paste0(covariates, collapse="+")))
  X.grid <- model.matrix(fmla, data.grid)
  
  # Point predictions of the CATE and standard errors 
  forest.pred <- predict(results_env_gen$forest, newdata = X.grid, estimate.variance=TRUE)
  tau.hat <- forest.pred$predictions
  tau.hat.se <- sqrt(forest.pred$variance.estimates)
  tau.hat.tstat <- tau.hat / tau.hat.se
  tau.hat.p <- 1.96 * (1 - pnorm(abs(tau.hat.tstat)))
  tau.hat.p.fdr <- p.adjust(tau.hat.p)
  quantiles <- c('Q0','Q1','Q2','Q3','Q4')
  
  univ.sim[i,c(2:6)] <- tau.hat
  univ.sim$trend[i] <- is.sorted(tau.hat)
}



############################## @ Sensitivity analysis: PGS ethnicity ##############################
ethn_vars <- c('subjectkey','ethnic_g')
ethn <- read.csv('ABCD_EUR+Multibased_PRScsx_PC1-10resid_raw.csv', header = T, sep = ",", dec = ".") %>%
  select(-all_of('X')) %>%
  select(all_of(ethn_vars))

ethn_eur <- ethn %>% filter(ethnic_g == 'EUR')

fin_data_env_eur <- fin_data_env[which(fin_data_env$subjectkey %in% ethn_eur$subjectkey),]
fin_data_env_gen_eur <- fin_data_env_gen[which(fin_data_env_gen$subjectkey %in% ethn_eur$subjectkey),]
fin_data_env_br_eur <- fin_data_env_br[which(fin_data_env_br$subjectkey %in% ethn_eur$subjectkey),]
fin_data_env_gen_br_eur <- fin_data_env_gen_br[which(fin_data_env_gen_br$subjectkey %in% ethn_eur$subjectkey),]

# no downsampling + imputation
fin_data_env_eur_nodown <- data_env[which(data_env$subjectkey %in% ethn_eur$subjectkey),] 
fin_data_env_gen_eur_nodown <- data_env_gen[which(data_env_gen$subjectkey %in% ethn_eur$subjectkey),]
fin_data_env_br_eur_nodown <- data_env_br[which(data_env_br$subjectkey %in% ethn_eur$subjectkey),]
fin_data_env_gen_br_eur_nodown <- data_env_gen_br[which(data_env_gen_br$subjectkey %in% ethn_eur$subjectkey),]

# no downsampling + no imputation
fin_data_env_eur_nodown_noimp <- data_env[which(data_env$subjectkey %in% ethn_eur$subjectkey),] %>% drop_na()
fin_data_env_gen_eur_nodown_noimp <- data_env_gen[which(data_env_gen$subjectkey %in% ethn_eur$subjectkey),] %>% drop_na()
fin_data_env_br_eur_nodown_noimp <- data_env_br[which(data_env_br$subjectkey %in% ethn_eur$subjectkey),] %>% drop_na()
fin_data_env_gen_br_eur_nodown_noimp <- data_env_gen_br[which(data_env_gen_br$subjectkey %in% ethn_eur$subjectkey),] %>% drop_na()

fin_data_env_br_eur_nodown_noimp <- left_join(smri_incl_subjects, fin_data_env_br_eur_nodown_noimp, by='subjectkey') %>% drop_na()
fin_data_env_gen_br_eur_nodown_noimp <- left_join(smri_incl_subjects, fin_data_env_gen_br_eur_nodown_noimp, by='subjectkey') %>% drop_na()


############################## @ Sensitivity analysis: brain PCA ##############################
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library(factoextra)

data <- fin_data_env_gen_br
not_numeric <- c('subjectkey','parent_adhd','site_id_l',
                 'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
                 'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
                 'preeclampsia','preterm_birth','planned_preg','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
                 'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2',
                 'demo_prnt_marital_v2_3','demo_prnt_marital_v2_4',
                 'demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
numeric_cols <- names(data)[!(names(data) %in% not_numeric)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% numeric_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}

data_brain <- data[, which(names(data) %in% names(smri_filt)[-1])]
data_brain.pca <- princomp(data_brain)
summary(data_brain.pca)
fviz_eig(data_brain.pca, addlabels = TRUE)
fviz_cos2(data_brain.pca, choice = "var", axes = 1:2)

ind <- which(names(fin_data_env_br) == 'smri_vol_cdk_banksstslh')
fin_data_env_br_pca <- data.frame(cbind(fin_data_env_br[,1:(ind-1)], data_brain.pca$scores[,1:26]))
ind <- which(names(fin_data_env_gen_br) == 'smri_vol_cdk_banksstslh')
fin_data_env_gen_br_pca <- data.frame(cbind(fin_data_env_gen_br[,1:(ind-1)], data_brain.pca$scores[,1:26]))



############################## @ Sensitivity analysis: brain & PRS PCA ##############################
data <- fin_data_env_gen_br
not_numeric <- c('subjectkey','parent_adhd','site_id_l',
                 'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
                 'race_ethnicity_1','race_ethnicity_2','race_ethnicity_3','race_ethnicity_4','race_ethnicity_5',
                 'preeclampsia','preterm_birth','planned_preg','alcohol','tobacco','marijuana','coc_crack','neo_illness','sex',
                 'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2',
                 'demo_prnt_marital_v2_3','demo_prnt_marital_v2_4',
                 'demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')
numeric_cols <- names(data)[!(names(data) %in% not_numeric)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% numeric_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}

data_prs <- data[, which(names(data) %in% names(PRS)[-1])]
data_prs.pca <- princomp(data_prs)
summary(data_prs.pca)
for(i in 1:length(colnames(data_prs.pca$scores))){
  colnames(data_prs.pca$scores)[i] <- paste0(colnames(data_prs.pca$scores)[i], '.prs')
}

ind <- which(names(fin_data_env_gen_br) == 'ADHDeur6')
fin_data_env_gen_br_pca2 <- data.frame(cbind(fin_data_env_gen_br[,1:(ind-1)], data_prs.pca$scores[,1:14], data_brain.pca$scores[,1:26]))