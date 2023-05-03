# load libraries

library(rio)
library(vegan)

# import data

key <- import("raw/id_conversion.txt")
pheno1 <- import("raw/scapis_data_20190215_uppsala.txt")
pheno2 <- import("raw/scapis_data_20190215_malmo.txt")
tax_pheno <- import("raw/pheno_MGS_shannon_bray_curtis_MGP_4839_upp_4980_malmo.tsv")
cv <- import("raw/SCAPIS-DATA-PETITION-454-20221011.csv", na.strings = "")
hem1 <- import("raw/scapis neu lpk.dta")
hem2 <- import("raw/LPK_tove.xlsx")
diet <- import("raw/AdSugMacro_ScapisMöUppJanuary2021_UE220907_tSergiOct2022_.sav")
met <- import("raw/scapis_metabolon_labels_batchnorm.tsv")
drugs <- import("raw/SCAPIS-REGISTERSU-1715-LMED-20220511-T_R_LMED_27947_2021.txt")
king <- import("raw/bacteria_cacs_list_related_kingtable.kin0", format = "tsv")
ds1 <- readRDS("raw/lungut03.dsMgsCounts.rds")
ds2 <- readRDS("raw/upugut03.dsMgsCounts.rds")
out<-read.table("raw/bacteria_cacs_list_related_1st.king.cutoff.out.id")
GMM<-import("/proj/nobackup/sens2019512/users/baldanzi/pa_gut/processed/GMMrelativeabundance_9819.tsv")

# clean

key$pheno_id <- key$export_id
key$pheno_id[which(grepl("2-", key$subject_id))] <- key$subject_id[which(grepl("2-", key$subject_id))]
  
colnames(pheno1)[1] <- "pheno_id"
colnames(pheno2)[1] <- "pheno_id"
pheno <- rbind(pheno1, pheno2)
pheno$scapis_id <- key$subject_id[match(pheno$pheno_id, key$pheno_id)]

colnames(tax_pheno)[2] <- "pheno_id"
tax_pheno$scapis_id <- key$subject_id[match(tax_pheno$pheno_id, key$pheno_id)]

colnames(cv)[1] <- "export_id"
cv$scapis_id <- key$subject_id[match(cv$export_id, key$export_id)]
cv <- cv[which(cv$CCTA_Valid_4_Proximal == "YES"), ]

hem1 <- hem1[, c("subjid","lpk_res","neut_res")]
colnames(hem2) <- tolower(colnames(hem2))
hem <- rbind(hem1, hem2)
colnames(hem)[1] <- "scapis_id"

colnames(diet)[1] <- "scapis_id"

colnames(drugs)[1] <- "scapis_id"

king <- king[which(king$KINSHIP >= 0.177), ]
king <- king[which(!king$ID1 %in% king$ID2), ]
king$scapis_id <- king$ID2

ds <- rbind(ds1, ds2)
rownames(ds) <- gsub("_b", "", rownames(ds))
rownames(ds) <- gsub("b", "", rownames(ds))
rownames(ds) <- gsub("c", "", rownames(ds))
rownames(ds) <- tax_pheno$scapis_id[match(rownames(ds), tax_pheno$sample.id)]

# match data

id <- pheno$scapis_id[which(pheno$scapis_id %in% tax_pheno$scapis_id)]
pheno <- pheno[match(id, pheno$scapis_id), ]
tax_pheno <- tax_pheno[match(id, tax_pheno$scapis_id), ]
cv <- cv[match(id, cv$scapis_id), ]
hem <- hem[match(id, hem$scapis_id),]
diet <- diet[match(id, diet$scapis_id),]
met <- met[match(id, met$scapis_id), ]
king <- king[match(id, king$scapis_id), ]
ds <- ds[match(id, rownames(ds)), ]

# phenotypes

scapis_id <- pheno$scapis_id
cacstot <- pheno$casctot
cacstot_cat <- pheno$casctot
cacstot_cat[which(pheno$casctot == 0)] <- "0"
cacstot_cat[which(pheno$casctot > 0 & pheno$casctot <= 100)] <- "1-100"
cacstot_cat[which(pheno$casctot > 100 & pheno$casctot <= 400)] <- "101-400"
cacstot_cat[which(pheno$casctot > 400)] <- ">400"
athero <- ifelse(cv$Atherosclerosis_Any == "YES", "Yes", "No")
stenosis <- ifelse(cv$Above50_Any == "YES", "Yes", "No")
duke <- cv$DukeIndex
duke[which(!is.na(duke))] <- paste0("DUKE_", duke[which(!is.na(duke))])
sis <- cv$SIS
sis[which(!is.na(sis))] <- paste0("SIS_", sis[which(!is.na(sis))])
left <- ifelse(pheno$LcaTotPlaques > 0, 1, 0)
right <- ifelse(pheno$RcaTotPlaques > 0, 1, 0)
plaque <- rowSums(cbind(left, right))
plaque[which(rowSums(cbind(left, right)) == 0)] <- "None"
plaque[which(rowSums(cbind(left, right)) == 1)] <- "Unilateral"
plaque[which(rowSums(cbind(left, right)) == 2)] <- "Bilateral"
age <- pheno$agev1
sex <- ifelse(pheno$Gender == 1, "Male", "Female")
country <- tax_pheno$q005a
country[which(country == "scandinavia")] <- "Scandinavia"
country[which(country == "europe")] <- "Non-Scandinavian Europe"
country[which(country == "asia")] <- "Asia"
country[which(country == "other")] <- "Other"
smoke <- pheno$smokestatus
smoke[which(pheno$smokestatus == 1)] <- "Never smoker"
smoke[which(pheno$smokestatus == 2)] <- "Former smoker"
smoke[which(pheno$smokestatus == 3)] <- "Current smoker"
smoke[which(pheno$smokestatus == 9)] <- NA
pa <- pheno$q134
pa[which(pheno$q134 == 0)] <- "Sedentary"
pa[which(pheno$q134 == 1)] <- "Moderate exercise"
pa[which(pheno$q134 == 2)] <- "Moderate but regular exercise"
pa[which(pheno$q134 == 3)] <- "Regular exercise and training"
pa[which(pheno$q134 == 4)] <- NA
tg <- log(pheno$tg_res)
ldl <- pheno$ldl_res
hdl <- pheno$hdl_res
chol <- pheno$chol_res
sbp <- pheno$bbps
dbp <- pheno$bbpd
bmi <- pheno$bmi
crp <- log(pheno$hscrp_res)
neut <- log(hem$neut_res)
leuk <- log(hem$lpk_res)
energy <- diet$Enerkcal
carb <- diet$CHOepro
protein <- diet$ProEpro
fat <- diet$FettEpro
fiber <- diet$FibPer1000
diab <- pheno$q030lx
diab[which(pheno$q030lx == 0)] <- "No"
diab[which(pheno$q030lx == 1)] <- "Yes"
diab[which(pheno$q030lx == 2)] <- NA
crohn <- pheno$q030wx
crohn[which(pheno$q030wx == 0)] <- "No"
crohn[which(pheno$q030wx == 1)] <- "Yes"
crohn[which(pheno$q030wx == 2)] <- NA
chol_med <- pheno$q030k
chol_med[which(pheno$q030kx == 0)] <- "No"
chol_med[which(pheno$q030k == 0)] <- "No"
chol_med[which(pheno$q030k == 1)] <- "Yes"
chol_med[which(pheno$q030k == 2)] <- NA
bp_med <- pheno$q030j
bp_med[which(pheno$q030jx == 0)] <- "No"
bp_med[which(pheno$q030j == 0)] <- "No"
bp_med[which(pheno$q030j == 1)] <- "Yes"
bp_med[which(pheno$q030j == 2)] <- NA
diab_med <- pheno$q030l
diab_med[which(pheno$q030lx == 0)] <- "No"
diab_med[which(pheno$q030l == 0)] <- "No"
diab_med[which(pheno$q030l == 1)] <- "Yes"
diab_med[which(pheno$q030l == 2)] <- NA
ppi <- ifelse(met$MET_100002725 == "measured" | met$MET_100002808 == "measured", "Yes", "No")
narrow <- drugs[which(drugs$ATC %in% c("J01CE02","J01CF05","J01EA01","J01FA01","J01FA06","J01FA09","J01FA10","J01FF01","J01XC01","J01XE01")), ]
narrow$visit <- as.Date(pheno$Day1VisitD, "%d-%B-%Y")[match(narrow$scapis_id, pheno$scapis_id)]
narrow$edatum <- as.Date(narrow$EDATUM)
narrow <- narrow[which(narrow$visit - narrow$edatum < 366 & narrow$visit - narrow$edatum >= 0), ]
narrow <- ifelse(scapis_id %in% narrow$scapis_id, "Yes", "No")
broad <- drugs[which(drugs$ATC %in% c("J01AA02","J01AA04","J01AA06","J01AA07","J01CA04","J01CA08","J01CR02","J01DB05","J01DD14","J01EE01","J01MA02","J01MA06","J01MA12","J01MA14","J01XX05","J01XX08")), ]
broad$visit <- as.Date(pheno$Day1VisitD, "%d-%B-%Y")[match(broad$scapis_id, pheno$scapis_id)]
broad$edatum <- as.Date(broad$EDATUM)
broad <- broad[which(broad$visit - broad$edatum < 366 & broad$visit - broad$edatum >= 0), ]
broad <- ifelse(scapis_id %in% broad$scapis_id, "Yes", "No")
site <- ifelse(pheno$SITEID == 5, "Uppsala", "Malmo")
plate <- tax_pheno$plate
site_plate <- paste(site, plate, sep = "_")
family <- king$ID1
family[which(is.na(family))] <- scapis_id[which(is.na(family))]
family1=scapis_id
family1[which(family1%in%out[,1]==T)]<-0
shannon <- diversity(ds, index = "shannon")
simpson <- diversity(ds, index = "simpson")
invsimpson <- diversity(ds, index = "invsimpson")
chao <- estimateR(ds)["S.chao1", ]
chao[which(is.na(shannon))] <- NA
cvd <- apply(pheno[, c("q030ax", "q030bx", "q030cx", "q030dx", "q030ex", "q030fx", "q030gx", "q030ix")], 1, function(x) any(x == 1, na.rm = T))

cvd2 <- pheno[, c("q030ax", "q030bx", "q030cx", "q030dx", "q030ex", "q030fx", "q030gx", "q030ix")]
names(cvd2)=c("MI","Angina","Atrial.Fibrillation","Heart.Failure","Valvular.disease","Bypass.Surgery","Revascularization","Stroke")

cvd2.1=ifelse(cvd2$MI==1 |  cvd2$Angina==1 | cvd2$Bypass.Surgery==1 | cvd2$Revascularization==1,
1,0)
cvd2.2=ifelse(cvd2.1==1 | cvd2$Heart.Failure==1 |cvd2$Atrial.Fibrillation==1 | cvd2$Valvular.disease==1 | cvd2$Stroke==1,1,0)

# mgs

mgs <- tax_pheno[, which(grepl("HG3A", colnames(tax_pheno)))]
colnames(mgs) <- gsub("^.*___", "", colnames(mgs))

# prepare the table

data <- data.frame(scapis_id, cacstot, cacstot_cat, athero, stenosis, duke, sis, plaque, age, sex, country, smoke, pa, tg, ldl, hdl, chol, sbp, dbp, bmi, crp, neut, leuk, energy, carb, protein, fat, fiber, diab, crohn, chol_med, bp_med, diab_med, ppi, narrow, broad, site, site_plate, family,family1, shannon, simpson, invsimpson, chao, mgs)
data2<- data.frame(scapis_id, cacstot, cacstot_cat, athero, stenosis, duke, sis, plaque, age, sex, country, smoke, pa, tg, ldl, hdl, chol, sbp, dbp, bmi, crp, neut, leuk, energy, carb, protein, fat, fiber, diab, crohn, chol_med, bp_med, diab_med, ppi, narrow, broad, site, site_plate, family,family1, shannon, simpson, invsimpson, chao, mgs,cvd2,cvd2.1,cvd2.2)

data <- data[which(complete.cases(data[, c("cacstot", "age", "sex", "country", "site_plate")]) & !cvd), ]
data2 <- data2[complete.cases(data2[, c("cacstot", "age", "sex", "country", "site_plate")]) , ]

#GMM

names(GMM)[1]="scapis_id"

data=merge(data,GMM,by="scapis_id")
data2=merge(data2,GMM,by="scapis_id")

#Export the table

export(data, "processed/data.tsv", na = "NA")
export(data2,"processed/data_CVD.tsv",na="NA")

sessionInfo()



