# load libraries

library(rio)

# import data

pheno <- import("raw/MODS CM Sergi.sav")
oral <- import("raw/malmos02_r1v1_040_MGS-abundance.xlsx")
gut <- import("raw/MOS_clin_mgs_n2262_210324.csv")
drugs <- import("raw/MOS_Drugvariables_PPI_AB_211013.sav")
diet <- import("raw/MOS_fiber_energy_220425.csv")
family <- import("raw/MOSlopnr_family_n2644.csv")
family_bmi <- import("raw/BMI_family_MODS.csv")

# clean

oral_info <- oral[, 1:11]
rownames(oral) <- oral$MGS
oral <- t(oral[, 12:ncol(oral)])
rownames(oral) <- pheno$lopnrMOS[match(rownames(oral), pheno$Sample_names)]

# match

id <- pheno$lopnrMOS[which(pheno$lopnrMOS %in% rownames(oral))]
pheno <- pheno[match(id, pheno$lopnrMOS), ]
oral <- oral[match(id, rownames(oral)), ]
gut <- gut[match(id, gut$ID_num), ]
drugs <- drugs[match(id, drugs$lopnrMOS), ]
diet <- diet[match(id, diet$lopnr_mos), ]
family <- family[match(id, family$lopnrMOS), ]
family_bmi <- family_bmi[match(id, family_bmi$lopnrMOS), ]

# phenotypes

lopnr <- pheno$lopnrMOS
caries <- pheno$DiS32_DS32
fillings <- pheno$FS_32
gingivitis <- pheno$BoP
sex <- ifelse(pheno$Sex_k == 1, "Female", "Male")
age <- pheno$Age
plate <- gut$Plate.x
country <- ifelse(gut$country_birth == 1, "Swedish", "Non-Swedish")
smoke <- pheno$Smoking_ctr
smoke[which(pheno$Smoking_ctr == 1)] <- "Never smoker"
smoke[which(pheno$Smoking_ctr == 2)] <- "Current smoker"
smoke[which(pheno$Smoking_ctr == 4)] <- "Former smoker"
education <- pheno$Education_ctr
education[which(pheno$Education_ctr == 2)] <- "Completed primary school"
education[which(pheno$Education_ctr == 3)] <- "Qualified vocational training"
education[which(pheno$Education_ctr == 4)] <- "University"
antibiotics <- ifelse(drugs$Antibiotics_All == 1, "Yes", "No")
ppi <- ifelse(drugs$PPI_All == 1, "Yes", "No")
lasthour <- pheno$LastHour
lasthour[which(pheno$LastHour == 9)] <- NA
eat <- ifelse(lasthour %in% c(1, 12, 123), "Yes", "No")
smoke2 <- ifelse(lasthour %in% c(2, 12, 23, 123), "Yes", "No")
brush <- ifelse(lasthour %in% c(13, 23, 123), "Yes", "No")
lasthour <- ifelse(eat == "Yes" | smoke2 == "Yes" | brush == "Yes", "Yes", "No")
plaque <- pheno$Plaque_score
bmi <- gut$BMI
bmi[which(is.na(bmi))] <- family_bmi$BMI[which(is.na(bmi))]
family <- family$family
family[which(is.na(family))] <- family_bmi$family[which(is.na(family))]

# mgs

gut <- apply(gut[, grepl("HG3A", colnames(gut))], 2, as.numeric)
colnames(gut) <- sapply(strsplit(colnames(gut), split = "[.]"), function(x) paste(x[1], x[2], sep = "."))

# make and export tables

data <- data.frame(lopnr, caries, fillings, gingivitis, sex, age, plate, country, smoke, education, antibiotics, ppi, eat, smoke2, brush, lasthour, plaque, gut, bmi, family, gut, oral)
data <- data[which(!is.na(age) & !is.na(sex)), ]

export(data, "processed/mods.tsv", na = "NA")
export(oral_info, "processed/mods_tax.tsv")

sessionInfo()
