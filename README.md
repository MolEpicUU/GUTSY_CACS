# GUTSY_CACS

### Code for reproducing the analyses in: " Streptococcus species abundance in the gut is linked to subclinical coronary atherosclerosis in 8,973 participants from the SCAPIS cohort".
***

#### Last update: 2022-June-09

### SCRIPTS (they need to be run in the following order):
* 0_Data/0_create_data.R : It creates the test datasets used in the subsequent analyses. You can also download the datasets directly from this folder.
* 1_simulations/create.random.R: It shuffles the metagenomics and the metadata files to test different models in a random sample maintaining the distributions. In the created file, the association between species and coronary artery calcium score (CACS) follows the null hypothesis.
* 1_simulations/linear_model_robust_random_data.R: It runs the association among the species and the CACS using a robust linear model. CACS and species are ln(x+1) transformed.
* 1_simulations/two_stage_model_random_data.R: It runs the association among the species and the CACS using two-stage model. The first stage is a logistic regression and the second a linear model on the non-zero values in the CACS. Counts on CACS are natural log transformed and species are ln(x+1) transformed.
* 1_simulations/spearman_model_random_data.R: It runs the association among the species and the CACS using Spearman correlation. 
* 1_simulations/nb_model_random_data.R: It runs the association among the species and the CACS using negative binomial model. Species are ln(x+1) transformed.
* 1_simulations/nb_model_clr_random_data.R: It runs the association among the species and the CACS using negative binomial model. Species are center log ratio transformed.
* 1_simulations/nb_hurdle_model_random_data.R: It runs the association among the species and the CACS using a hurdle model. The first part uses a logistic regression and the second part applies a negative binomial model. Species are ln(x+1) transformed.
* 1_simulations/nb_bootsrap_model_random_data.R: It runs the association among the species and the CACS using negative binomial model and bootstrapping the standard errors. Species are ln(x+1) transformed.
* 1_simulations/boot_model_random_data.R: It runs the association among the species and the CACS using linear model and bootstrapping the standard errors. CACS and species are ln(x+1) transformed.
* 1_simulations/boot_resi_random_data.R: It runs the association among the species and the CACS using linear model and bootstrapping the residual values. CACS and species are ln(x+1) transformed.
* 1_simulations/ordinal_model_random_data.R: It runs the association among the species and the CACS using ordinal model. Species are ln(x+1) transformed.
* 1_simulations/linear_model_clr_random_data.R: It runs the association among the species and the CACS using linear model. CACS are ln(x+1) transformed and the species are center log transformed.
* 1_simulations/linear_model_random_data.R: It runs the association among the species and the CACS using linear model. CACS and species are ln(x+1) transformed.
* 1_simulations/qqplot_random_data.R: It runs the lambda estimation for each model and plot the qqplot.
* 2_diversity/linear_model_cacs_shannon.R: It runs the association between Shannon diversity index and CACS adjusting the models for the variables in the basic and full model respectively.
* 2_diversity/permanova_dsbray_curtis_mod1.R: It performs PERMANOVA analysis between the CACS categories and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_mod2.R: It performs PERMANOVA analysis between the CACS categories and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_0_1.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 1) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_0_2.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 2) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_0_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_1_2.R: It performs the pairwise PERMANOVA analysis between the CACS categories (1 and 2) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_1_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (1 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod1_2_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (2 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the basic model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_0_1.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 1) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_0_2.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 2) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_0_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (0 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_1_2.R: It performs the pairwise PERMANOVA analysis between the CACS categories (1 and 2) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_1_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (1 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 2_diversity/permanova_dsbray_curtis_pairwise.mod2_2_3.R: It performs the pairwise PERMANOVA analysis between the CACS categories (2 and 3) and the Bray-Curtis dissimilarity matrix adjusting for the variables in the full model.
* 3_species_level/lm_cacs_CI_mod1.R: It runs the association among the species and CACS using linear model adjusting for the variables in the basic model.
* 3_species_level/lm_cacs_CI_mod2.R: It runs the association among the species and CACS using linear model adjusting for the variables in the full model.
* 4_taxonimc_enrichment/4_taxonomic_enrichment.R: It runs gene-set enrichment analysis on the genus taxonomical level using the ranked p-values from the associations between species and CACS adjusting for the variables in the basic model and stratified by the direction of the regression coefficient.
* 5_sex_strat_species/5_lm_cacs_CI_sex_mod1.R: It runs the associations among the species and CACS stratifying by sex using linear models adjusting for the variables in the basic model.
* 5_sex_strat_species/5_lm_cacs_CI_sex_mod2.R: It runs the associations among the species and CACS stratifying by sex using linear models adjusting for the variables in the full model.
* 5_sex_strat_species/5_heterogeneity_test.R: It runs the Cochran's Q-test to determine if the estimates in the female- and male- stratified analysis are different.
* 6_carotid/ ordinal_carotidplaque_CI_mod1.R: It runs the associations among species and carotid plaques using ordinal regression adjusting for the variables in the basic model.
* 6_carotid/ ordinal_carotidplaque_CI_mod2.R: It runs the associations among species and carotid plaques using ordinal regression adjusting for the variables in the full model.
* 7_validation/7_validation.R: It runs logistic regression between CACS-related species from the basic model and disease status in a case-control study described in Jie, Z. et al. The gut microbiome in atherosclerotic cardiovascular disease. Nat. Commun. 8, 845 (2017).
* 8_GMM_enrichment/8_GMM_enrichment.R: It runs gene-set enrichment analysis on the functional Gut Metabolic Modules (GMM) using the ranked p-values from the associations between species and CACS adjusting for the variables in the basic model and stratified by the direction of the regression coefficient.
* 9_metabolites/metabolites_spearman_mod1.R: It runs partial Spearman correlations among the Streptococcus spp. and the plasma metabolites adjusting for the variables in the basic model.
* 9_metabolites/metabolites_spearman_mod2.R: It runs partial Spearman correlations among the Streptococcus spp. and the plasma metabolites adjusting for the variables in the full model.
* 9_metabolites/subpathway_enrichment.R: It runs gene-set enrichment analysis on the metabolic sub-pathway using the ranked p-values from the associations between Streptococcus spp. and plasma metabolites adjusting for the variables in the full model and stratified by the direction of the regression coefficient.
* 9_metabolites/heterogeniety_test_androgenic_metabolites.R: First, it runs partial Spearman correlations among the Streptococcus spp. and the androgenic plasma metabolites adjusting for the variables in the full model stratifying the analysis by sex. Second, it runs the Cochran's Q-test to determine if the estimates in the female- and male- stratified analysis are different.
* 9_metabolites/linear_model_CACS_metabolites.R: It runs the association among Streptococcus spp. and CACS using linear models adjusting for the variables in the full model and the metabolites involved in the enriched metabolic sub-pathways.
* 10_inflammation_infection/10_infl_infl_lm.R: It runs the association among the Streptococcus spp. and three markers of inflammation and infection. These markers are high-sensitivity C-reactive protein (CRP), neutrophils counts and leukocytes counts.
* 11_gut_oral_microbiome/11_gut_oral_microbiome.R: It runs the association between the Streptococcus spp. in the gut and in the oral cavity using partial Spearman correlation with cluster standard errors (family id as a cluster) adjusting for age, sex, country of birth and metagenomics extraction plate used in the fecal samples.
* 12_oral_health/12_oral_health_full_model.R: It runs the association between the oral Streptococcus spp. and three phenotypes of oral health using ordinal regression with cluster standard errors (family id as a cluster). The models were adjusting for age, sex, smoking, education, oral hygiene, activity realized the hour before attending to the dental examination, and oral Shannon diversity index. The three oral health phenotypes are filled surfaces, gingival inflammation and caries.
* 12_oral_health/12_oral_health_full_bmi_ppi_antib_model.R: It runs the association between the oral Streptococcus spp. and three phenotypes of oral health using ordinal regression with cluster standard errors (family id as a cluster). The models were adjusting for age, sex, smoking, education, oral hygiene, activity realized the hour before attending to the dental examination, oral Shannon diversity index, body mass index, and medication use of proton-pump inhibitors and antibiotics. The three oral health phentoypes are filled surfaces, gingival inflammation and caries.
* functions: Folder that contains functions to be used in the different analyses. These scripts do not need to be run separately. They will be summoned from other scripts.

### Test datasets (all the test datasets are simulated data).
*	0_Data/pheno_2021_06_21.tsv: it contains 100 observations (rows), phenotypic information (columns) and 50 metagenomics species (columns).
*	0_Data/upugut_batch1.mgsComp_tax.txt: it contains 40% of the participants in the pheno_2021_06_21.tsv (columns) file and all the 50 metagenomics species (rows). (6_enrichment_GMM/enrich_GMM.R).
*	0_Data/bray_curtis.tsv: it contains a bray-Curtis dissimilarity matrix estimated from the species in the pheno_2021_06_21.tsv file
*	0_Data/HG3.A.7_tax.txt: It contains the taxonomic information for each specie.
*	0_Data/pheno_validation.tsv: it contains 100 observations (rows) and phenotypic information.
*	0_Data/mgs_validation.tsv: it contains 100 observations (columns) and the relative abundances of the species (rows).
*	0_Data/modules.RData: It contains the information for each MGS to which Gut Metabolic Module belong.
*	0_Data/pheno_2021_09_06_metab.tsv: it contains 100 observations (rows), phenotypic information (columns), 50 metagenomics species (columns), and metabolites (columns).
*	0_Data/metabolic_subpathways.txt: It contains the metabolic sub-pathway information for each metabolite.
*	0_Data/MOS.csv: it contains 100 observations (rows), phenotypic information (columns), 50 gut metagenomics species (columns), and 50 oral metagenomics species (columns).
*	0_Data/pheno_oral.tsv: it contains 100 observations (rows), phenotypic information (columns) and 50 oral metagenomics species (columns).

### Installation
To run the entire code, the R software (https://www.r-project.org/) is needed. The following R packages should be installed: vegan, scales, FRGEpistasis, MASS, boot, rio, car, caret, robust, BiocParallel, sandwich, lmtest, microbiome, pscl, ordinal, reshape, mosaic, ppcor, dplyr, fastDummies, and fgsea.

### Demo
In R: Run the scripts from the first script: "0_Data/0_create_data.R" to the last script "2_oral_health/12_oral_health_full_bmi_ppi_antib_model.R"
