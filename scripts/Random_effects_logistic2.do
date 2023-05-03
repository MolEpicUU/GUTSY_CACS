/*Code for running logistic mixed models for athero, stenosis.
These models are run for significant species (n=64)
*/

import delimited "/proj/nobackup/sens2019512/Projects/cacs_mgsv2/results/main_significant.tsv",clear //creating a local (`VList') with names of the sig species
replace id=lower(id)
replace id=subinstr(id,".","",.)
drop if id=="hg3a1967" | id=="hg3a1800" | id=="hg3a0270"
local VList ""
count
forv i=1/`=r(N)'{
	local test=id[`=`i'']
	local VList `VList' `test'
}

import delimited "/proj/nobackup/sens2019512/Projects/cacs_mgsv2/processed/data.tsv",clear //dataset with all species
keep `VList' sis stenosis athero plaque cacstot_cat age carb protein fiber sbp dbp ///
		chol hdl ldl tg bmi sex country smoke diab diab_med pa chol_med ///
		bp_med site_plate* family

*Pre-processing*********************
foreach var of varlist tg-dbp carb-fiber athero-plaque pa smoke diab diab_med chol_med bp_med{ 
	replace `var'="" if `var'=="NA"
	replace `var'=subinstr(`var',"DUKE_","",.)
	replace `var'=subinstr(`var',"SIS_","",.)
}
replace plaque="0" if plaque=="None"
replace plaque="1" if plaque=="Unilateral"
replace plaque="2" if plaque=="Bilateral"
replace stenosis="1" if stenosis=="Yes"
replace stenosis="0" if stenosis=="No"
replace athero="1" if athero=="Yes"
replace athero="0" if athero=="No"
replace cacstot_cat="1" if cacstot_cat=="1-100"
replace cacstot_cat="2" if cacstot_cat=="101-400"
replace cacstot_cat="3" if cacstot_cat==">400"
destring *,replace

foreach var of varlist sex-pa diab chol_med-family{
	encode `var',gen(`var'_num)
	drop `var'
	rename `var'_num `var'
}
xtset family

foreach var of varlist hg*{
	replace `var'=log(`var'+1)
	egen `var't=std(`var')
	drop `var'
	rename `var't `var'
}
****************************

*************Running the models
foreach v2 of varlist stenosis athero{
	cap postclose myfile
	postfile myfile str20 Variable est_mixed pval_mixed lci_mixed hci_mixed using "/home/ulfha881/Desktop/Ulf//`v2'.dta",replace

	*Main models (stored in "stenosis_MGS.dta" and "athero_MGS.dta")
	foreach var of varlist hg*{
		cap xtlogit `v2' `var' age carb protein fiber sbp dbp ///
		chol hdl ldl tg bmi i.sex i.country i.smoke i.diab i.diab_med i.pa i.chol_med ///
		i.bp_med i.site_plate*,or
		if _rc==0{
			matrix A=r(table)
			post myfile ("`var'") (A[1,1]) (A[4,1]) (A[5,1]) (A[6,1])
		}
		if _rc!=0 post myfile ("`var'") (.) (.) (.) (.)
	}
	postclose myfile
}
