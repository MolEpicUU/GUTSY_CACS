//combine the results from the other files to a single Excel-file
sysdir set PLUS "/home/ulfha881/Desktop/plus"
sysdir set PERSONAL "/home/ulfha881/Desktop/plus"

foreach var in sis stenosis plaque athero duke cacstot_catfull{
	use "/home/ulfha881/Desktop/Ulf//`var'.dta",clear
	noi di "`var'"
	count
	local v2 "`var'"
	if "`var'"=="cacstot_catfull" local v2 "cacstot"
	export excel "/home/ulfha881/Desktop/Ulf/Results.xlsx",firstrow(variable) sheet("`v2'",replace)
}
