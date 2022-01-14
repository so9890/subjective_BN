clear

global maindir "C:/Users/glv2/Dropbox/Dropbox/HSV/taxation/QJE_Final/Data/Section_8"


//SORT DATASETS AND HARMONIZE COUNTRY NAMES

//World Tax Database
use ${maindir}/WTI.dta
replace country = "Macao" if country == "Macao Special Administrative Region of China"
replace country = "Iran" if country == "Iran, Islamic Republic of"
replace country = "Lao" if country == "Lao People's Democratic Republic"
replace country = "Lybia" if country == "Libyan Arab Jamahiriya"
replace country = "Syria" if country == "Syrian Arab Republic"
replace country = "St. Vincent and the Grenadines" if country == "Saint Vincent and the Grenadines"
replace country = "USSR" if country == "Union of Soviet Socialist Republics [former]"
replace country = "Tanzania" if country == "United Republic of Tanzania"
drop if country == "Serbia and Montenegro" //Keep Serbia and Montenegro separately

gen prog1 = 1-(1-mr_y/100)/(1-ar_y/100)
gen prog2 = 1-(1-mr_2y/100)/(1-ar_2y/100)
gen prog3 = 1-(1-mr_3y/100)/(1-ar_3y/100)
gen prog4 = 1-(1-mr_4y/100)/(1-ar_4y/100)
gen prog_avg = (prog1+prog2+prog3+prog4)/4
drop mr_y mr_2y mr_3y mr_4y mr_y ar_y ar_2y ar_3y ar_4y prog1 prog2 prog3 prog4

sort country year
save ${maindir}/WTI_NEW.dta, replace

//Standardized World Income Inequality Database
use ${maindir}/SWIID.dta
replace country = "Kyrgyzstan" if country == "Kyrgyz Republic"
replace country = "Saint Lucia" if country == "St. Lucia"
replace country = "Yemen" if country == "Yemen, Republic of"
replace country = "Turks and Caicos Islands" if country == "Turks and Caicos"
drop if country == "Serbia and Montenegro" //Keep Serbia and Montenegro separately
sort country year
save ${maindir}/SWIID_NEW.dta, replace


//World Income Database
use ${maindir}/WID.dta
sort country year
save ${maindir}/WID.dta, replace

//Penn World Tables
use ${maindir}/PWT.dta
replace country = "Macao" if country == "China: Macao SAR"
replace country = "Iran" if country == "Iran (Islamic Republic of)"
replace country = "Lao" if country == "Lao People's Democratic Republic"
replace country = "Syria" if country == "Syrian Arab Republic"
replace country = "China" if country == "China, People's Republic of"
replace country = "Hong Kong" if country == "China: Hong Kong SAR"
replace country = "Moldova" if country == "Republic of Moldova"
replace country = "Slovakia" if country == "Slovak Republic"
replace country = "Macedonia" if country == "The Former Yugoslav Republic of Macedonia"
replace country = "Tanzania" if country == "United Republic of Tanzania: Mainland"
sort country year
save ${maindir}/PWT_NEW.dta, replace

// MERGE WTI AND GINI DATA
use ${maindir}/WTI_NEW.dta, clear
sort country year
save ${maindir}/datamerge0.dta, replace

merge country year using ${maindir}/SWIID_NEW.dta
sort country year
rename _merge _merge_gini
drop avg_gini_market_string
replace avg_gini_market = avg_gini_market/100 
save ${maindir}/datamerge1.dta, replace

// MERGE ALSO WITH TOPINCOME
merge country year using ${maindir}/WID.dta
rename _merge _merge_top
gen paretocoeff2 = paretocoeff^2
sort country year
save ${maindir}/datamerge2.dta, replace

//MERGE ALSO WITH PWT
merge country year using ${maindir}/PWT_NEW.dta
rename _merge _merge_pwt
gen popint = ceil((pop))
sort country year
save ${maindir}/datamerge3.dta, replace


//FINAL DATASET
use ${maindir}/datamerge3.dta
sort country year

**SELECTION;
drop if prog_avg == . | mrp_all == .
drop if _merge_pwt ~=3
drop if _merge_gini ~=3
keep if year>=1990

** MERGE WITH DATA ON VARIANCE OF LOG-NORMAL PIECE;
merge country year using ${maindir}/DATA_SIGMA.dta
rename _merge _merge_jon
save ${maindir}/datamerge4.dta, replace
sort country year

**DROP DUPLICATES;
duplicates report country year
duplicates tag country  year, gen(isdup)
quietly by country year:  gen dup = cond(_N==1,0,_n)
drop if dup>1
		
**DESTRING REGION AND COUNTRY;
encode region, generate(regionnum)
encode country, generate(countrynum)
gen wb_classifnum = .
replace wb_classifnum =1 if wb_classif=="Low"
replace wb_classifnum =2 if wb_classif=="Lower middle"
replace wb_classifnum =3 if wb_classif=="Upper middle"
replace wb_classifnum =4 if wb_classif=="High"

**WEIGHTS;
gen sqrtgdp   = (gdp)^0.5
gen sqrtpop   = (pop)^0.5

** COUNT OBSERVATIONS
by country: gen countcountry = _n
by country: gen totcountry = _N

**REGRESSIONS AND TABLES
eststo clear
eststo: reg mrp_all g_y_pwt avg_gini_market i.year [pweight=sqrtgdp]
eststo: reg mrp_all g_y_pwt paretocoeff paretocoeff2 sigma i.year   [pweight=sqrtgdp]
eststo: reg mrp_all g_y_pwt avg_gini_market i.year i.regionnum [pweight=sqrtgdp]
eststo: reg mrp_all g_y_pwt paretocoeff paretocoeff2 sigma i.year i.regionnum [pweight=sqrtgdp]
eststo: reg mrp_all g_y_pwt avg_gini_market i.year i.regionnum  i.wb_classifnum [pweight=sqrtgdp]
eststo: reg mrp_all g_y_pwt paretocoeff paretocoeff2 sigma i.year i.regionnum i.wb_classifnum [pweight=sqrtgdp]
esttab, ar2 drop(*_cons*  *year* *region* *classif*) compress

eststo clear
eststo: reg prog_avg g_y_pwt avg_gini_market i.year [pweight=sqrtgdp]
eststo: reg prog_avg g_y_pwt paretocoeff paretocoeff2 sigma i.year   [pweight=sqrtgdp]
eststo: reg prog_avg g_y_pwt avg_gini_market i.year i.regionnum [pweight=sqrtgdp]
eststo: reg prog_avg g_y_pwt paretocoeff paretocoeff2 sigma i.year i.regionnum [pweight=sqrtgdp]
eststo: reg prog_avg g_y_pwt avg_gini_market i.year i.regionnum  i.wb_classifnum [pweight=sqrtgdp]
eststo: reg prog_avg g_y_pwt paretocoeff paretocoeff2 sigma i.year i.regionnum i.wb_classifnum [pweight=sqrtgdp]
esttab, ar2 drop(*_cons*  *year* *region* *classif*) compress

