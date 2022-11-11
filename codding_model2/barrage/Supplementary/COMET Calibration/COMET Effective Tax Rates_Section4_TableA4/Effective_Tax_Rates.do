*****************************************************************************************************************
*** Barrage (2019) "Optimal Dynamic Carbon Taxes in a Climate-Economy Model with Distortionary Fiscal Policy" ***
******************************************************************************************************************
*Author: Lint Barrage (lint.barrage@gmail.com)

****************************************
***	Baseline Effective Tax Rates ***
****************************************

*This do-file processes IMF Government Finance Statistics and IMF International Financial Statistics data in order to produce paper Table 2
*Input 1: IMF_GFS_raw.csv, IMF Government Finance Statistics, "Expense" Database (see: https://www.imf.org/en/Data, pulled for all available countries, accessed July 2012)
*Input 2: IMF_IFS_raw.csv, IMF International Financial Statistics, "Government Finance" Indicators, "Expense" variables (see: https://www.imf.org/en/Data, pulled for all available countries, accessed July 2012)
*Input 3: WB_GDP_2005_raw.csv, World Bank $2005 PPP GDP Data (World Development Indicators, https://datacatalog.worldbank.org/dataset/world-development-indicators, accessed July 2012)
*Output: Table2.xlsx

clear
set more off
global dir = "C:\Users\lintb\Desktop\COMET\Supplementary\COMET Calibration\COMET Effective Tax Rates_Section4_TableA4"

*********************************
***Step 0: Prepare WB GDP Data***
*********************************

insheet using "$dir\WB_GDP_2005_raw.csv", names clear
rename countryname country
rename v5 gdp

keep country gdp
sort country
save "$dir\gdp2005_ETR.dta", replace


***********************************
***Step 1: Process ETR Estimates***
***********************************

import excel using "$dir\Effective_Tax_Rates_raw.xlsx", first clear

*Basic cleaning:
replace tl = trim(tl)
replace tl = subinstr(tl,"…","",.)
destring tl, replace
rename Country country
bysort country: gen n = _n
count if n==1
drop if mi(country)

*Adjust country names to match World Bank GDP data:
replace country = trim(country)
replace country = "Bosnia and Herzegovina" if country=="Bosnia"
replace country = "Dominican Republic" if country=="Dominican R."
replace country = "Egypt, Arab Rep." if country=="Egypt"
replace country = "Hong Kong SAR, China" if country=="Hong Kong"
replace country = "Iran, Islamic Rep." if country=="Iran"
replace country = "Korea, Rep." if country=="Korea"
replace country = "Macedonia, FYR" if country=="Macedonia"
replace country = "Russian Federation" if country=="Russia"
drop if country=="Taiwan"
replace country="Trinidad and Tobago" if country=="Trinidad"
replace country = "Venezuela, RB" if country=="Venezuela"
replace country = "Slovak Republic" if country=="Slovakia"

*Merge on World Bank GDP data:
sort country
merge country using "$dir\gdp2005_ETR.dta", nokeep
tab _merge	/*looks good*/


*********************************
***Step 2: Calculate Averages ***
*********************************

***Country-Level Averages -> Online Appendix Table A4
*************************
bysort country: egen avg_tl = mean(tl)
*br country tl avg_tl
bysort country: egen avg_tk = mean(tk)
*br country tk avg_tk
bysort country: egen avg_tc = mean(tc)
*br country tc avg_tc
bysort country: egen avg_tcl = mean(t_c_l)

preserve
keep country avg_tl avg_tk avg_tc 
duplicates drop
export excel using "$dir\TableA4.xlsx", firstrow(var) replace
restore


***Global GDP-Weighted Averages -> Paper Section 4.5
*******************************
keep country avg_* gdp
duplicates drop

*Calculate GDP weights by tax rate*
gen y_tl = gdp if !missing(avg_tl)
gen y_tk = gdp if !missing(avg_tk)
gen y_tc = gdp if !missing(avg_tc)
gen y_tcl = gdp if !missing(avg_tcl)

egen total_y_tl = total(y_tl)
egen total_y_tk = total(y_tk)
egen total_y_tc = total(y_tc)
egen total_y_tcl = total(y_tcl)

gen share_tl = gdp/total_y_tl
gen share_tk = gdp/total_y_tk
gen share_tc = gdp/total_y_tc
gen share_tcl = gdp/total_y_tcl

*** Counts:
count if !missing(avg_tl)
count if !missing(avg_tk)
count if !missing(avg_tc)
count if !missing(avg_tcl)

*Weighted averages:
gen wavg_tl = share_tl*avg_tl
gen wavg_tk = share_tk*avg_tk
gen wavg_tc = share_tc*avg_tc
gen wavg_tcl = share_tcl*avg_tcl

*Global GDP-Weighted averages:
egen tao_l = total(wavg_tl)
egen tao_k = total(wavg_tk)
egen tao_c = total(wavg_tc)
egen tao_cl = total(wavg_tcl)
gen tao_cl_2 = ((tao_l/100)+(1-(tao_l/100))*(tao_c/100))*100

*Summary statistics:
sum tao_l tao_k tao_c tao_cl tao_cl_2

*Save:
keep tao_l tao_k tao_c tao_cl_2
duplicates drop

export excel using "$dir\Effective_Tax_Rates_Global.xlsx", firstrow(var) replace

