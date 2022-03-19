******************************************************************************************************************
*** Barrage (2019) "Optimal Dynamic Carbon Taxes in a Climate-Economy Model with Distortionary Fiscal Policy" ***
******************************************************************************************************************
*Author: Lint Barrage (lint.barrage@gmail.com)

*************************************************
***	TABLE 2: Government Expenditure Shares***
*************************************************

*This do-file processes IMF Government Finance Statistics and IMF International Financial Statistics data in order to produce paper Table 2
*Input 1: IMF_GFS_raw.csv, IMF Government Finance Statistics, "Expense" Database (see: https://www.imf.org/en/Data, pulled for all available countries, accessed July 2012)
*Input 2: IMF_IFS_raw.csv, IMF International Financial Statistics, "Government Finance" Indicators, "Expense" variables (see: https://www.imf.org/en/Data, pulled for all available countries, accessed July 2012)
*Input 3: WB_GDP_2005_raw.csv, World Bank $2005 PPP GDP Data (World Development Indicators, https://datacatalog.worldbank.org/dataset/world-development-indicators, accessed July 2012)
*Output: Table2.xlsx

clear
set more off
global dir = "C:\Users\lintb\Desktop\COMET\Supplementary\COMET Calibration\COMET Government Expenditures_Table2"

*********************************
***Step 0: Prepare WB GDP Data***
*********************************
insheet using "$dir\WB_GDP_2005_raw.csv", names clear
rename countryname country
rename v5 gdp
replace gdp = gdp/1000000000
label var gdp "GDP in 2005$ billions"

**Adjust country names to match IMF data:
replace country = "Bahrain, Kingdom of" if country=="Bahrain"
replace country = "China, P.R.: Hong Kong" if country=="Hong Kong SAR, China"
replace country = "China, P.R.: Macao" if country=="Macao SAR, China"
replace country = "Congo, Republic of" if country=="Congo, Rep."
replace country = "Egypt" if country=="Egypt, Arab Rep."
replace country = "Iran, Islamic Republic of" if country=="Iran, Islamic Rep."
replace country = "Yemen, Republic of" if country=="Yemen, Rep."

keep country gdp
sort country
save "$dir\gdp2005.dta", replace


*******************************
***Step 1: Insheet IMF Data****
*******************************
insheet using "$dir\IMF_GFS_raw.csv", clear
save "$dir\temp_imf_raw1.dta", replace
insheet using "$dir\IMF_IFS_raw.csv", clear
append using "$dir\temp_imf_raw1.dta"
drop datasource status frequency
duplicates drop

*****************************************
***Step 2: Clean and Process IMF Data****
*****************************************
**Calibration base year: 2005
keep if time==2005

*Drop sub-categories so as to avoid double-counting (see IMF GFS Manual 2001, Table 6.1):
drop if concept=="Social Assistance Benefits, Cash (Gen. Govt.)"			/*Sub-category of Social Benefits*/
drop if concept=="Social Security Benefits, Cash (Gen. Govt.)"				/*Sub-category of Social Benefits*/
drop if concept=="Employer Social Benefits, Cash (Gen. Govt.)"				/*Sub-category of Social Benefits*/
drop if concept=="To Non-residents, Foreign Interest, Cash (Gen. Govt.)"	/*Sub-category of Interest*/
drop if concept=="Property Expense Other than Interest, Cash (Gen. Govt.)"	/*Sub-category of Interest*/
drop if concept=="Miscellaneous Other Expense, Cash (Gen. Govt.)"			/*Sub-category of Other expense*/

*Drop unnecessary variables:
drop if concept=="Gross Domestic Product, Real"
drop if concept=="Gross Domestic Product, Deflator"

*Some countries have duplicates up to slight rounding error across the two databases; eliminate to avoid double-counting:
bysort concept country time unit: egen max = max(value)	
replace value = max
duplicates drop

*Generate accounting basis indicator:
gen test = word(concept, -3)
gen cash = (test=="Cash")
gen noncash = (test=="Noncash")
drop test

*Adjust currency units into billions:
replace value = value/1000000000 if unit~="Percent of GDP"
label var value "billions"
replace value = value/100 if unit=="Percent of GDP"

*Prepare denominator variables for share calculations:
egen num = group(country)
gen temp = value if concept=="Gross Domestic Product, Nominal"
bysort num: egen gdp_nom = max(temp)
drop temp

*Compute shares:
gen gdp_share = value/gdp_nom if unit~="Percent of GDP"
replace gdp_share = 0 if value==0
replace gdp_share = value if unit=="Percent of GDP"	

*Due to rounding there are some slight differences between recorded and calculated shares -> Harmonize:
bysort concept country time: egen max2 = max(gdp_share)
replace gdp_share = max2
drop max max2 unit value
duplicates drop

*Drop redundant totals:
drop if concept=="Total Expenditure, Cash (Gen. Govt.)"
drop if concept=="Expense, Cash (Gen. Govt.)"
drop if concept=="Total Expenditure, Noncash (Gen. Govt.)"
drop if concept=="Expense, Noncash (Gen. Govt.)"
drop if concept=="Gross Domestic Product, Nominal"
drop if concept=="Total expenditure, Cash (Gen. Govt.)"

*Check number of records per country:
gen counter = 1 
bysort country time: egen zahl = total(counter)
tab zahl	/*Looks good - up to 15 is normal (7 cash + 8 noncash indicators)*/
drop counter zahl

*Eliminate consumption of fixed capital:
sum gdp_share if concept=="Consumption of Fixed Capital, Noncash (Gen. Govt.)" /*Mean is 1.86% with std. dev 1.06% and min (0.00642%) and max (4.9%) -> Small*/
drop if concept =="Consumption of Fixed Capital, Noncash (Gen. Govt.)"

*Simplify "Concept" descriptions:
split concept, p(", ")
replace concept = concept1
drop concept1 concept2

*For countries with both cash and noncash basis information, keep the former:
gen counter = 1
bysort country time concept: egen test = total(counter)
drop if test==2 & noncash==1

*Clean:
drop counter test gdp_nom
duplicates drop
drop if country=="Euro Area"

*** Merge in World Bank PPP GDP Data***
***************************************
sort country
merge country using "$dir\gdp2005.dta", uniqusing
tab _merge

*Check GDP share of countries covered by IMF data:
replace _merge=999 if country=="World"
bysort _merge: egen gdp_total = total(gdp)
tab gdp_total if _merge==2	/*Miss 16.882 trillion dollars out of 56.4 trillion -> Account for 71% of world GDP*/

*Calculate Data-World GDP Shares
sum gdp if country=="World"
local world = r(mean)
sum gdp_total if _merge==2
local data_world = r(mean)
local value = `world'-`data_world'
gen data_world_gdp_share = gdp/`value'

*Calculate Weighted Shares*
***************************
gen w_gdp_share = gdp_share*data_world_gdp_share

***********
**TOTALS***
***********
drop if _merge==2	/*Countries not in IMF data*/
drop if _merge==999 /*World aggregate*/

bysort concept time: egen weighted_global_gdp_share = total(w_gdp_share)
tabulate concept, summarize(weighted_global_gdp_share)

keep concept weighted_global_gdp_share
duplicates drop

*Subtract Interest from totals (endogenous in model):
drop if concept=="Interest"

*Aggregate into Government Consumption vs. Transfers:
gen concept2 = "Government Consumption"
replace concept2 = "Social Benefits" if concept=="Social Benefits"
replace concept2 = "" if concept=="Interest"
bysort concept2: egen weighted_global_gdp_share2 = total(weighted_global_gdp_share) if concept!="Interest"

*Compute Total:
egen total_gdp_share = total(weighted_global_gdp_share) 
sum total_gdp_share

*Compute expenditure shares:
gen gov_spend_share = (weighted_global_gdp_share2/total_gdp_share)


*********************************
*** Step 3: Make PAPER TABLE 2***
*********************************

keep concept2 weighted_global_gdp_share2 gov_spend_share total_gdp_share
duplicates drop

rename weighted pct_GDP

expand 2 in 1
replace concept2 = "Total" if _n==3
sum total
local temp = r(mean)
replace pct_GDP = `temp' if concept2=="Total"
replace gov_spend_share = . if concept2=="Total"
drop total

export excel "$dir\Table2.xlsx", replace

**********************************
**********************************
**********************************


