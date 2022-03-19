*** Barrage (2019) "Optimal Dynamic Carbon Taxes in a Climate-Economy Model with Distortionary Fiscal Policy" ***
******************************************************************************************************************
*Author: Lint Barrage (lint.barrage@gmail.com)

*****************************************************************
***	Damage Function Input: Malaria Labor Productivity Impacts ***
*****************************************************************

*This do-file combines World Bank data on malaria prevalence and populations with Tol's (2008) 
*estimates of climate-induced changes in malaria incidence in order to compute regional TFP losses
*from reduced adult labor productivity associated with increases in malarious childhoods

*Input 1: "WB_Malaria_Pop_raw.csv", World Bank malaria prevalence, population, and GDP (World Development Indicators, https://datacatalog.worldbank.org/dataset/world-development-indicators)
*Input 2:  Tol (2008) estimates of climate change-induced increases in malaria morbidity by region (entered manually below)
		  *Tol, Richard "Climate, Development, and Malaria: An Application of FUND" Climatic Change (2008), 88:21-34.
*Output:  "Malaria_TFP.xlsx" - regional and aggregate global estimates of TFP losses associated with 2.5C warming

clear
set more off
global dir = "C:\Users\lintb\Desktop\COMET\Supplementary\COMET Calibration\COMET Damages"

*************************************
***Step 1: Process World Bank Data***
*************************************
insheet using "$dir\WB_Malaria_Pop_raw.csv", clear names

rename v5 y2008
keep countryname indicatorname y2008 /*Use 2008 baseline to match Tol (2008)*/
replace indicator = "malaria" if indicator=="Notified cases of malaria (per 100,000 people)"
replace indicator = "population" if indicator=="Population, total"
replace indicator = "gdp" if indicator=="GDP (constant 2000 US$)"
keep if indicator=="population" | indicator=="malaria" | indicator=="gdp"

reshape wide y2008, i(country) j(indicator) string

foreach var of varlist y2008* {
	local temp = subinstr("`var'","y2008","",.)
	rename `var' `temp'
	}
	*
drop if country=="World"


************************************
***Step 2: Add Tol (2008) Results***
************************************

***Step 2.1: Match countries to FUND regions
gen region = ""
replace region = "USA" if country=="United States"
replace region = "CAN" if country=="Canada"
replace region = "WEU" if (country=="Andorra"|country=="Austria"|country=="Belgium"|country=="Cyprus"|country=="Denmark"|country=="Finland"|country=="France"|country=="Germany"|country=="Greece"|country=="Iceland"|country=="Ireland"|country=="Italy"|country=="Liechtenstein"|country=="Luxembourg"|country=="Malta"|country=="Monaco"|country=="Netherlands"|country=="Norway"|country=="Portugal"|country=="San Marino"|country=="Spain"|country=="Sweden"|country=="Switzerland"|country=="United Kingdom")
replace region = "JPK" if (country=="Japan"|country=="Korea, Rep.")
replace region = "ANZ" if (country=="Australia"|country=="New Zealand")
replace region = "CEE" if (country=="Albania"|country=="Bosnia and Herzegovina"|country=="Bulgaria"|country=="Croatia"|country=="Czech Republic"|country=="Hungary"|country=="Macedonia, FYR"|country=="Poland"|country=="Romania"|country=="Slovak Republic"|country=="Slovenia"|country=="Kosovo"|country=="Montenegro"|countr=="Serbia")
replace region = "FSU" if (country=="Armenia"|country=="Azerbaijan"|country=="Belarus"|country=="Estonia"|country=="Georgia"|country=="Kazakhstan"|country=="Latvia"|country=="Lithuania"|country=="Moldova"|country=="Russian Federation"|country=="Tajikistan"|country=="Turkmenistan"|country=="Ukraine"|country=="Uzbekistan"|country=="Kyrgyz Republic")
replace region = "MDE" if (country=="Bahrain"|country=="Iran, Islamic Rep."|countr=="Iraq"|country=="Israel"|country=="Jordan"|country=="Kuwait"|country=="Lebanon"|countr=="Oman"|country=="Qatar"|country=="Saudi Arabia"|country=="Syrian Arab Republic"|countr=="Turkey"|country=="United Arab Emirates"|country=="West Bank and Gaza"|country=="Yemen, Rep.")
replace region = "CAM" if (countr=="Belize"|country=="Costa Rica"|countr=="El Salvador"|country=="Guatemala"|country=="Honduras"|country=="Mexico"|countr=="Nicaragua"|countr=="Panama")
replace region = "SAM" if (country=="Argentina"|country=="Bolivia"|country=="Brazil"|country=="Chile"|country=="French Guiana"|country=="Guyana"|country=="Paraguay"|country=="Peru"|country=="Suriname"|country=="Uruguay"|country=="Venezuela, RB"|country=="Colombia"|country=="Ecuador")
replace region = "SAS" if (country=="Afghanistan"|country=="Bangladesh"|country=="Bhutan"|country=="India"|country=="Nepal"|countr=="Pakistan"|country=="Sri Lanka")
replace region = "SEA" if (countr=="Brunei Darussalam"|countr=="Cambodia"|countr=="Timor-Leste"|country=="Indonesia"|country=="Lao PDR"|country=="Malaysia"|country=="Myanmar"|country=="Papua New Guinea"|country=="Philippines"|country=="Singapore"|countr=="Taiwan"|country=="Thailand"|country=="Vietnam")
replace region = "CHI" if (country=="China"|country=="Hong Kong SAR, China"|country=="Korea, Dem. Rep."|country=="Macao SAR, China"|country=="Mongolia")
replace region = "NAF" if (country=="Algeria"|country=="Egypt, Arab Rep."|country=="Libya"|country=="Morocco"|country=="Tunisia"|country=="Western Sahara")
replace region = "SSA" if (country=="Angola"|country=="Benin"|country=="Botswana"|country=="Burkina Faso"|country=="Burundi"|country=="Cameroon"|country=="Cape Verde"|country=="Central African Republic"|country=="Chad"|country=="Congo, Dem. Rep."|country=="Congo, Rep."|country=="Cote d'Ivoire"|country=="Djibouti"|country=="Equatorial Guinea"|country=="Eritrea"|country=="Ethiopia"|country=="Gabon"|country=="Gambia, The"|country=="Ghana"|country=="Guinea"|country=="Guinea-Bissau"|country=="Kenya"|country=="Lesotho"|country=="Liberia"|country=="Madagascar"|country=="Malawi"|country=="Mauritania"|country=="Mozambique"|countr=="Namibia"|country=="Niger"|countr=="Nigeria"|country=="Rwanda"|country=="Senegal"|country=="Sierra Leone"|country=="Somalia"|country=="South Africa"|country=="Sudan"|countr=="Swaziland"|country=="Tanzania"|countr=="Togo"|countr=="Uganda"|country=="Zambia"|country=="Zimbabwe"|country=="Mali"|country=="South Sudan")
replace region = "SIS" if missing(region)
*check:
tab country if region=="SIS" /*looks good*/

***Step 2.2: Enter Tol (2008) baseline Malaria morbidity rates (in in years of life diseased per million people)
gen cmorb = 0
replace cmorb = 350 if region=="MDE"
replace cmorb = 140 if region=="CAM"
replace cmorb = 140 if region=="SAM"
replace cmorb = 500 if region=="SAS"
replace cmorb = 370 if region=="SEA"
replace cmorb = 40 if region=="CHI"
replace cmorb = 350 if region=="NAF"
replace cmorb = 5300 if region=="SSA"
replace cmorb = 140 if region=="SIS"

***Step 2.3: Enter Tol (2008) climate-induced Malaria change estimates (in thousands of life years per degree C)
gen change = 0
replace change = 5.2 if region=="MDE"
replace change = 1.3 if region=="CAM"
replace change = 3.4 if region=="SAM"
replace change = 45 if region=="SAS"
replace change = 13.1 if region=="SEA"
replace change = 4.1 if region=="CHI"
replace change = 3.2 if region=="NAF"
replace change = 218.5 if region=="SSA"
replace change = 0.5 if region=="SIS"

*******************************************
***Step 3: Calculate Productivity Losses***
*******************************************

*Convert country population to millions:
replace population = population/1000000
label var population "Pop in millions"

*Calculate country total baseline malaria incidence in life-years
gen total_cmorb = cmorb*population
label var total_cmorb "Baseline malaria incidence in life-years"

*Aggregate to regional total malaria icidence
bysort region: egen total_cmorb_region = total(total_cmorb)

*Calculate percent change in malaria incidence due to climate change:
replace change = change*1000
gen pct_mal_change = change/total_cmorb_region
replace pct_mal_change = 0 if missing(pct_mal_change)

*Calculate baseline malaria population incidence 
gen pop2 = population*10
label var pop2 "Population in 100,000s"
gen total_malaria_cases = malaria*pop2
bysort region: egen total_malaria_cases_region = total(total_malaria_cases)
gen pop3 = population*1000000
label var pop3 "Population"
bysort region: egen pop3_region = total(pop3)
gen malaria_fraction_affected = total_malaria_cases_region/pop3_region

*Calculate the fraction of the population additionally sick per degree C
gen share_affected = malaria_fraction_affected*pct_mal_change

*** Main Calculation***
***********************
*Bleakely (2003) estimates that a malarious childhood decreases adult wages by 15%
*We first translate this labor productivity loss into the equivalent TFP loss based on the COMET Cobb-Douglas specification
*The marginal product of labor for a household with malaria is 0.85 of a healthy household's, implying:
*(z)^(1-alpha-v) = .85
local alpha = 0.3
local v = 0.03
local z = (.85)^(1/(1-`alpha'-`v'))
local z_loss = 1-`z'

*Change in z due to 2.5C warming:
gen temp1 = (`z_loss')*share_affected*(2.5)
sum temp1

*Corresponding residual TFP factor:
gen temp2 = (1-temp1)^(1-`alpha'-`v')

*Corresponding loss in TFP:
gen TFP_loss_region = 1-temp2

*Convert to percentage:
replace TFP_loss_region = TFP_loss_region*100
tabulate region, summarize(TFP_loss_region)

*Calculate global aggregate GDP-weighted impact
bysort region: egen regdp = total(gdp)
egen global_gdp = total(gdp)
gen regshare = regdp/global
bysort region: gen m = _n
gen weighted_loss = TFP_loss_region*regshare if m==1
egen global_loss = total(weighted_loss)

************************************************************
***Step 4: Aggregate to DICE Region-Level via GDP Weights***
************************************************************

*Step 4.1: Map into DICE regions
gen DICE_region = ""
replace DICE = "OHI" if region=="ANZ"
replace DICE = "LA" if region=="CAM"
replace DICE = "OHI" if region=="CAN"
replace DICE = "EE/FSU" if region=="CEE"
replace DICE = "China" if region=="CHI"
replace DICE = "Russia" if region=="FSU"
replace DICE = "Japan" if region=="JPK"
replace DICE = "MidEast" if region=="MDE" | region=="NAF"
replace DICE = "LA" if region=="SAM" | region=="SIS"
replace DICE = "India" if region=="SAS"
replace DICE = "OthAsia" if region=="SEA"
replace DICE = "SSA" if region=="SSA"
replace DICE = "US" if region=="USA"
replace DICE = "WE/Euro" if region=="WEU"
count if mi(DICE)

*Step 4.2: Compute GDP shares for multi-FUND to single-DICE matches
bysort DICE: egen DICE_gdp = total(gdp)
bysort region: egen reg_gdp = total(gdp)
gen gdp_share = reg_gdp/DICE_gdp

*Step 4.3: Compute GDP-share-weighted malaria TFP losses by DICE region
*To avoid confusion, keep only regions
keep DICE* region gdp_share reg_gdp TFP_loss_region global_loss
duplicates drop
*br /*looks good*/
gen temp = gdp_share*TFP_loss_region 
bysort DICE_r: egen TFP_loss_DICE_region = total(temp)
*br /*looks good*/


************
*** SAVE ***
************

keep DICE_r TFP_loss_DICE global_loss
duplicates drop
export excel "$dir\Malaria_TFP.xlsx", replace firstrow(varlabels)


