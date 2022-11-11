*****************************************************;
*****************************************************;
*** OVERVIEW OF EXTERNAL PROGRAMS                 ***;
*****************************************************;
*****************************************************;

This file is an organized and combined set of external programs that supported work at the Center for Economic Studies (CES) for 
the environment project. 

The first set of programs run a sequence to consolidate many sources of information about patent technology codes for environment
work, ranging from OECD structures, Popp and Newell (2012), code word searches by research assistants, and internal CES observation 
of clean behavior firms. These program culminate in the env-data-prep1c2.dta file that is provided. This data file was uploaded to 
the Census Bureau and serves several function there as recorded in the CES program set. At the end of this sequence, some simple 
tabulations with citations data of incumbents versus entrants provide the background for the alpha parameter assumptions.

The last program at the end of this consolidated file provides the IPUMS calculation used to model the supply of scientists and 
engineers.

************************************************************************************************;
************************************************************************************************;
*** env-data-prep*.do                                                                        ***;
*** These programs document the sequence used to consolidate our tech code data and model    ***;
*** the alpha parameter. The sequence runs prep1a->1c1,1c2->1b->anal1a.                      ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cap n log close; 
log using env-data-prep1a.log, replace; 

* William Kerr, TGG;
* Describe basic file sets and test merger;
* cd /export/projects/wkerr_h1b_project/kerr/environment/programs;
* Last Modified: July 2012;

clear all; set matsize 5000; set more off;

***********************************;
*** PREPARE BASE ALPHA TEST FILE **;
***********************************;

*** Prepare subclass data;
use ../../recomb/data/USPTO_classes.dta, clear; 
des; sum; sort patent primary;
gen techn=class+"-"+subclass;
egen tech=group(techn); codebook tech;
egen startyr=min(gyear), by(tech);
egen techct=count(primary), by(patent);
keep if primary==1;
keep patent tech techn startyr techct;
sort patent; drop if patent==patent[_n-1];
save ../data/USPTO_classes_env1, replace;

*** Define citation counts;
use patent assignee using ../../mainwork-2009/svdata/p-work-pat_2009, clear;
sort patent assignee; drop if patent==patent[_n-1]; save temp1-pat, replace;
use  ../../mainwork-2009/rawdata/citations, clear;
sort patent; merge patent using temp1-pat; tab _m; drop if _m==2; drop _m; 
for var patent assignee: ren X SX; ren citation patent; 
sort patent; merge patent using temp1-pat; tab _m; drop if _m==2; drop _m;
gen citetot=1; gen citenof=1 if (assignee!=Sassignee | assignee==.);
collapse (sum) cite*, by(patent) fast;
sort patent; save temp1-cites, replace; erase temp1-pat.dta;

*** Prepare patent core;
use ../../mainwork-2009/svdata/p-work-pat_2009, clear;
sort patent; drop if patent==patent[_n-1];
sort patent; merge patent using ../data/USPTO_classes_env1, nok;
tab _m; keep if _m==3; drop _m;
sort patent; merge patent using temp1-cites, nok; erase temp1-cites.dta;
tab _m; table ayear _m; drop _m; for var cite*: replace X=0 if X==.;

*** Save working file;
sort patent; save ../working/env-data-prep1a, replace;
des; sum; 

*** End of program;
log close; 


#delimit;
cap n log close; 
cd /export/projects/wkerr_h1b_project/kerr/environment/programs;
log using env-data-prep1c1.log, replace; 

* William Kerr, TGG;
* do /export/projects/wkerr_h1b_project/kerr/environment/programs/env-data-prep1c1.do;
* Quick reference clean/dirty for patents, part 1 of development, used in prep1b;
* Last Modified: Aug 2012;

clear all; set matsize 5000; set more off;

*****************************************************************;
*** PREPARE AB BROAD TECHNIQUE                                 **;
*****************************************************************;

*** Open working data file;
use patent tech* uspc assignee
    using ../working/env-data-prep1a, clear; des;

*** Prepare rough tech class for merger;
gen temp1=strpos(techn,"-"); replace temp1=temp1+1; gen subclass=substr(techn,temp1,.); drop temp1;
destring subclass, generate(subclass2) ignore("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z") force; 
tab subclass if subclass2==.; drop if subclass2==.; replace subclass2=int(subclass2); compress; gen code=subclass2;
gen class2=uspc; tostring class2 subclass2, replace; 
ren tech ztech; gen tech=class2+"-"+subclass2; drop class2 subclass2;

*** Merge AB main file;  
*** FYI, AB IS THE INITIALS OF THE MAIN RA WORKING ON THE PROJECT ***
*** THIS WORK CODED TECHNOLOGIES AND ALSO PROVIDED CONFIDENCE LEVELS IN ASSIGNMENTS ***
sort tech; merge tech using ./../data/patent-upload-120803b, nok; tab _m; gen ab=(_m==3); drop _m;

*** Consolidate AB codes;
* 3 is most confident;
log off;
gen dirty=0; gen dirtyconf=0;
for any Coal FossF NatGas Oil Petr Shale:
\ replace dirty=1 if cat1=="En-X" \ replace dirtyconf=conf1 if cat1=="En-X" & conf1>dirtyconf & conf1!=.
\ replace dirty=1 if cat2=="En-X" \ replace dirtyconf=conf2 if cat2=="En-X" & conf2>dirtyconf & conf2!=.
\ replace dirty=1 if cat3=="En-X" \ replace dirtyconf=conf3 if cat3=="En-X" & conf3>dirtyconf & conf3!=.
\ replace dirty=1 if cat4=="En-X" \ replace dirtyconf=conf4 if cat4=="En-X" & conf4>dirtyconf & conf4!=.
\ replace dirty=1 if cat5=="En-X" \ replace dirtyconf=conf5 if cat5=="En-X" & conf5>dirtyconf & conf5!=.
\ replace dirty=1 if cat6=="En-X" \ replace dirtyconf=conf5 if cat6=="En-X" & conf6>dirtyconf & conf6!=.;
gen clean=0; gen cleanconf=0;
for any BasicGeo ConUtil GeoTh Nucl Other Solar Wind:
\ replace clean=1 if cat1=="En-X" \ replace cleanconf=conf1 if cat1=="En-X" & conf1>cleanconf & conf1!=.
\ replace clean=1 if cat2=="En-X" \ replace cleanconf=conf2 if cat2=="En-X" & conf2>cleanconf & conf2!=.
\ replace clean=1 if cat3=="En-X" \ replace cleanconf=conf3 if cat3=="En-X" & conf3>cleanconf & conf3!=.
\ replace clean=1 if cat4=="En-X" \ replace cleanconf=conf4 if cat4=="En-X" & conf4>cleanconf & conf4!=.
\ replace clean=1 if cat5=="En-X" \ replace cleanconf=conf5 if cat5=="En-X" & conf5>cleanconf & conf5!=.
\ replace clean=1 if cat6=="En-X" \ replace cleanconf=conf5 if cat6=="En-X" & conf6>cleanconf & conf6!=.;
gen pol=0; gen polconf=0;
for any Air Other SWaste Water:
\ replace pol=1 if cat1=="Pol-X" \ replace polconf=conf1 if cat1=="Pol-X" & conf1>polconf & conf1!=.
\ replace pol=1 if cat2=="Pol-X" \ replace polconf=conf2 if cat2=="Pol-X" & conf2>polconf & conf2!=.
\ replace pol=1 if cat3=="Pol-X" \ replace polconf=conf3 if cat3=="Pol-X" & conf3>polconf & conf3!=.
\ replace pol=1 if cat4=="Pol-X" \ replace polconf=conf4 if cat4=="Pol-X" & conf4>polconf & conf4!=.
\ replace pol=1 if cat5=="Pol-X" \ replace polconf=conf5 if cat5=="Pol-X" & conf5>polconf & conf5!=.
\ replace pol=1 if cat6=="Pol-X" \ replace polconf=conf5 if cat6=="Pol-X" & conf6>polconf & conf6!=.;
log on;
tab cat1 if dirty+clean+pol==0;
for any BasicGeo ConUtil Other: tab uspc if cat1=="En-X";

*** Perform AB mega assignments and corrections;
* 405 a bit iffy;
for any 166 175 196 405 507: 
\ sum dirty clean pol if uspc==X
\ replace dirty=1 if dirty==0 & uspc==X
\ replace dirtyconf=2 if dirtyconf==0 & uspc==X
\ replace ab=1 if uspc==X;
for any 208: 
\ sum dirty clean pol if uspc==X & (code!=12 & (code<24 | code>38)) 
\ replace dirty=1 if dirty==0 & uspc==X & (code!=12 & (code<24 | code>38)) 
\ replace dirtyconf=2 if dirtyconf==0 & uspc==X & (code!=12 & (code<24 | code>38))
\ replace ab=1 if uspc==X & (code!=12 & (code<24 | code>45));
for any 376 976: 
\ sum dirty clean pol if uspc==X
\ replace clean=1 if clean==0 & uspc==X
\ replace cleanconf=2 if cleanconf==0 & uspc==X
\ replace ab=1 if uspc==X;
for any 588: 
\ sum dirty clean pol if uspc==X
\ replace pol=1 if pol==0 & uspc==X
\ replace polconf=2 if polconf==0 & uspc==X
\ replace ab=1 if uspc==X;
for any 250 60 323 44 D13: cap n sum dirty clean pol if uspc==X;

*** Summarize current point;
cap n assert dirty+clean+pol>=1 if ab==1;
table dirty clean;
table dirty clean if pol==1;
tab cat1;
tab cat1 if dirty==1;
tab cat1 if clean==1;
tab cat1 if pol==1;
for var dirty clean pol: egen TX=sum(X), by(assignee);
sort assignee;
for var dirty clean pol: sum T* if TX>=10;
keep patent assignee techn uspc ab dirty dirtyconf clean cleanconf pol polconf;
for var dirty dirtyconf clean cleanconf pol polconf: ren X AB1X; ren ab AB1;
sum;

*****************************************************************;
*** MERGE WITH AB MANUAL SEARCHES                              **;
*****************************************************************;

*** Merge in AB manual searches;
ren techn tech; sort tech; merge tech using ../data/patent-upload-120803x, keep(type cer) nok; 
tab _m; gen AB2=(_m==3); drop _m; ren tech techn; for var type cer: ren X AB2X;

*****************************************************************;
*** MERGE WITH POPP GROUPINGS                                  **;
*****************************************************************;

*** Merge in Popp data;
sort techn; merge techn using ../data/patent-popp-120803, keep(alt clean dirty) nok;
tab _m; gen POPP=(_m==3); drop _m; for var alt clean dirty: ren X PX;

*****************************************************************;
*** SAVE WORKING FILE TO ADJUST LATER                          **;
*****************************************************************;

*** Save final file;
compress; sort patent; 
save ../working/env-data-prep1c1, replace; 
des; sum;

*** End of program;
log close; 

#delimit;
cap n log close; 
cd /export/projects/wkerr_h1b_project/kerr/environment/programs;
log using env-data-prep1c2.log, replace; 

* William Kerr, TGG;
* do /export/projects/wkerr_h1b_project/kerr/environment/programs/env-data-prep1c2.do;
* Quick reference clean/dirty for patents, part 2 of development, used in prep1b;
* Last Modified: Aug 2012;

clear all; set matsize 5000; set more off;

*****************************************************************;
*** PREPARE PRIORITY OF AB AND POPP GROUPINGS - AGGRESSIVE     **;
*****************************************************************;

*** Use working file;
use ../working/env-data-prep1c1, clear; des; sum;
keep if AB1==1 | AB2==1 | POPP==1; sum;

*** Priority 1: Popp;
for any alt clean dirty: ren PX X;
replace alt=. if alt==0;
su alt clean dirty;

*** Priority 2: AB Manual Searches;
for var alt clean dirty: tab AB2type if X==1;
drop if AB2type=="n" & POPP==0;
replace clean=1 if AB2type=="c" & clean==. & dirty==.;
replace dirty=1 if AB2type=="d" & clean==. & dirty==.;
replace dirty=1 if AB2type=="p" & clean==. & dirty==.;
drop AB2*;
su alt clean dirty;

*** Priority 3: AB Broad Search;
replace clean=1 if AB1clean==1 & clean==. & dirty==.;
replace dirty=1 if AB1dirty==1 & clean==. & dirty==.;
replace dirty=1 if AB1pol==1 & clean==. & dirty==.;
drop AB1*;
su alt clean dirty;

*** Narrow and save temp file;
keep if (alt==1 | clean==1 | dirty==1);
compress; sort patent; 
save ../working/env-data-prep1c2, replace; des; sum;

*** List features;
use ../working/env-data-prep1c2, clear; gen ct=clean; collapse (sum) ct, by(techn) fast; sum, d; gsort -ct; list in 1/50;
use ../working/env-data-prep1c2, clear; gen ct=dirty; collapse (sum) ct, by(techn) fast; sum, d; gsort -ct; list in 1/50;
use ../working/env-data-prep1c2, clear; gen ct=alt; collapse (sum) ct, by(techn) fast; sum, d; gsort -ct; list in 1/50;

*** End of program;
log close; 

#delimit;
cap n log close; 
log using env-data-prep1b.log, replace; 

* William Kerr, TGG;
* Describe basic file sets and test merger;
* cd /export/projects/wkerr_h1b_project/kerr/environment/programs;
* Last Modified: Aug 2012;

clear all; set matsize 5000; set more off;

***********************************;
*** PREPARE LEAD PATENT DESIGN   **;
***********************************;

*** Merge US citations data together;
use  ../../mainwork-2009/rawdata/citations, clear;
sort patent; merge patent using ../working/env-data-prep1a, keep(assignee ayear ctryn tech uspc scat); 
tab _m; keep if _m==3; keep if ctryn=="US"; keep if (ayear>=1975 & ayear<=2004); drop ctryn _m; 
for var patent tech uspc scat assignee ayear: ren X SX; ren citation patent;
sort patent; merge patent using ../working/env-data-prep1a, keep(assignee ayear ctryn tech uspc scat); 
tab _m; keep if _m==3; keep if ctryn=="US"; keep if (ayear>=1975 & ayear<=2004); drop ctryn _m; 
for var patent tech uspc scat assignee ayear: ren X DX;

*** Collapse citations to patent level by year of citation;
gen ucitetot=1; 
gen ucitenof1=1 if (Dassignee!=Sassignee | Sassignee==. | Dassignee==.);
gen ucitenof2=1 if (Dassignee!=Sassignee | Sassignee==. | Dassignee==.) & (Suspc==Duspc & Suspc!=. & Duspc!=.);
collapse (sum) ucite* (mean) Duspc Dayear, by(Dpatent Sayear) fast;
for var ucite*: replace X=0 if X==.; sum; cap n assert Sayear>=Dayear;

*** Merge in patent designation;
ren Dpatent patent; sort patent; merge patent using ../working/env-data-prep1c2, nok; tab _m; drop _m; ren patent Dpatent;

*** Define entrants and incumbents;
gen entry=(Sayear-Dayear)>=0 & (Sayear-Dayear)<=3; gen incumb=(Sayear-Dayear)>=5 & (Sayear-Dayear)<=10;
tab entry; table Sayear Dayear if entry==1; table Sayear Dayear if entry==1 & clean==1; table Sayear Dayear if entry==1 & dirty==1;
tab incumb; table Sayear Dayear if incumb==1; table Sayear Dayear if incumb==1 & clean==1; table Sayear Dayear if incumb==1 & dirty==1;

*** Define indicators for high entry citations - general;
for var ucite*:
\ gen temp0=X if incumb==1 \ egen temp1=pctile(temp0), by(Duspc Sayear) p(90) \ sum temp1 \ count if temp1==. \ gen T90X=X>temp1 & X>1 & X!=. \ drop temp*;
for var ucitenof1:
\ gen temp0=X if incumb==1 \ egen temp1=pctile(temp0), by(Duspc Sayear) p(95) \ sum temp1 \ count if temp1==. \ gen T95X=X>temp1 & X>1 & X!=. \ drop temp*
\ gen temp0=X if incumb==1 \ egen temp1=pctile(temp0), by(Duspc Sayear) p(99) \ sum temp1 \ count if temp1==. \ gen T99X=X>temp1 & X>1 & X!=. \ drop temp*
\ gen temp0=X if incumb==1 \ egen IX=count(temp0), by(Duspc Sayear) \ drop temp*;
for var T*: gen EX=X if entry==1; drop T90* T95* T99*; sum;

*** Define indicators for high entry citations - energy, 90th (removed Duspc category to allow overlap);
for var ucitenof1: 
\ gen temp0B=X if incumb==1 & (clean==1 | dirty==1) \ egen temp1B=pctile(temp0B), by(Sayear) p(90) \ sum temp1B \ count if temp1B==.
\ gen temp0C=X if incumb==1 & clean==1 \ egen temp1C=pctile(temp0C), by(Sayear) p(90) \ sum temp1C \ count if temp1C==.
\ gen temp0D=X if incumb==1 & dirty==1 \ egen temp1D=pctile(temp0D), by(Sayear) p(90) \ sum temp1D \ count if temp1D==.
\ gen TB90X=X>temp1B & X>1 & X!=. \ gen ETB90X=TB90X if entry==1 \ egen IBX=count(temp0B), by(Sayear)
\ gen TC90X=X>temp1C & X>1 & X!=. \ gen ETC90X=TC90X if entry==1 \ egen ICX=count(temp0C), by(Sayear)
\ gen TD90X=X>temp1D & X>1 & X!=. \ gen ETD90X=TD90X if entry==1 \ egen IDX=count(temp0D), by(Sayear)
\ drop temp*;

*** Collapse data;
collapse (max) E* I* alt clean dirty, by(Dpatent) fast;
ren Dpatent patent; sort patent; save temp1-lead, replace; sum;

*** Merge back into base data;
use if (ayear>=1985 & ayear<=1999) & ctryn=="US"
    using ../working/env-data-prep1a, clear;
sort patent; merge patent using temp1-lead; erase temp1-lead.dta; 
tab _m; drop if _m==2; drop _m;
for var E*: replace X=0 if X==.;
for var ETB* ETC* ETD* IC* ID*: replace X=. if (clean!=1 & dirty!=1); sum;
sort patent; save ../working/env-data-prep1b, replace;

*** End of program;
log close; 

#delimit;
cap n log close; 
log using env-data-anal1a.log, replace; 

* William Kerr, TGG;
* Describe basic file sets and test merger;
* cd /export/projects/wkerr_h1b_project/kerr/environment/programs;
* Last Modified: Aug 2012;

clear all; set matsize 5000; set more off;

***********************************;
*** ANALYSIS - FULL              **;
***********************************;

*** Open working data file;
use patent E* I* ayear scat startyr tech* cite* uspc alt dirty clean
    using ../working/env-data-prep1b, clear; des;
sort uspc; merge uspc using ../data/uspc-codes; tab _m; drop if _m==2; drop _m;

*** Analyze overall probabilities;
sum E*; 
sum ET90ucitenof1 if citenof>0;

*** Analyze energy probabilities;
keep if clean==1 | dirty==1 | alt==1;
sum E*; 
sum ET90ucitenof1 if citenof>0;
sum E* if clean==1;
sum ET90ucitenof1 if clean==1 & IC>10 & IC!=.;
sum E* if dirty==1;
sum ET90ucitenof1 if dirty==1 & ID>10 & ID!=.;
sum E* if alt==1; 
sum E* if alt==1 | (uspc==376 | uspc==976); 

tab uspcn if clean==1;
tab uspcn if clean==1 & ETC90==1;
tab uspcn if clean==1 & ETD90==1;
tab uspcn if dirty==1;
tab uspcn if dirty==1 & ETC90==1;
tab uspcn if dirty==1 & ETD90==1;
tab uspcn if alt==1;

*** End of program;
log close;

************************************************************************************************;
************************************************************************************************;
*** c-env01.do                                                                               ***;
*** Performs simple tabulations on IPUMS 5% Samples for Energy-Related Industries            ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cd; cd unix/census/programs/environment;
* do /export/home/faculty/wkerr/unix/census/programs/environment/c-env01.do;
cap n log close; log using c-env01.log, replace;

* William Kerr, TGG;
* 2000 Ed Shares in Env;
* Last Updated Aug 2012;

clear all; set matsize 1000; set more off;

!gunzip ../../data/2000-5sp5.dta;
use if (age>=20 & age <=65) & gqtyped==0 & (occ!=0 & occ!=.) & educ99!=.
    using ../../data/2000-5sp5.dta, clear;

* Consolidate education codes;
gen ed="1-nohs" if educ<=9; 
replace ed="2-hs" if educ==10;
replace ed="3-somecol" if educ==11 | educ==12;
replace ed="4-col" if educ==14;
replace ed="5-adv" if educ>14 & educ!=.;

* Consolidate occupations;
gen se=0;
replace se=1 if (occ>=132 & occ<=176);
replace se=2 if (occ>=100 & occ<=131);
replace se=3 if (occ>=190 & occ<=196);

* Define energy;
for var ind*: tab X;
gen energy=0;
for num 37/38 207 308 318 57/59: replace energy=1 if ind==X;
keep if energy==1;

* Total;
tab ed [fw=perwt];
tab se;
tab se [fw=perwt];

* Same with educated sE;
drop se; gen se=0;
replace se=1 if (occ>=132 & occ<=176) & educ>=14;
replace se=3 if (occ>=190 & occ<=196) & educ>=14;
tab se;
tab se [fw=perwt];

save energy-extract.dta;
cap n erase energy-extract.dta.gz;
!gzip energy-extract.dta;

*** END OF PROGRAM;
log close;