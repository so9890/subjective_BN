*****************************************************;
*****************************************************;
*** OVERVIEW OF CENTER OF ECONOMIC STUDY PROGRAMS ***;
*****************************************************;
*****************************************************;

This file is an organized and combined set of programs used at the Center for Economic Studies for the environment project. 

The programs have been redacted according to CES guidelines to remove 1) file directory root information and 2) any 
firm-specific information in the programs. In several places, this second part regarding firm-specific work is substantial.
There is also one program (env-firm1-anal5-questions.do) that is references below but not disclosed in any form as it is solely
devoted to processing appropriately the firm-level records and correctly designating the firms to be included in the sample. Full 
programs are available in the BR0598 project archive for those with appropriate security clearance and data access points.

The data required to execute these programs include Longitudinal Business Database, the R&D Survey, and patent data. Patent data 
are connected to Census firm through firm-name matching algorithms. The algorithms and final match dataset are available through 
the Census Bureau by requesting the BR0598 patent match. We separately provide our clean/dirty technology code concordance for 
environment-related patents.

************************************************************************************************;
************************************************************************************************;
*** env-firm1-anal5-regr1.do                                                                 ***;
*** This file generates all of our reported regressions and moments. It also serves to       ***;
*** prepare the final sample that was used for Census Bureau disclosure testing. The BR598   ***;
*** folders contain archives of earlier variants of this program and extensions.             ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cd @ROOT/programs/reall/files/energy1;
cap n log close; log using env-firm1-anal5-regr1.log, replace;

* William Kerr, TGG;
* Environment Regressions;
* Last Modified: Oct 2012;

clear all; *set matsize 11000; *set maxvar 15000; set more off;

***********************************************************;
*** Prepare citation weighted patents using uspc        ***;
***********************************************************;

*** Unique patent merge file; *** BELOW DRAWS UPON RAW PATENT DATA FILES ***
use patent assignee uspc ayear using @ROOT/trips/data-patents/p-work-pat_2008_CES, clear;
sort patent; drop if patent==patent[_n-1]; save temp1-pat1, replace;

*** Merge assignees and categories onto patents;
use @ROOT/ghosh/cluster/1_PATENT/citations, clear;
des; sum; cap n assert cit_appyear<=cit_year if cit_appyear!=.;
sort patent; merge patent using temp1-pat1, keep(assignee uspc ayear) nok; 
tab _m; for var _m patent assignee uspc ayear: ren X SX; ren citation patent;
sort patent; merge patent using temp1-pat1, keep(assignee uspc ayear) nok; erase temp1-pat1.dta;
tab _m; for var _m patent assignee uspc ayear: ren X DX;

*** Restrict and collapse patent citation counts;
keep if (Dayear>=1975 & Dayear<=2005); keep if D_m==3;
gen cite_tot=1; gen cite_nof=1 if (Sassignee!=Dassignee | Sassignee==.);
collapse (sum) cite_* (mean) Duspc Dayear, by(Dpatent) fast;
ren Dpatent patent; for var cite_*: egen temp1=mean(X), by(Duspc Dayear) \ gen nX=X/temp1 \ drop temp1 \ replace X=0 if X==.; format ncite* %5.2f;
sort patent; save env-firm1-anal5-prep1-cite1, replace; sum; for var cite* ncite*: tab Dayear, s(X);

***********************************************************;
*** Prepare group of patent tech for env firm sample    ***;
***********************************************************;

*** Prepare new patent cut - flow and stock;
use if (assignee!=. & ctryn=="US") using @ROOT/trips/data-patents/p-work-pat_2008_CES, clear;
gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (ayear>=X-2 & ayear<=X+2); 
sort patent; merge patent using env-data-prep1c2; tab _m; drop if _m==2; gen ab=(_m==3); drop _m;
sort patent; merge patent using env-firm1-anal5-prep1-cite1, nok keep(cite* ncite*); tab _m; drop _m; for var cite* ncite*: replace X=0 if X==.;

*** Merge in assignee-firm; *** THIS IS THE PATENT-CENSUS BRIDGE ***
sort assignee; merge assignee using @ROOT/trips/data-patents/BRDC01-MI-comb-best, nok; 
keep if _m==3; gen ct=1; 

*** FIRM SPECIFIC PIECES EXCLUDED ***

*** Generate stock values 1992;
for any 1992:
\ egen STctX=sum(ayear>=1975 & ayear<=X+2), by(assignee)
\ egen STabX=sum(ayear>=1975 & ayear<=X+2 & ab==1), by(assignee)
\ egen STclX=sum(ayear>=1975 & ayear<=X+2 & clean==1), by(assignee)
\ egen STdiX=sum(ayear>=1975 & ayear<=X+2 & dirty==1), by(assignee);

*** Generate stock values 2002;
for any 2002:
\ egen STctX=sum(ayear>=1975 & ayear<=X+2), by(assignee)
\ egen STabX=sum(ayear>=1975 & ayear<=X+2 & ab==1), by(assignee)
\ egen STclX=sum(ayear>=1975 & ayear<=X+2 & clean==1), by(assignee)
\ egen STdiX=sum(ayear>=1975 & ayear<=X+2 & dirty==1), by(assignee);

*** Generate cite-weighted stock values 2002;
for any 2002:
\ gen temp1=ncite_nof if ayear>=1975 & ayear<=2004
\ egen STWctX=sum(temp1), by(assignee)
\ egen STWabX=sum(temp1*ab), by(assignee)
\ egen STWclX=sum(temp1*clean), by(assignee)
\ egen STWdiX=sum(temp1*dirty), by(assignee);
sum ST*; tab ayear if yr==.; drop if yr==.;

*** Collapse data;
collapse (sum) ct ab alt clean dirty POPP cite* ncite* (mean) ST*, by(assignee firm yr) fast;
collapse (sum) ct ab alt clean dirty POPP cite* ncite* ST*, by(firm yr) fast;
for var ct ab alt clean dirty POPP cite* ncite* ST*: replace X=0 if X==.; ren ct pat4;
sort firm yr; compress; save env-firm1-anal5-base-pat, replace; sum;
collapse (sum) ab, by(firm); keep if ab>2;
keep firm; duplicates drop; sort firm; save temp1-firm-pat, replace;

***********************************************************;
*** Build Operating and R&D Dataset                     ***;
***********************************************************;

*** Build long R&D file; *** BELOW DRAWS UPON R&D SURVEY APPENDED ACROSS YEARS ***
use @ROOT/programs/reall/working/rad-pat/unbpanel, clear;
keep id year dns dne set rdtot rdfed; egen temp1=rsum(dns dne set rdtot); keep if temp1>0 & temp1!=.; drop temp1;
for num 2/5:
\ append using @ROOT/ghosh/data/rd3/rd200X
\ egen temp1=rsum(dns dne set rdtot)
\ keep if temp1>0 & temp1!=.
\ keep id year dns dne set rdtot rdfed;
for var dns dne rdtot rdfed set: replace X=. if X==0 \ tab year, s(X); replace dns=. if dns==99999999;
ren id firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; duplicates drop; drop if firm==""; 
collapse (sum) dns dne rdtot rdfed set, by(firm year) fast; 
gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2); tab year if yr==.; drop if yr==.;
collapse (mean) dns dne rdtot rdfed set, by(firm yr) fast; 
egen temp1=sum(rdtot), by(firm); drop if (temp1==0 | temp1==.); drop temp1;
compress; sort firm yr; save temp1-firm, replace; 
append using temp1-firm-pat; erase temp1-firm-pat.dta; 
keep firm; duplicates drop; sort firm; save temp1-firm-list, replace;

*** Build long LBD file; *** BELOW DRAWS UPON LBD FILES APPENDED OVER YEARS ***
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd1976a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m; compress; save temp1-lbd, replace;
forvalues YR1=1977(1)2001 {;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd`YR1'a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd.dta; compress; save temp1-lbd, replace;
};
use lbdnum yr emp cfn firmid bestsic if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd2002a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
replace firm=firmid; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; ren best sic; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd.dta; compress; drop if firm==""; save temp1-lbd, replace;
ren emp LB_emp; tab yr, s(LB_emp);

*** Collapse and merge R&D;
gen sic4=substr(sic,1,4); drop sic; destring sic4, replace force; gen LB_est=1;
collapse (sum) LB_emp LB_est, by(firm sic4 yr) fast; 
gen sic3=int(sic4/10); replace sic3=999 if sic3==.; gen I_tot4=1;
collapse (sum) LB_emp LB_est I_tot4, by(firm sic3 yr) fast; gen I_tot=1;
collapse (sum) LB_emp LB_est I_tot4 I_tot, by(firm yr) fast; 
ren yr year; gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2);
collapse (mean) LB_emp LB_est I_tot I_tot4, by(firm yr) fast; 
sort firm yr; merge firm yr using temp1-firm; tab _m; drop _m;
for var rdtot set: replace X=0 if X==.;
pwcorr LB_emp dne; sum LB_emp dne;

*** Save base file and clean-up;
compress; sort firm yr; save env-firm1-anal5-base-ces, replace; sum; sum, d;
for any temp1-firm temp1-firm-list temp1-lbd: erase X.dta;

***********************************************************;
*** Prepare group of patent tech for env firm sample    ***;
***********************************************************;

*** Open CES base and define year range;
use env-firm1-anal5-base-ces, clear; 
for var LB_emp LB_est I_tot I_tot4 rdtot rdfed dns dne: replace X=int(X);
format LB_emp LB_est I_tot I_tot4 rdtot rdfed dns dne %8.0f; sum;

*** Baseline work on variables;
for var dne LB_emp: replace X=1.5e6 if X>1.5e6 & X!=.;
ren dne zdne; egen dne=rmean(zdne LB_emp); 
replace dne=zdne if LB_emp==. & dne==.;
replace dne=LB_emp if zdne==. & dne==.;
drop if dne==.; pwcorr dne zdne LB_emp;
replace I_tot=1 if (I_tot==. | I_tot==0) & dne>0 & dne!=.; 
replace I_tot4=1 if (I_tot4==. | I_tot4==0) & dne>0 & dne!=.; sum; 

*** Merge in patent data to core;
sort firm yr; merge firm yr using env-firm1-anal5-base-pat; tab _m; drop if _m==2; drop _m;
for var pat4 alt clean dirty POPP cite* ncite*: replace X=0 if X==.; 

*** Impute small R&D values where set and pat are known;
sum pat4 set rdtot if rdtot>0 & rdtot<1000;
sum pat4 set rdtot if rdtot==0; 
gen rdproxy=264.1/4.9*set if rdtot==0 & set>0; 
replace rdproxy=264.1/1.5*pat4 if rdproxy==. & rdtot==0 & pat4>0;
replace rdproxy=1000 if rdproxy>1000 & rdproxy!=.;
replace rdproxy=. if rdtot!=0;
replace rdtot=rdproxy if rdtot==0 & rdproxy!=.;

*** Apply new conditional rule of entry with R&D;
egen yrmin=min(yr), by(firm); 
gen flag=(yr==yrmin & rdtot==0); tab flag; tab flag if yr==yrmin;
egen zdrop=max(flag), by(firm); count; drop if zdrop==1; count;
drop yrmin zdrop flag;

***********************************************************;
*** Prepare Variables for Estimations                   ***;
***********************************************************;

*** Define environmental sample;
egen temp1=sum(ab), by(firm); egen temp2=sum(pat4), by(firm); gen absh=temp1/temp2; 
gen temp3=(temp1>2 & absh>0.10 & absh!=.); egen envsamp=max(temp3), by(firm); 
gen temp4=(temp1>15 & absh>0.25 & absh!=.); egen envsamp2=max(temp4), by(firm); drop temp*;
egen temp1=sum(POPP), by(firm); egen temp2=sum(pat4), by(firm); gen poppsh=temp1/temp2; 
gen temp3=(temp1>2 & poppsh>0.10 & absh!=.); egen envsamp3=max(temp3), by(firm); drop temp*;

*** Define group types;
gen temp1=(dirty/(dirty+clean)>0.75); egen envcoal=max(temp1), by(firm yr); egen temp2=mean(temp1), by(firm); gen envcoalx=(temp2>0.5); drop temp*;
gen temp1=(alt/(dirty+clean)>0.15); egen envalt=max(temp1), by(firm yr); egen temp2=mean(temp1), by(firm); gen envaltx=(temp2>0.5); drop temp*;
for var envcoal* envalt*: replace X=. if envsamp==0; 
for var rdtot rdfed: table envcoalx yr if rdproxy==., c(sum X) row col f(%8.0f);
gen rdfedsh=rdfed/rdtot; for num 0 1: sum rdfed rdfedsh if envcoalx==X, d;
tab firm if rdfed>500000 & rdfed!=. & envcoalx==0;
tab firm if rdfed>500000 & rdfed!=. & envcoalx==1;

*** Variable prep;
egen patmin=min(pat4), by(firm); sum patmin, d;
gen rddns=rdtot/dns; gen rdemp=rdtot/dne; gen rdict=rdtot/I_tot; gen rd4ict=rdtot/I_tot4; gen rdest=rdtot/LB_est;
gen setdns=set/dns; gen setemp=set/dne; gen setict=set/I_tot; gen set4ict=set/I_tot4; gen setest=set/LB_est;
gen patdns=pat4/dns; gen patemp=pat4/dne; gen patict=pat4/I_tot; gen pat4ict=pat4/I_tot4; gen patest=pat4/LB_est;
gen ncitedns=ncite_nof/dns; gen nciteemp=ncite_nof/dne; gen nciteict=ncite_nof/I_tot; gen ncite4ict=ncite_nof/I_tot4; gen nciteest=ncite_nof/LB_est;
for var set* dns dne LB_emp LB_est I_tot I_tot4 pat* rd* ncite*:
\ sort firm yr \ gen lX=ln(X) \ gen DlX=lX-lX[_n-1] if firm==firm[_n-1] & yr==yr[_n-1]+5;
for var rdtot rddns rdemp rdict setemp: sort firm yr \ gen LlX=ln(X[_n-1]) if firm==firm[_n-1];
for var lrdtot lLB_emp LB_emp: egen mn_X=mean(X), by(firm); 

*** Final details;
egen cyr=group(envcoal yr); egen ayr=group(envalt yr);
egen fgroup=group(firm); gen wt1=1; gen wt2=mn_lLB_emp; set seed 1;

***********************************************************;
*** Show Distributions and Eta Elasticity Calculations  ***;
***********************************************************;

*** OLS FD - Full;
for any emp: areg DlpatX DlrdX [aw=wt1], a(cyr) cl(fgroup);

*** Restrict panel;
keep if envsamp==1;
do env-firm1-anal5-questions1.do;  *** THIS IS THE EXCLUDED FIRM-SPECIFIC PREP PROGRAM ***
count; codebook firm; egen yrct=count(yr),by(firm);
egen temp1=count(yr), by(firm); sort firm; tab temp1 if firm!=firm[_n-1];
tab temp1 if firm!=firm[_n-1] & patmin>0; drop temp1;

*** Descriptive tabulations by size;
tab envcoal; tab envcoalx;
gen temp1=(dirty/(dirty+clean)); sum temp1, d; 
count; count if (temp1>0.25 & temp1<0.75); sum temp1 if envcoal==0, d; drop temp1;

********************;
*** Min Yr Count ***;
********************;

preserve;
gen temp1=(rdproxy!=. & rdproxy!=0); egen temp2=max(temp1), by(firm); drop if temp2>0; drop temp*;
keep if yrct>=3;

*** OLS XSECT - Env Sample;
for any rd set: 
\ areg lnciteict lXict [aw=wt1], a(yr) cl(fgroup)
\ areg lncite4ict lX4ict [aw=wt1], a(yr) cl(fgroup)
\ areg lncitedns lXdns [aw=wt1], a(yr) cl(fgroup)
\ areg lnciteest lXest [aw=wt1], a(yr) cl(fgroup);
for any rd set: 
\ areg lpat4ict lX4ict [aw=wt1], a(yr) cl(fgroup)
\ areg lncite4ict lX4ict [aw=wt1] if absh>0.3, a(yr) cl(fgroup)
\ areg lncite4ict lX4ict [aw=wt1] if envsamp3==1, a(yr) cl(fgroup);

*** OLS FD - Env Sample;
for any rd set: 
\ areg Dlnciteict DlXict [aw=wt1], a(yr) cl(fgroup)
\ areg Dlncite4ict DlX4ict [aw=wt1], a(yr) cl(fgroup)
\ areg Dlncitedns DlXdns [aw=wt1], a(yr) cl(fgroup)
\ areg Dlnciteest DlXest [aw=wt1], a(yr) cl(fgroup);
for any rd set: 
\ areg Dlpat4ict DlX4ict [aw=wt1], a(yr) cl(fgroup)
\ areg Dlncite4ict DlX4ict [aw=wt1] if absh>0.3, a(yr) cl(fgroup)
\ areg Dlncite4ict DlX4ict [aw=wt1] if envsamp3==1, a(yr) cl(fgroup);

*** FE Poisson Model Prep;;
set seed 1; egen fgrp=group(firm); xtset fgrp yr; 
cons de 1 lrdtot+lI_tot=1; cons de 2 lset+lI_tot=1; 
cons de 3 lrdtot+lI_tot4=1; cons de 4 lset+lI_tot4=1;
cons de 5 lrdtot+lLB_est=1; cons de 6 lset+lLB_est=1;

*** FE Poisson Models by Type;
xi i.yr;
tab envcoalx;
xtpoisson ncite_nof lrdtot lI_tot4 _Iy*, re const(3) iter(50) vce(boot);
xtpoisson ncite_nof lrdtot lLB_est _Iy*, re const(5) iter(50) vce(boot);
xtpoisson ncite_nof lrdtot lI_tot4 _Iy*, fe const(3) iter(50) vce(boot);
xtpoisson ncite_nof lrdtot lLB_est _Iy*, fe const(5) iter(50) vce(boot);

xtpoisson ncite_nof lset lI_tot4 _Iy*, re const(4) iter(50) vce(boot);
xtpoisson ncite_nof lset lLB_est _Iy*, re const(6) iter(50) vce(boot);
xtpoisson ncite_nof lset lI_tot4 _Iy*, fe const(4) iter(50) vce(boot);
xtpoisson ncite_nof lset lLB_est _Iy*, fe const(6) iter(50) vce(boot);

xtpoisson pat4 lrdtot lI_tot4 _Iy*, fe const(3) iter(50) vce(boot);
xtpoisson pat4 lset lI_tot4 _Iy*, fe const(4) iter(50) vce(boot);
xtpoisson ncite_nof lrdtot lI_tot4 _Iy*, fe iter(50) vce(boot);
xtpoisson ncite_nof lset lI_tot4 _Iy*, fe iter(50) vce(boot);

restore;

********************;
*** Full Yr Count **;
********************;

preserve;
gen temp1=(rdproxy!=. & rdproxy!=0); egen temp2=max(temp1), by(firm); drop if temp2>0; drop temp*;

*** OLS XSECT - Env Sample;
for any rd set: 
\ areg lncite4ict lX4ict [aw=wt1], a(yr) cl(fgroup);

*** OLS FD - Env Sample;
for any rd set: 
\ areg Dlncite4ict DlX4ict [aw=wt1], a(yr) cl(fgroup);

restore;

***********************************************************;
*** Save Files                                          ***;
***********************************************************;

*** Save working files;
sort firm yr; save env-firm1-anal5.dta, replace;
keep firm; duplicates drop;
sort firm; save env-firm1-anal5-firm.dta, replace;

***********************************************************;
*** Moments                                             ***;
***********************************************************;

*** Base moments;
use env-firm1-anal5, clear;
keep if envsamp==1; *keep if yr>=1992;
egen yrmin=min(yr); egen yrmax=max(yr);
gen CM_tvs=dns; gen CM_dtvs=dns; 
for var CM_dtvs: 
\ replace X=X*(113.6/60.6) if yr==1977
\ replace X=X*(113.6/96.5) if yr==1982
\ replace X=X*(113.6/140.3) if yr==1992
\ replace X=X*(113.6/160.5) if yr==1997
\ replace X=X*(113.6/179.9) if yr==2002;
gen CM_dtvsdne=CM_dtvs/dne; gen lCM_dtvs=ln(CM_dtvs);
cap n ren yr year;

*** Table values;
table envcoalx year, row col;
for var rdtot dne dns CM_dtvs pat4 ab:
table envcoalx year, c(sum X) row col;

*** Aggregate sales growth per worker;
egen temp1a=sum(CM_dtvs*(year==1977)); egen temp1b=sum(dne*(year==1977)); gen temp3_77=temp1a/temp1b; 
egen temp2a=sum(CM_dtvs*(year==2002)); egen temp2b=sum(dne*(year==2002)); gen temp3_02=temp2a/temp2b; 
gen agggr=(temp3_02-temp3_77)/((temp3_02+temp3_77)/2);
gen agggr2=(1+agggr)^(1/(25))-1;
gen agggr3=(1+agggr)^(1/(30))-1;
su agggr*; drop agggr* temp*;

*** Generate entry and exit rates;
sort firm year;
gen ZZentry=(firm!=firm[_n-1]);
gen ZZexit=(firm!=firm[_n+1]);
replace ZZentry=. if year==yrmin;
replace ZZexit=. if year==yrmax;
for var ZZentry ZZexit: gen TX=X/5; 
sum TZZentry TZZexit; drop TZZentry TZZexit;

*** Generate entry and exit employment shares;
gen temp0=dne if ZZentry!=.; egen temp1=sum(temp0), by(ZZentry); egen temp2=sum(temp0); replace temp1=temp1/temp2; tab ZZentry, s(temp1); drop temp*;
gen temp0=dne if ZZexit!=.; egen temp1=sum(temp0), by(ZZexit); egen temp2=sum(temp0); replace temp1=temp1/temp2; tab ZZexit, s(temp1); drop temp*;

*** Generate R&D/Sales ratio;
gen zwt=ln(1+CM_dtvs);
gen ZZrdtvs2=rdtot/dns;
gen temp1=10; gen ZZrdtvs3=ZZrdtvs2; replace ZZrdtvs3=temp1 if ZZrdtvs3>temp1 & ZZrdtvs3!=.; drop temp1;
egen temp1=pctile(ZZrdtvs2), p(99); gen ZZrdtvs4=ZZrdtvs2; replace ZZrdtvs4=temp1 if ZZrdtvs4>temp1 & ZZrdtvs4!=.; drop temp1;
sum ZZrdtvs*, d;
sum ZZrdtvs* [aw=CM_dtvs];
sum ZZrdtvs* [aw=zwt];

*** Generate DHS growth (QQ);
sort firm year; gen QQgremp=2*(dne[_n+1]-dne)/(dne[_n+1]+dne) if firm==firm[_n+1] & year==year[_n+1]-5;
for var QQ*: sort firm year \ assert X==. if firm!=firm[_n+1] \ replace X=-2 if firm!=firm[_n+1] & year<yrmax;

*** Growth over size distribution;
log off; 
for num 1/9: egen tempX0=pctile(dne), p(X0) by(year);
gen size=0; for num 2 4 6 8: replace size=X if dne>=tempX0 & dne!=.; tab size;
replace size=. if dne==.; log on;
table size if QQgremp!=-2, c(mean QQgremp); drop size temp*;

*** Entrant size;
log off;
for var dne CM_dtvs CM_dtvsdne:
\ gen temp1a=X if ZZentry==1 \ egen temp2a1=pctile(temp1a), p(45) \ egen temp2a2=pctile(temp1a), p(55) 
\ gen temp2a3=temp1a if temp1a>=temp2a1 & temp1a<=temp2a2 \ egen temp3a=mean(temp2a3)
\ gen temp1b=X if ZZentry==0 \ egen temp2b1=pctile(temp1b), p(45) \ egen temp2b2=pctile(temp1b), p(55)
\ gen temp2b3=temp1b if temp1b>=temp2b1 & temp1b<=temp2b2 \ egen temp3b=mean(temp2b3)
\ gen RX=temp3a/temp3b \ drop temp*; log on;
sum R*;

*** Tabulate growth ranges with unconditional declines;
log off;
for var QQgremp:
\ gen ESh99=X==-2
\ gen ESh75=X<(0.25-1)*2/(1+0.25)
\ gen ESh50=X<(0.50-1)*2/(1+0.50)
\ gen ESh25=X<(0.75-1)*2/(1+0.75)
\ gen EGr25=X>(1.25-1)*2/(1+1.25)
\ gen EGr50=X>(1.50-1)*2/(1+1.50)
\ gen EGr75=X>(1.75-1)*2/(1+1.75)
\ gen EGr99=X>(2.00-1)*2/(1+2.00);
for var ESh* EGr*: replace X=. if QQgremp==.; log on;
sum ESh* EGr*;

*** End of Program;
cap n log close;








************************************************************************************************;
************************************************************************************************;
*** env-firm1-anal4-distr1.do                                                                ***;
*** This file generates the SIC3 distribution used for initial values in simulation          ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cd @ROOT/programs/reall/files/energy1;
cap n log close; log using env-firm1-anal4-distr1.log, replace;

* William Kerr, TGG;
* Environment tests using joint AB+Popp designations;
* Panel Extended to Non-Mfg Industries;
* Last Modified: August 2012;

clear all; *set matsize 11000; *set maxvar 15000; set more off;

***********************************************************;
*** Prepare SIC4 Distribution of Firm                   ***;
***********************************************************;

*** Build long LBD file and collapse on SIC3 and SIC4;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd1976a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using env-firm1-anal4-firm, nok; tab _m; keep if _m==3; drop _m; compress; save temp1-lbd-distr, replace;
forvalues YR1=1977(1)2001 {;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd`YR1'a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using env-firm1-anal4-firm, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd-distr.dta; compress; save temp1-lbd-distr, replace;
};
use lbdnum yr emp cfn firmid bestsic if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd2002a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
replace firm=firmid; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; ren best sic; 
sort firm; merge firm using env-firm1-anal4-firm, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd-distr.dta; compress; drop if firm==""; save temp1-lbd-distr, replace;
ren emp LB_emp; tab yr, s(LB_emp);

*** Collapse on SIC4;
gen sic4=substr(sic,1,4); drop sic; destring sic4, replace force;
gen sic3=int(sic4/10); replace sic3=999 if sic3==.; gen LB_est=1;
collapse (sum) LB_est LB_emp, by(firm sic4 yr) fast; gen I_tot=1;
ren yr year; gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2);
collapse (mean) LB_est LB_emp, by(firm sic4 yr) fast; 
for var LB_est LB_emp: replace X=int(X);
sort firm sic4 yr; save env-firm1-anal4-base-sic4yr, replace; erase temp1-lbd-distr.dta;

***********************************************************;
*** Prepare SIC Distribution of Firm                    ***;
***********************************************************;

*** Pull out patent stocks per firm;
use env-firm1-anal4-base-pat, clear; 
collapse (max) STct2002, by(firm) fast;
sort firm; save temp1-firmpat, replace;

*** Pull out clean/dirty per firm;
use env-firm1-anal4, clear; 
tab envcoalx; drop clean; gen clean=(envcoalx==0);
collapse (mean) clean, by(firm) fast; tab clean;
sort firm; save temp1-firmclean, replace;

*** Create SIC distribution of dirty and clean;
use env-firm1-anal4-base-sic4yr, clear; replace sic4=9999 if sic4==. | sic4==0;
collapse (mean) LB_est LB_emp, by(firm sic4) fast;
gen int sic3=int(sic4/10);
collapse (sum) LB_est LB_emp, by(firm sic3) fast;
for var LB_est LB_emp: egen temp1=sum(X), by(firm) \ gen Xp=X/temp1\ drop temp1; sum;
sort firm; merge firm using temp1-firmpat; erase temp1-firmpat.dta; tab _m; keep if _m==3; drop _m;
sort firm; merge firm using temp1-firmclean; erase temp1-firmclean.dta; tab _m; keep if _m==3; drop _m;
for var STct2002:
\ gen STCest=X*LB_estp if clean==1 \ gen STCemp=X*LB_empp if clean==1
\ gen STDest=X*LB_estp if clean==0 \ gen STDemp=X*LB_empp if clean==0;
collapse (sum) ST*, by(sic3) fast;
for var ST*: replace X=0 if X==.; format ST* %8.0f;
gen SH=(STCemp+STDemp)/STct*1000; sum SH, d; format SH %3.0f;
count; gen cln1=(STCest>STDest); gen cln2=(STCemp>STDemp); gen ct=STCemp+STDemp;
gen sic3x=int(sic3/10)*10;
gen z=(sic3==sic3x);
sum cln* ct;
sum cln* ct if z==1;
sum cln* ct if z==0;
sum cln* ct if z==0 & ct>=120;
sum cln* ct if z==0 & ct<120;
sort sic3; list sic3 STCemp STDemp SH  if z==0 & SH>=5 & ct>=120;
sort sic3; list sic3 STCemp STDemp SH  if z==0 & SH<5 & ct>=120;
sort sic3; save env-firm1-anal4-distr1, replace;

keep if z==0 & ct>=120;
for any 16 271 302 482 536 538 558 566: drop if sic3==X;
sort sic3; list sic3 STCest STDest STCemp STDemp , clean noobs;

*** End of Program;
cap n log close;








************************************************************************************************;
************************************************************************************************;
*** env-firm1-anal5-distr1.do                                                                ***;
*** This file generates the SIC4 distribution used for initial values in simulation          ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cd @ROOT/programs/reall/files/energy1;
cap n log close; log using env-firm1-anal5-distr1.log, replace;

* William Kerr, TGG;
* Environment Initial Conditions;
* Last Modified: Oct 2012;

clear all; set more off;

***********************************************************;
*** Prepare SIC4 Distribution of Firm                   ***;
***********************************************************;

*** Build long LBD file and collapse on SIC3 and SIC4;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd1976a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using env-firm1-anal5-firm, nok; tab _m; keep if _m==3; drop _m; compress; save temp1-lbd-distr, replace;
forvalues YR1=1977(1)2001 {;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd`YR1'a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using env-firm1-anal5-firm, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd-distr.dta; compress; save temp1-lbd-distr, replace;
};
use lbdnum yr emp cfn firmid bestsic if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd2002a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
replace firm=firmid; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; ren best sic; 
sort firm; merge firm using env-firm1-anal5-firm, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd-distr.dta; compress; drop if firm==""; save temp1-lbd-distr, replace;
ren emp LB_emp; tab yr, s(LB_emp);

*** Collapse on SIC4;
gen sic4=substr(sic,1,4); drop sic; destring sic4, replace force; gen LB_est=1;
collapse (sum) LB_est LB_emp, by(firm sic4 yr) fast; gen I_tot=1;
ren yr year; gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2);
collapse (mean) LB_est LB_emp, by(firm sic4 yr) fast; 
for var LB_est LB_emp: replace X=int(X);
sort firm sic4 yr; save env-firm1-anal5-base-sic4yr, replace; erase temp1-lbd-distr.dta;

***********************************************************;
*** Prepare Distribution of Initial Conditions          ***;
***********************************************************;

*** Pull out patent stocks per firm;
use env-firm1-anal5-base-pat, clear; 
collapse (max) STct2002 STWct2002, by(firm) fast;
sort firm; save temp1-firmpat, replace;

*** Pull out clean/dirty per firm;
use env-firm1-anal5, clear; 
tab envcoalx; drop clean; gen clean=(envcoalx==0);
collapse (mean) clean, by(firm) fast; tab clean;
sort firm; save temp1-firmclean, replace;

*** Create PAT-SIC distribution of dirty and clean;
use env-firm1-anal5-base-sic4yr, clear; replace sic4=9999 if (sic4==. | sic4==0);
collapse (mean) LB_est LB_emp, by(firm sic4) fast; 
sort sic4; merge sic4 using env-firm1-anal5-base-xsic.dta, keep(bad); tab _m; keep if _m==3; drop _m;
	*** THIS IS DUE TO DISCLOSURE ISSUES AROUND TOO SMALL CELL COUNTS ***
gen int sic3=int(sic4/10); keep if (sic3>=200 & sic3<=399) | (sic3>=491 & sic3<=499);
*Toggle; gen zsic=sic4;
collapse (sum) LB_est LB_emp (max) bad, by(firm zsic) fast;
for var LB_est LB_emp: egen temp1=sum(X), by(firm) \ gen Xp=X/temp1\ drop temp1; sum;
sort firm; merge firm using temp1-firmpat; erase temp1-firmpat.dta; tab _m; keep if _m==3; drop _m;
sort firm; merge firm using temp1-firmclean; erase temp1-firmclean.dta; tab _m; keep if _m==3; drop _m;
for var STct2002:
\ gen STCest=X*LB_estp if clean==1 \ gen STCemp=X*LB_empp if clean==1
\ gen STDest=X*LB_estp if clean==0 \ gen STDemp=X*LB_empp if clean==0;
for var STWct2002:
\ gen STWCest=X*LB_estp if clean==1 \ gen STWCemp=X*LB_empp if clean==1
\ gen STWDest=X*LB_estp if clean==0 \ gen STWDemp=X*LB_empp if clean==0;

*** Collapse on SIC and generate comparisons;
collapse (sum) ST* (max) bad, by(zsic) fast;
for var ST*: replace X=0 if X==.; format ST* %8.0f;
gen SH=(STCemp+STDemp)/STct*1000; sum SH, d; format SH %3.0f;
gen SHW=(STWCemp+STWDemp)/STWct*1000; sum SHW, d; format SHW %3.0f;
gen cln1=(STCest>STDest); gen cln2=(STCemp>STDemp); gen ct=STCemp+STDemp;
gen cln1w=(STWCest>STWDest); gen cln2w=(STWCemp>STWDemp); gen ctw=STWCemp+STWDemp;
gen gap=STDemp-STCemp; gen gapw=STWDemp-STWCemp; format gap* %8.0f; gen rgapw=gapw/ctw;
gen zsicx=int(zsic/10)*10; gen z=(zsic==zsicx);

*** Describe broad properties;
gen szC0=(STWCemp>=0 & STWCemp<=20);
gen szC1=(STWCemp>20 & STWCemp<=100);
gen szC2=(STWCemp>100 & STWCemp<=500);
gen szC3=(STWCemp>500 & STWCemp!=.);
gen szD0=(STWDemp>=0 & STWDemp<=20);
gen szD1=(STWDemp>20 & STWDemp<=100);
gen szD2=(STWDemp>100 & STWDemp<=500);
gen szD3=(STWDemp>500 & STWDemp!=.);
sum;
sum if bad==1;
keep if bad==0; 
*** SEE ABOVE RE "BAD" ***
sum;
sum gapw rgapw if cln2w==1;
sum gapw rgapw if cln2w==0;

*** Save file;
keep zsic ST* gap*; drop STct2002 STWct2002;
sort zsic; list, clean noobs; list zsic gapw;
sort zsic; save env-firm1-anal5-distr1-output, replace;

*** List extreme cases;
gsort -gapw; list zsic gap gapw in 1/20; 
gsort gapw; list zsic gap gapw in 1/20; 

*** Provide a graph;
replace gapw=gapw/100;
tw kdensity gapw, ytitle("Kernel Density of Product Lines") subtitle("Figure 1: Density Function of Initial Conditions") xtitle("Initial Dirty-Energy Patent Stock Advantage");
graph save density, replace;

*** End of Program;
cap n log close;








************************************************************************************************;
************************************************************************************************;
*** env-firm1-anal4.do                                                                       ***;
*** This file generates some early data elements used in above programs                      ***;
*** Analytical elements are now superceded and summarized in the anal5-distr1 codes          ***;
************************************************************************************************;
************************************************************************************************;

#delimit;
cd @ROOT/programs/reall/files/energy1;
cap n log close; log using env-firm1-anal4.log, replace;

* William Kerr, TGG;
* Environment tests using joint AB+Popp designations;
* Panel Extended to Non-Mfg Industries;
* Last Modified: August 2012;

clear all; *set matsize 11000; *set maxvar 15000; set more off;

***********************************************************;
*** Prepare group of patent tech for env firm sample    ***;
***********************************************************;

*** Prepare new patent cut - flow and stock;
use if (assignee!=. & ctryn=="US") using @ROOT/trips/data-patents/p-work-pat_2008_CES, clear;
gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (ayear>=X-2 & ayear<=X+2); 
sort patent; merge patent using env-data-prep1c2; tab _m; drop if _m==2; gen ab=(_m==3); drop _m;

*** Merge in assignee-firm;
sort assignee; merge assignee using @ROOT/trips/data-patents/BRDC01-MI-comb-best, nok; 
keep if _m==3; gen ct=1; 

*** FIRM SPECIFIC PIECES EXCLUDED ***

*** Generate stock values 1992;
for any 1992:
\ egen STctX=sum(ayear>=X-15 & ayear<=X+2), by(assignee)
\ egen STabX=sum(ayear>=X-15 & ayear<=X+2 & ab==1), by(assignee)
\ egen STclX=sum(ayear>=X-15 & ayear<=X+2 & clean==1), by(assignee)
\ egen STdiX=sum(ayear>=X-15 & ayear<=X+2 & dirty==1), by(assignee);

*** Generate stock values 2002;
for any 2002:
\ egen STctX=sum(ayear>=1975 & ayear<=X+2), by(assignee)
\ egen STabX=sum(ayear>=1975 & ayear<=X+2 & ab==1), by(assignee)
\ egen STclX=sum(ayear>=1975 & ayear<=X+2 & clean==1), by(assignee)
\ egen STdiX=sum(ayear>=1975 & ayear<=X+2 & dirty==1), by(assignee);
sum ST*; tab ayear if yr==.; drop if yr==.;
sum ST* if firm=="0951643307";

*** Collapse data;
collapse (sum) ct ab alt clean dirty POPP (mean) ST*, by(assignee firm yr) fast;
collapse (sum) ct ab alt clean dirty POPP ST*, by(firm yr) fast;

*** FIRM SPECIFIC PIECES EXCLUDED ***

for var ct ab alt clean dirty POPP ST*: replace X=0 if X==.; ren ct pat4;
sort firm yr; compress; save env-firm1-anal4-base-pat, replace; sum;
collapse (sum) ab, by(firm); keep if ab>2;
keep firm; duplicates drop; sort firm; save temp1-firm-pat, replace;

***********************************************************;
*** Build New Operating and R&D Dataset                 ***;
***********************************************************;

*** Build long R&D file;
use @ROOT/programs/reall/working/rad-pat/unbpanel, clear;
keep id year dns dne set rdtot rdfed; egen temp1=rsum(dns dne set rdtot); keep if temp1>0 & temp1!=.; drop temp1;
for num 2/5:
\ append using @ROOT/ghosh/data/rd3/rd200X
\ egen temp1=rsum(dns dne set rdtot)
\ keep if temp1>0 & temp1!=.
\ keep id year dns dne set rdtot rdfed;
for var dns dne rdtot rdfed set: replace X=. if X==0 \ tab year, s(X); replace dns=. if dns==99999999;
ren id firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; duplicates drop; drop if firm==""; 
collapse (sum) dns dne rdtot rdfed set, by(firm year) fast; 
gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2); 
tab year if yr==.; drop if yr==.;
collapse (mean) dns dne rdtot rdfed set, by(firm yr) fast; 
egen temp1=sum(rdtot), by(firm); drop if (temp1==0 | temp1==.); drop temp1;
compress; sort firm yr; save temp1-firm, replace; 
append using temp1-firm-pat; erase temp1-firm-pat.dta; 
keep firm; duplicates drop; sort firm; save temp1-firm-list, replace;

*** Build long LBD file;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd1976a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m; compress; save temp1-lbd, replace;
forvalues YR1=1977(1)2001 {;
use lbdnum yr emp cfn sic* if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd`YR1'a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd.dta; compress; save temp1-lbd, replace;
};
use lbdnum yr emp cfn firmid bestsic if emp>0 & emp!=. using @ROOT/data/lbd11/data/lbd2002a, clear;
ren cfn firm; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; 
replace firm=firmid; replace firm=substr(firm,1,6)+"0000" if substr(firm,1,1)!="0"; ren best sic; 
sort firm; merge firm using temp1-firm-list, nok; tab _m; keep if _m==3; drop _m;
append using temp1-lbd.dta; compress; drop if firm==""; save temp1-lbd, replace;
ren emp LB_emp; tab yr, s(LB_emp);

*** Collapse and merge R&D;
gen sic4=substr(sic,1,4); drop sic; destring sic4, replace force;
gen sic3=int(sic4/10); replace sic3=999 if sic3==.; *keep if (sic3>=200 & sic3<400); gen LB_est=1;
collapse (sum) LB_emp LB_est, by(firm sic4 yr) fast; 
gen sic3=int(sic4/10); replace sic3=999 if sic3==.; gen I_tot4=1;
collapse (sum) LB_emp LB_est I_tot4, by(firm sic3 yr) fast; gen I_tot=1;
collapse (sum) LB_emp LB_est I_tot4 I_tot, by(firm yr) fast; 
ren yr year; gen yr=.; for any 1977 1982 1987 1992 1997 2002: replace yr=X if (year>=X-2 & year<=X+2);
collapse (mean) LB_emp LB_est I_tot I_tot4, by(firm yr) fast; 
sort firm yr; merge firm yr using temp1-firm; tab _m; drop _m;
for var rdtot set: replace X=0 if X==.;
pwcorr LB_emp dne; sum LB_emp dne;

*** Save base file and clean-up;
compress; sort firm yr; save env-firm1-anal4-base-ces, replace; sum; sum, d;
for any temp1-firm temp1-firm-list temp1-lbd: erase X.dta;

***********************************************************;
*** Prepare group of patent tech for env firm sample    ***;
***********************************************************;

*** Open CES base and define year range;
use env-firm1-anal4-base-ces, clear; 
for var LB_emp LB_est I_tot I_tot4 rdtot rdfed dns dne: replace X=int(X);
format LB_emp LB_est I_tot I_tot4 rdtot rdfed dns dne %8.0f; sum;
*keep if (yr>=1987 & yr<=1997); 

*** Baseline work on variables;
for var dne LB_emp: replace X=1.5e6 if X>1.5e6 & X!=.;
ren dne zdne; egen dne=rmean(zdne LB_emp); 
replace dne=zdne if LB_emp==. & dne==.;
replace dne=LB_emp if zdne==. & dne==.;
drop if dne==.; pwcorr dne zdne LB_emp;
replace I_tot=1 if (I_tot==. | I_tot==0) & dne>0 & dne!=.; 
replace I_tot4=1 if (I_tot4==. | I_tot4==0) & dne>0 & dne!=.; sum; 

*** Merge in patent data to core;
sort firm yr; merge firm yr using env-firm1-anal4-base-pat; tab _m; drop if _m==2; drop _m;
for var pat4 alt clean dirty POPP: replace X=0 if X==.; 

*** Impute small R&D values where set and pat are known;
sum pat4 set rdtot if rdtot>0 & rdtot<1000;
sum pat4 set rdtot if rdtot==0; 
gen rdproxy=264.1/4.9*set if rdtot==0 & set>0; 
replace rdproxy=264.1/1.5*pat4 if rdproxy==. & rdtot==0 & pat4>0;
replace rdproxy=1000 if rdproxy>1000 & rdproxy!=.;
replace rdproxy=. if rdtot!=0;
replace rdtot=rdproxy if rdtot==0 & rdproxy!=.;

*** Apply new conditional rule of entry with R&D;
egen yrmin=min(yr), by(firm); 
gen flag=(yr==yrmin & rdtot==0); tab flag; tab flag if yr==yrmin;
egen zdrop=max(flag), by(firm); count; drop if zdrop==1; count;
drop yrmin zdrop flag;

***********************************************************;
*** Prepare Variables for Estimations                   ***;
***********************************************************;

*** Define environmental sample;
egen temp1=sum(ab), by(firm); egen temp2=sum(pat4), by(firm); gen absh=temp1/temp2; 
gen temp3=(temp1>2 & absh>0.10 & absh!=.); egen envsamp=max(temp3), by(firm); 
gen temp4=(temp1>15 & absh>0.25 & absh!=.); egen envsamp2=max(temp4), by(firm); drop temp*;

*** FIRM SPECIFIC PIECES EXCLUDED ***

*** Define group types;
gen temp1=(dirty/(dirty+clean)>0.75); egen envcoal=max(temp1), by(firm yr);
egen temp2=mean(temp1), by(firm); gen envcoalx=(temp2>0.5); drop temp*;
gen temp1=(alt/(dirty+clean)>0.15); egen envalt=max(temp1), by(firm yr);
egen temp2=mean(temp1), by(firm); gen envaltx=(temp2>0.5); drop temp*;
for var envcoal* envalt*: replace X=. if envsamp==0; 
for var rdtot rdfed: table envcoalx yr if rdproxy==., c(sum X) row col f(%8.0f);
gen rdfedsh=rdfed/rdtot; for num 0 1: sum rdfed rdfedsh if envcoalx==X, d;
tab firm if rdfed>500000 & rdfed!=. & envcoalx==0;
tab firm if rdfed>500000 & rdfed!=. & envcoalx==1;

*** Variable prep;
egen patmin=min(pat4), by(firm); sum patmin, d;
gen rddns=rdtot/dns; gen rdemp=rdtot/dne; gen rdict=rdtot/I_tot; gen rd4ict=rdtot/I_tot4; gen rdest=rdtot/LB_est;
gen setdns=set/dns; gen setemp=set/dne; gen setict=set/I_tot; gen set4ict=set/I_tot4; gen setest=set/LB_est;
gen patdns=pat4/dns; gen patemp=pat4/dne; gen patict=pat4/I_tot; gen pat4ict=pat4/I_tot4; gen patest=pat4/LB_est;
for var set* dns dne LB_emp LB_est I_tot I_tot4 pat* rd*:
\ sort firm yr \ gen lX=ln(X) \ gen DlX=lX-lX[_n-1] if firm==firm[_n-1] & yr==yr[_n-1]+5;
for var rdtot rddns rdemp rdict setemp: sort firm yr \ gen LlX=ln(X[_n-1]) if firm==firm[_n-1];
for var lrdtot lLB_emp LB_emp: egen mn_X=mean(X), by(firm); 

*** Winsorize;
log off;
*for var Dl*:
\ egen temp1=pctile(X), p(02) \ replace X=temp1 if X<temp1 & X!=. \ drop temp1
\ egen temp1=pctile(X), p(98) \ replace X=temp1 if X>temp1 & X!=. \ drop temp1;
*for var lrdtot lrdemp Dlrdemp lrdict Dlrdict: gen XO=X*(1-envcoal) \ gen XA=X*(envalt);
log on;

*** Final details;
egen cyr=group(envcoal yr); egen ayr=group(envalt yr);
egen fgroup=group(firm); gen wt1=1; gen wt2=mn_lLB_emp; set seed 1;

***********************************************************;
*** Show Distributions and Eta Elasticity Calculations  ***;
***********************************************************;

*** PIECES EXCLUDED, SUPERCEDED BY ANAL5-REGR1 ABOVE ***

***********************************************************;
*** Save Files and Extra Checks                         ***;
***********************************************************;

*** Save working files;
sort firm yr; save env-firm1-anal4.dta, replace;
keep firm; duplicates drop;
sort firm; save env-firm1-anal4-firm.dta, replace;

*** FIRM SPECIFIC PIECES AND OBSOLETE CODE EXCLUDED ***

*** End of Program;
cap n log close;