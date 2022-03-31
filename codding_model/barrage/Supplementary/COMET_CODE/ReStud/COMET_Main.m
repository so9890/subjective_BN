%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%COMET Model M-File%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Lint Barrage, Brown University
%%Version: August 2018

clear
%cd C:\Users\lintb\Desktop\COMET\ReStud
cd '/home/sonja/Documents/projects/Overconsumption/codding_model/barrage/Supplementary/COMET_CODE/ReStud'
%M-File Outline%
%%%%%%%%%%%%%%%%
%General Notes
%Section 1: Select Fiscal Scenario
%Section 2: Set Parameters
%Section 3: Solve for Optimal Allocation
%Section 4: Compute Implementing Tax Rates and Outcomes
            %Generates main results for Tables 4, 5, 7, A5, A6
%Section 5: Welfare Calculations
%Section 6: Saved Results and Figures
            %Generates paper figures from saved results
            
%General Notes:
%%%%%%%%%%%%%%%
%For a complete description of parameter sources, please see Online Appendix.
%This file needs to be run SEPARATELY for each fiscal or parameter scenario.
%That is, it requires manually changing the scenario and re-running to produce 
%and save ONE set of model results at a time. Saved results from all model
%runs referenced in the paper and Online Appendix are included and used
%below to re-generate graphs and welfare calculations, as initial guesses, 
%and as calibration inputs (e.g., carbon tax constraints to first-best level).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 1: Select Fiscal Scenario        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 25;  %Direct optimization period time horizon, 1 period = 10 years

%Note: A value MUST be selected for each of the following:

%1.1) Distortionary taxes or lump-sum taxes?
distortionary = 1;              %Revenue must be raised through distortionary taxes (Table 4 Scenarios 1a,1b,2,3a,4a,3b,4b,5; Table 5 Scenarios 4b',4a'; Table 6 Scenarios 4a,4b,5; Table 7 All Scenarios; Online Appendix Table A5 All Scenarios, Online Appendix Table A6 Scenarios 5,4a,4b)
%distortionary = 0;             %Lump-sum taxation permitted (Table 4 Scenario 6; Table 6 Scenario 6; Online Appendix Table A6 Scenario 6)

%1.2) Restrict labor or capital tax levels (wedges):
if distortionary==1
   %%Labor income tax (labor-consumption wedge) level constrained?
         tao_l_fix = [];        %No constraint (Table 4 Scenarios 1a,2,3b,4b,5; Table 5 Scenario 4b'; Table 6 Scenarios 5,4b; Table 7 Scenarios 1a,3b,4b; Online Appendix Table A5 4b,4b',4b''; Table A6 Scenarios 5,4b)
%       tao_l_fix = .3825;     %Constrained to 38.25% (Table 4 Scenarios 1b,3a,4a; Table 5 Scenario 4a'; Online Appendix Table A5 Scenarios 4a,4a',4a'')
%       tao_l_fix = .39;       %Constrained to 39.00% (Table 6, Frisch=2, Scenario 4a; Table 7, Frisch=2,  Scenarios 3a,4a)
%       tao_l_fix = .4075;     %Constrained to 40.75% (Table 6, Sigma=2,  Scenario 4a; Table 7, Sigma=2,   Scenarios 3a,4a)
%       tao_l_fix = .35;       %Constrained to 35.00% (Table 6, Sigma=1.1,Scenario 4a; Table 7, Sigma=1.1, Scenarios 3a,4a)
%       tao_l_fix = .32;       %Constrained to 32.00% (Table 6, G-10%,    Scenario 4a; Table 7, G-10%,     Scenarios 3a,4a)
%       tao_l_fix = .45;       %Constrained to 45.00% (Table 6, G+10%,    Scenario 4a; Table 7, G+10%,     Scenarios 3a,4a)
%       tao_l_fix = .3775      %Constrained to 37.75% (Online Appendix Table A6 Scenario 4a)
    
    %%%Capital income tax (intertemporal wedge) level constrained?
        tao_k_fix = [];         %No constraint (Table 4 Scenarios 1b,2,3a,4a,5; Table 5 Scenario 4a'; Table 6 Scenarios 5,4a; Table 7 Scenarios 1a,3a,4a; Online Appendix Table A5 4a,4a',4a''; Table A6 Scenarios 5,4a)
%      tao_k_fix = .3457;      %Constrained to 34.57% (Table 4 Scenarios 1a,3b,4b; Table 5 Scenario 4b'; Online Appendix Table A5 Scenarios 4b,4b',4b'')
%      tao_k_fix = .409;       %Constrained to 40.09% (Table 6, Frisch=2, Scenario 4b; Table 7, Frisch=2,  Scenarios 3b,4b)
%      tao_k_fix = .394;       %Constrained to 39.40% (Table 6, Sigma=2,  Scenario 4b; Table 7, Sigma=2,   Scenarios 3b,4b)
%      tao_k_fix = .332;       %Constrained to 33.20% (Table 6, Sigma=1.1,Scenario 4b; Table 7, Sigma=1.1, Scenarios 3b,4b)
%      tao_k_fix = .296;       %Constrained to 29.60% (Table 6, G-10%,    Scenario 4b; Table 7, G-10%,     Scenarios 3b,4b)
%      tao_k_fix = .38;        %Constrained to 38.00% (Table 6, G+10%,    Scenario 4b; Table 7, G+10%,     Scenarios 3b,4b)
%      tao_l_fix = .3457       %Constrained to 34.57% (Online Appendix Table A6 Scenario 4b)
else
    tao_l_fix = [];            %Table 4 Scenario 6; Table 6 Scenario 6; Online Appendix Table A6 Scenario 6
    tao_k_fix = [];            %Table 4 Scenario 6; Table 6 Scenario 6; Online Appendix Table A6 Scenario 6
end
    
%1.3) Restrict labor or consumption tax flexibility over time?
if distortionary==1
    %%%Labor income tax constant over time?
      tao_l_const = 0;          %Flexible labor wedge over time (Table 4 Scenarios 2,5, also select for 1b,3a,4a to avoid unnecessary constraints; Table 5 Scenario 4b'; Table 6 Scenario 5 and also select for 4a; Table 7 Scenario 5; Online Appendix Table A5 Scenario 4b'; Table A6 Scenario 5)
      %tao_l_const = 1;         %Constrained to be constant over time (Table 4 Scenarios 1a,3b,4b; Table 6 Scenario 4b; Table 7 Scenarios 1a,3b,4b; Online Appendix Table A5 4b,4b''; A6 Scenarios 4b)

   %%%Capital income tax constant over time?
      tao_k_const = 0;           %Flexible capital wedge over time (Table 4 Scenarios 2,5, also select 1a,3b,4b to avoid unnecessary constraints; Table 5 Scenario 4a'; Table 6 Scenario 5 and also select for 4b; Table 7 Scenario 5; Online Appendix Table A5 Scenario 4a'; Table A6 Scenario 5)
      %tao_k_const = 1;          %Constrained to be constant over time  (Table 4 Scenarios 1b,3a,4a; Table 6 Scenario 4a; Table 7 Scenarios 3a,4a; Online Appendix Table A5 4a,4a''; Table A6 Scenario 4a)
else
    tao_l_const = [];          %Table 4 Scenario 6; Table 6 Scenario 6; Online Appendix Table A6 Scenario 6
    tao_k_const = [];          %Table 4 Scenario 6; Table 6 Scenario 6; Online Appendix Table A6 Scenario 6
end    
   
%1.4) Restrict carbon prices (direct and via intermediate energy good tax levels)? 
%Note: Both the direct emissions tax taoE and intermediate energy goods tax (via energy_wedge = taoE+taoI) 
%need to be selected and adjusted in tandem to ensure overall wedges are at desired constraint levels (e.g., no carbon price).
  tao_E_fix = [];                         % No constraint (optimized) (Table 4 Scenarios 4a,4b,5,6; Table 5 Scenarios 4a',4b'; Table 6 All Scenarios; Table 7 Scenarios 4b,4a,5; Online Appendix Table A5 All Scenarios; Table A6 All Scenarios)
  energy_wedge_fix = [];                  % No constraint (optimized) (Table 4 Scenarios 4a,4b,5,6; Table 5 Scenarios 4a',4b'; Table 6 All Scenarios; Table 7 Scenarios 4b,4a,5; Online Appendix Table A5 All Scenarios; Table A6 All Scenarios)
% tao_E_fix = 0.01*ones(T,1);             % BAU carbon price (Table 4 Scenarios 1a,1b,2)
% energy_wedge_fix = 0.01*ones(T,1);      % BAU carbon price (Table 4 Scenarios 1a,1b,2)
% load('MAC_LS','MAC_LS')                 %`Wrong' first-best carbon price (Table 4 Scenarios 3a,3b)
% tao_E_fix = MAC_LS;                     %`Wrong' first-best carbon price (Table 4 Scenarios 3a,3b)
% load('EWedge_LS','EWedge_LS')           %`Wrong' first-best carbon price (Table 4 Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS;           %`Wrong' first-best carbon price (Table 4 Scenarios 3a,3b)
% load('MAC_LS_Fr2','MAC_LS_Fr2')         %`Wrong' first-best carbon price (Table 7, Frisch=2, Scenarios 3a,3b)
% MAC_fix = MAC_LS_Fr2;                   %`Wrong' first-best carbon price (Table 7, Frisch=2, Scenarios 3a,3b)
% load('EWedge_LS_Fr2','EWedge_LS_Fr2')   %`Wrong' first-best carbon price (Table 7, Frisch=2, Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS_Fr2;       %`Wrong' first-best carbon price (Table 7, Frisch=2, Scenarios 3a,3b)
% load('MAC_LS_sig2','MAC_LS_sig2')       %`Wrong' first-best carbon price (Table 7, Sigma=2, Scenarios 3a,3b)
% MAC_fix = MAC_LS_sig2;                  %`Wrong' first-best carbon price (Table 7, Sigma=2, Scenarios 3a,3b)
% load('EWedge_LS_sig2','EWedge_LS_sig2') %`Wrong' first-best carbon price (Table 7, Sigma=2, Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS_sig2;      %`Wrong' first-best carbon price (Table 7, Sigma=2, Scenarios 3a,3b)
% load('MAC_LS_sig11','MAC_LS_sig11')     %`Wrong' first-best carbon price (Table 7, Sigma=1.1, Scenarios 3a,3b)
% MAC_wedge_fix = MAC_LS_sig11;           %`Wrong' first-best carbon price (Table 7, Sigma=1.1, Scenarios 3a,3b)
% load('EWedge_LS_sig11','EWedge_LS_sig11')%`Wrong' first-best carbon price (Table 7, Sigma=1.1, Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS_sig11;      %`Wrong' first-best carbon price (Table 7, Sigma=1.1, Scenarios 3a,3b)
% load('MAC_LS_g09','MAC_LS_g09')          %`Wrong' first-best carbon price (Table 7, G-10%, Scenarios 3a,3b)
% MAC_wedge_fix = MAC_LS_g09;              %`Wrong' first-best carbon price (Table 7, G-10%, Scenarios 3a,3b)
% load('EWedge_LS_g09','EWedge_LS_g09')    %`Wrong' first-best carbon price (Table 7, G-10%, Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS_g09;        %`Wrong' first-best carbon price (Table 7, G-10%, Scenarios 3a,3b)
% load('MAC_LS_g11','MAC_LS_g11')          %`Wrong' first-best carbon price (Table 7, G+10%, Scenarios 3a,3b)
% MAC_wedge_fix = MAC_LS_g11;              %`Wrong' first-best carbon price (Table 7, G+10%, Scenarios 3a,3b)
% load('EWedge_LS_g11','EWedge_LS_g11')    %`Wrong' first-best carbon price (Table 7, G+10%, Scenarios 3a,3b)
% energy_wedge_fix = EWedge_LS_g11;        %`Wrong' first-best carbon price (Table 7, G+10%, Scenarios 3a,3b)
if isempty(tao_E_fix)==0
     T_tao_E_fix = 10;                    %Constrain carbon price, energy wedge for 10 periods (100 years) (Table 4 Scenarios 1a,1b,2; Table 7 Scenario 1a)
%       T_tao_E_fix = T;                  %Constrain carbon price, energy wedge for T periods (250 years) (Table 4 Scenarios 3a,3b; Table 7 Scenario 3a,3b)
    else 
    T_tao_E_fix = [];
end
    
%1.5) Restrict intermediate energy input taxes/subsidy (marginal energy product-price wedge above and beyond carbon tax)?
    no_interm_tax = 0;                   %Allow intermediate energy input tax/subsidy (ALL Scenarios except Online Appendix Table A5 Scenarios 4a'', 4b''))
   %no_interm_tax = 1;                   %Prohibit intermediate energy input tax/subsidy (except through carbon emissions price) (Online Appendix Table A5 Scenarios 4a'', 4b'')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 2: Parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calibration Base Year Tax Rates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 tao_k_0 = .3340;           %Base year capital income tax rate
 tao_l_0 = .3609;           %Base year labor income tax

%%% Time %%%
%%%%%%%%%%%%
periods = 10;               %Simulation periods after time T
y = (1:1:T);                %Time in calendar years
for i = 2:1:T+periods;
     y(i) = 2015+((i-1)*10);
end
y(1) = 2015;

%%% Households %%%
%%%%%%%%%%%%%%%%%%
beta = (.985)^10;            %Decadal pure rate of social time preference
 sigma = 1.5;                %Inverse of intertemporal elasticity of substitution
% sigma = 2;                 %Sigma=2   - Tables 6, 7
% sigma = 1.1;               %Sigma=1.1 - Tables 6, 7

%%% Population %%%
%%%%%%%%%%%%%%%%%%
N0 = 6.411;                 %Base-year population, billions of people
NMax = 8.7;                 %Asymptotic population, billions of people
N = ones((T+periods),1);    %Population vector
arPop = 0.485;              %Population growth rate parameter taken from 2010-DICE
N(1) = N0*((NMax/N0)^(arPop));
    for i = 2:1:T;
        N(i)= N(i-1)*((NMax/N(i-1))^(arPop));
    end
   for i = 1:1:periods;     
        N(T+i) = N(T);
   end
gPop = ones(T+periods,1);   %Annual population growth rate
gPop0 = (N(1)/N0)-1;
    for i = 2:1:T;
        gPop(i-1) = (N(i)/N(i-1))-1;
    end
    gPop(T) = 0;
    for j = 1:1:periods;
        gPop(T+j) = gPop(T);
    end
 
%%% Carbon Cycle and Climate System %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [From 2010-DICE Model]
lambda = 3.2;               %Climate sensitivity parameter
Sbar = 596.4;               %Year 1750 atmospheric carbon concentration in GtC
S_t0 = 787;                 %Base year (2005) atmospheric concentrations in GtC
ELand0 = 1.6;               %Exogenous land-based emissions in GtC per year
ELand0 = ELand0*10;         %Exogenous land-based emissions in GtC per decade
gELand = 0.8;               %Land emissions growth parameter
    ELand = zeros(T+periods,1);
    for i = 1:1:(T+periods);
        ELand(i) = ELand0*(gELand)^(i-1);
        
    end
Fx2000 = 0.83;              %Year 2000 exogenous forcings, watts/m^2
Fx2100 = .3;                %Year 2100 exogenous forcings, watts/m^2
Fx = (ones((T+periods),1))*Fx2100; 
for i = 1:1:10;
    Fx(i) = Fx2000+0.1*(Fx2100-Fx2000)*(i);
end
Qt0 = 0.00680;              %Lower ocean temperature change between 1900 and base year 2005, degrees C
Zt0 = 10010.493;            %Base year 2005 lower ocean carbon concentrations, GtC
X0 = 1600;                  %Base year 2005 upper ocean/biosphere carbon concentrations, GtC
FCO22x = 3.8;               %Temperature-forcing parameter in watts/meter^2
ksi1 = .208;                %Speed of adjustment for atmospheric temperature
ksi2 = (FCO22x/lambda);     
ksi3 = .31;                 %Coefficient of heat loss from atmosphere to oceans
ksi4 = 0.05;                %Coefficient of heat gain by deep oceans
phi11 = 88/100;             %Atmosphere to atmosphere carbon cycle transition coefficient (% per decade)
phi21 = 4.704/100;          %Upper ocean/biosphere to atmosphere carbon cycle transition coefficient (% per decade)
phi12 = 12/100;             %Atmosphere to upper ocean/biosphere carbon cycle transition coefficient (% per decade)
phi22 = 94.796/100;         %Upper ocean/biosphere to upper ocean/biosphere carbon cycle transition coefficient(% per decade)
phi32 = 0.075/100;          %Lower ocean to upper ocean/biosphere carbon cycle transition coefficient(% per decade)
phi23 = 0.5/100;            %Upper ocean/biosphere to lower ocean carbon cycle transition coefficient(% per decade)
phi33 = 99.925/100;         %Lower ocean to lower ocean carbon cycle transition coefficient(% per decade)
eta = FCO22x;
TC0 = .83;                  %Mean atmospheric surface temperature change in 2005 in deg. Celsius over 1900 levels
    

%%% Climate Change Production Damages %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tot_Damages_2p5 = 1.744043233;              %Total output-equivalent losses due to 2.5C warming (%GDP)
Y_Damages_Share = 0.73848211;               %Share of climate impacts at 2.5C warming in the BAU scenario affecting production possibilities
%Y_Damages_Share = 1;                       %100% Ouptut Damages  - Figure 4
%Y_Damages_Share = 0.5;                     %50%  Ouptut Damages  - Figure 4
%Y_Damages_Share = 0.25;                    %25%  Ouptut Damages  - Figure 4
%Y_Damages_Share = 0;                       %0%   Ouptut Damages  - Figure 4
ydamage_2p5C = Y_Damages_Share*Tot_Damages_2p5; %Output damages at 2.5C warming in BAU scenario (%GDP)
theta2 = 2;                                 %Damage function curvature
%Solve for output damage parameter theta1 via: (1/(1+theta1*(2.5C)^theta2)=(1-ydamage(2.5C))
theta1 = (((1/(1-(ydamage_2p5C/100)))-1)/(2.5^theta2));


%%% Production Functions & Initial Capital, Labor Allocations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = .3;                                 %Final goods sector capital expenditure share
v = 0.03;                                   %Final goods sector energy expenditure share
alphaE = 0.403;                             %Energy sector labor expenditure share
E0 = 7.990;                                 %Base year carbon-energy inputs, bmtC/year [year 2005 initial value from DICE 2010]
delta = .1;                                 %Annual capital depreciation rate
Delta = 1-((1-delta)^10);                   %Decadal capital depreciation rate
Ysm1 = 55340;                               %Base year 2005 annual GDP, bil. of int. 2005 PPP dollars
rsm1 = 0.05;                                %Base year annual net rate of return on capital
K0 = (alpha*Ysm1)/(rsm1+delta);             %Total base year capital stock, bil. of int. 2005 PPP dollars
L0 = .2272;                                 %Base year 2005 time share spent on work 
pi_l_0 = (1-alpha-v)/(alphaE*v+1-alpha-v);  %Base year labor share in final goods production consistent with profit maximization and energy production E0
pi_k_0 = alpha/(((1-alphaE)*v)+alpha);      %Base year capital share in final goods production consistent with profit maximization and energy production E0
K0_FG = K0*(pi_k_0);                        %Base year final goods production capital stock, bil. of int. 2005 PPP dollars
KE0 = K0*(1-pi_k_0);                        %Base year energy sector capital stock, bil. of int. 2005 PPP dollars
LE0 = L0*N0*(1-pi_l_0);                     %Base year energy sector labor input

%%% Productivity Growth Rates %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gZ0 = 0.160023196685654;                    %Base year TFP growth rate (% per decade)
dgz0 = 0.00942588385340332;                 %Rate of decline in TFP growth rate (% per year)
ddgz0 = 0.00192375245926376;                %Rate of decline in rate of decline in TFP growth rate (% per year)
    gZt = (ones(T+periods,1))*gZ0;          %Vector of decadal TFP growth rates
    gXt = zeros(T+periods,1);               %Vector of decadal labor productivitiy growth rates
    gX0 = ((1+gZ0))^(1/(1-alpha-v))-1;
    for i = 1:1:T;
        gZt(i) = gZ0*exp(((-dgz0)*10*(i))*exp((-ddgz0)*10*(i)));
        gXt(i) = ((1+gZt(i))^(1/(1-alpha-v)))-1;
    end
for j = 1:1:periods;
    gXt(T+i) = gXt(T);
    gZt(i) = gZt(T);
end

%%%Sectoral Productivity Levels%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Final Goods Sector:
Z0 = (Ysm1*10)/((((1+theta1*(TC0)^2))^(-1))*(((K0_FG)^(alpha))*((L0*pi_l_0*N0)^(1-alpha-v))*((E0*10)^(v))));  %Base year final goods sector initial TFP (decadal)
Z = (ones(T,1));    
    Z(1) = Z0*(1/(1-gZ0));
    for i = 2:1:T+periods;
        Z(i) = Z(i-1)*(1/(1-gZt(i-1)));
    end
%Energy Sector:
A_E0 = (E0*10)/((KE0^(1-alphaE))*(LE0^alphaE)); %Base year energy sector TFP (decadal)
A_E = zeros(T+periods,1);
A_E(1) = A_E0*((1+gX0)^alphaE);
for i = 2:1:T+periods;
  A_E(i) = A_E(i-1)*((1+gXt(i-1))^alphaE);
end

%%% Model Time Zero Capital %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s0 = 24.554/100;                        %Base year 2005 savings rate (DICE)
K1 = K0*(1-Delta)+s0*(Ysm1*10);         %Model time zero (2015) total capital stock, bil. of int. 2005 PPP dollars
    
%%% Government %%%
%%%%%%%%%%%%%%%%%%
 GovSpenfr = .311;                       %Base year 2005 government expenditure/GDP ratio
% GovSpenfr = .311*(0.9);                %-10% Government Spending - Tables 6, 7
% GovSpenfr = .311*(1.1);                %+10% Government Spending - Tables 6, 7
trns_gdp = .1332;                        %Base year 2005 government transfers/GDP ratio
trans_share = trns_gdp/GovSpenfr;        %Transfer share of government expenditures
gov_c_fr = (GovSpenfr-trns_gdp);         %Government consumption share of government expenditures
Gsm1 = (GovSpenfr)*Ysm1;                 %Base year government expenditure per year, bil. of int. 2005 PPP dollars
G = ones(T+periods,1);
Gct = ones(T+periods,1);                 %Government consumption per decade, bil. of int. 2005 PPP dollars
G(1) = Gsm1*exp(gPop0+gX0);              %Government expenditures are assumed to grow at rates of population and labor productivity growth
for i = 2:1:T+periods;
    G(i) = G(i-1)*exp((gPop(i-1)+gXt(i-1)));  
end
for i = 1:1:T+periods;
    Gct(i) = G(i)*10*(1-trans_share);
end
B0 = (.6263)*Ysm1;                       %Base year government debt, bil. of int 2005 PPP dollars 

%%% Abatement / Clean Energy Technology Costs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pc0 = 1.26;                              %Base year price of backstop technology, thousands of dollars per mtC ($2005)
Pc = zeros((periods+T),1);
YearComp = 2250;                         %Year when the backstop technology becomes cost competitive (no additional cost for clean energy)
CostlyYears = YearComp - 2010;
CostlyPeriods = CostlyYears/10;
vau1 = 2;                               
vau2 = 0.05;
for i = 1:1:(CostlyPeriods);
    Pc(i) = Pc0*(vau1-1+exp((-vau2)*(i)))*(1/vau1);
end
 gamma = 0.9245;                         %="a-bar" in paper equation (44)
 a1 =  49.9096;
 a2 =  15.5724;
 a3 =  3.4648;
 a4 =  1.0995;
 b1 = 13.0210;
 b2 = 1.5725;
 b3 = 0.0921;
%Pc = Pc*(0.8);                         %-20% Abatement Cost - Online Appendix Table A6
%Pc = Pc*(1.2);                         %+20% Abatement Cost - Online Appendix Table A6


%%Labor Preferences%%
%%%%%%%%%%%%%%%%%%%%%
Frisch = 0.78;                          %Frisch elasticity of labor supply
%Frisch = 2;                            %Frisch elasticity = 2 - Tables 6, 7           
hh_cons_exp = 1-s0-gov_c_fr;            %Base year fraction of output consumed
y_pc0 = ((Ysm1*10)/N0)/10000;           %Base year implied per capita output per decade, $10,000s int. 2005 PPP dollars
c_05 = hh_cons_exp*y_pc0;               %Base year implied per capita consumption per decade, $10,000s int. 2005 PPP dollars
w_05 = (1-alpha-v)*((Ysm1*10)/(L0*N0))*(1/10000);   %Base year implied wage (marginal product of labor) per decade, $10,000s int. 2005 PPP dollars
temp_term = (1-sigma-(((1-sigma)^2)/(-sigma)));
%Calibrating to the Frisch elasticity and to rationalize base year labor supply yields utility parameters (see paper for details):
if distortionary==1
    phi_labor = (((1-alpha-v)*(1-tao_l_0)*y_pc0*temp_term)+(c_05/Frisch))/(L0*(1-alpha-v)*(1-tao_l_0)*y_pc0*temp_term+c_05*L0+L0*(1/Frisch)*c_05);
else
     phi_labor = (((1-alpha-v)*y_pc0*temp_term)+(c_05/Frisch))/(L0*(1-alpha-v)*y_pc0*temp_term+c_05*L0+L0*(1/Frisch)*c_05);
end
gamma_labor = ((((1-phi_labor*L0)/(phi_labor*L0))*(-1/Frisch))+1)/(1-sigma-(((1-sigma)^2)/(-sigma)));


%%% Climate Change Utility Impacts %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Udam_2p5 = ((1-Y_Damages_Share)*Tot_Damages_2p5)/100;   %Percent of output loss-equivalent of utility damages at 2.5C
alpha1 = 2;                                             %Curvature on temperature change in utility
%To calibrate willingness-to-pay, need consumption at calibration point (2.5C in BAU, between 2045-55):
N_2p5 = 8577049991;                                     %Population at 2.5C in 2010-DICE BAU (weighted avg. between 2045 and 2055)
Y_g_2p5 = 224573314409384;                              %Gross output in int. 2005 PPP dollars per year at 2.5C in 2010-DICE BAU (weighted avg. between 2045 and 2055)
s_2p5 = (22.02190532/100);                              %Gross savings rate at 2.5C in 2010-DICE BAU (weighted avg. between 2045 and 2055)
y_g_2p5 = Y_g_2p5/N_2p5;                                %Gross output per capita in int. 2005 PPP dollars per year at 2.5C in 2010-DICE BAU (weighted avg. between 2045 and 2055)
c_g_2p5 = (1-s_2p5-gov_c_fr)*y_g_2p5*10*(1/10000);      %Per capita consumption (net of savings and government consumption)
n_2p5 = L0;                                             %Baseline labor supply fraction
temp_term2 = (c_g_2p5^(1-sigma))*((1-phi_labor*n_2p5)^(gamma_labor*(1-sigma)));
temp_term3 = (((1-Udam_2p5)^(1-sigma))-1);
alpha0 = (((temp_term2*temp_term3+1)^(-1/(1-sigma)))-1)/(2.5^alpha1);   %Climate change disutility parameter

%%%Miscellaneous%%%
multip = 0;     %Placeholder for multiplier used in welfare calculations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 3: Solve for Optimal Allocation        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preview: Structure of vector of allocations x:
% C = x(1:T);               Per capita consumption ($10,000's)
% L = x(T+1:2*T);           Total labor supply (fraction of time)
% E = x(2*T+1:3*T);         Total energy inputs (GtC-eq)
% pi1_l = x(3*T+1:4*T);     Fraction of labor devoted to final goods sector
% K1t = x(4*T+1:5*T);       Capital in final goods production ($10,000's per capita)
% ECleanPct = x(5*T+1:6*T); Fraction clean energy / abatement 
% K2t = x(6*T+1+1:7*T+1);   Capital in energy production ($10,000's per capita)
% sT = x(6*T+1);            Continuation savings rate (fraction)

%%% Set Bounds %%%
%%%%%%%%%%%%%%%%%%
lb = zeros((7*T)+1,1);
ub = ones((7*T)+1,1);
ub = Inf*ub;
for j = 0:1:T-1;
   ub(T+1+j) = (1/phi_labor)-0.01;
   ub(3*T+1+j) = 0.9999;
   ub(5*T+1+j) = 1;
   ub(6*T+1) = 0.9999;
   lb(1+j) = 0.1;
   lb(T+1+j) = 0.01;
   lb(6*T+1+1+j) = 0.001;
end

%%% Initial Guess %%%
%%%%%%%%%%%%%%%%%%%%%
%Note: Good initial guess is critical.
%It is recommended to load and use whichever previous scenario
%comes closest to the scenario one is looking to analyze.
%For larger parameter changes, it may also help to run intermediate
%scenarios at intermediate values to produce better initial guesses.
%Listed below are only the scenarios from Table 4, but the full set of
%prior scenario and sensitivity results are included further below to be
%imported here as well if applicable.

%Table 4 Scenarios:
if distortionary==0
     load('LS_Benchmark','x')                           %Lump-Sum Taxation (first-best) - Table 4 Scenario 6
     x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==1 && isempty(tao_l_fix)==1 && isempty(tao_E_fix)==1 &&  no_interm_tax==0 && tao_k_const==0 && tao_l_const==0 
    load('Opt_Benchmark','x')                           %Optimized Distortionary - Table 4 Scenario 5
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==1 && isempty(tao_l_fix)==1 && isempty(tao_E_fix)==0 && mean(tao_E_fix(1:T_tao_E_fix))<0.02 &&  no_interm_tax==0 && tao_k_const==0 && tao_l_const==0 
    load('Opt_MAC_0_EWedge_0','x')                      %Optimized Distortionary & No Carbon Price -  Table 4 Scenario 2
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==0 && isempty(tao_l_fix)==1 && isempty(tao_E_fix)==1 &&  no_interm_tax==0 && tao_l_const==1
    load('BAU_taoK3457_taoLconst','x')                  %BAU Income Taxes - taoK=34.57%, taoL variable but constant over time - Table 4 Scenario 4b
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==0 && isempty(tao_l_fix)==1 && isempty(tao_E_fix)==0 && mean(tao_E_fix(1:T_tao_E_fix))>0.01 &&  no_interm_tax==0 && tao_l_const==1
    load('MAC_LS_EWedge_LS_taoK3457_taoLconst','x')     %BAU Income Taxes (taoK=34.57%, taoL constant) & `Wrong' Carbon pricing - Table 4 Scenario 3b
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==0 && isempty(tao_l_fix)==1 && isempty(tao_E_fix)==0 && mean(tao_E_fix(1:T_tao_E_fix))<0.02 &&  no_interm_tax==0 && tao_l_const==1
    load('MAC_0_EWedge_0_taoK3457_taoLconst','x')       %BAU Income Taxes (taoK=34.57%, taoL constant) & No Carbon Price - Table 4 Scenario 1a
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==1 && isempty(tao_l_fix)==0 && isempty(tao_E_fix)==1 &&  no_interm_tax==0 && tao_k_const==1
    load('BAU_taoKconst_taoL3825','x')                  %BAU Income Taxes - taoK variable but constant over time, taoL=38.25% - Table 4 Scenario 4a
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==1 && isempty(tao_l_fix)==0 && isempty(tao_E_fix)==0 && mean(tao_E_fix(1:T_tao_E_fix))>0.01 &&  no_interm_tax==0 && tao_k_const==1
    load('MAC_LS_EWedge_LS_taoKconst_taoL3825','x')     %BAU Income Taxes (taoK constant, taoL=38.25%) & `Wrong' Carbon pricing - Table 4 Scenario 3a
    x0 = x;
elseif distortionary==1 && isempty(tao_k_fix)==1 && isempty(tao_l_fix)==0 && isempty(tao_E_fix)==0 && mean(tao_E_fix(1:T_tao_E_fix))<0.02 &&  no_interm_tax==0 && tao_k_const==1
    load('MAC_0_EWedge_0_taoKconst_taoL3825','x')       %BAU Income Taxes (taoK constant, taoL=38.25%) & No Carbon Price - Table 4 Scenario 1b
    x0 = x;
end
%Table 5 Scenarios:
if isempty(tao_k_fix)==0 && tao_l_const==0;
    load('BAU_taoK3457_taoLflex','x')                    %Table 5 Scenario 4b'
    x0 = x;
elseif isempty(tao_l_fix)==0 && tao_k_const==0;
    load('BAU_taoKflex_taoL3825','x')                    %Table 5 Scenario 4a'
    x0 = x;
end

%%% Test Constraints and Objective Function %%%
[f] = COMET_Objective(x0,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip)
[c,ceq] = COMET_Constraints(x0,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,b3,B0,energy_wedge_fix)


%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

 options = optimset('algorithm','sqp','Tolfun',1e-9,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(@(x)COMET_Objective(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip),x0,[],[],[],[],lb,ub,@(x)COMET_Constraints(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,b3,B0,energy_wedge_fix),options);

x0 = x;
options = optimset('algorithm','active-set','Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(@(x)COMET_Objective(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip),x0,[],[],[],[],lb,ub,@(x)COMET_Constraints(x,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,b3,B0,energy_wedge_fix),options);

%Note: Check exit flags to ensure convergence. If program stops for other reasons
%(e.g., max iterations exceeded), it needs to be re-started, either from
%the current guess, or, in case of failure, from a revised guess.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 4: Comptute Implementing Policies and Outcomes        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 C = x(1:T);
 L = x(T+1:2*T);
 E = x(2*T+1:3*T);
 pi1_l = x(3*T+1:4*T);
 K1t = x(4*T+1:5*T);
 ECleanPct = x(5*T+1:6*T);
 K2t = x(6*T+1+1:7*T+1);
 sT = x(6*T+1);
 
%%% Compute Temperature Change %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TempE = E;
for i = 1:1:T;
    TempE(i) = (1-ECleanPct(i))*TempE(i);
    TempE(i) = TempE(i)+ELand(i);
end
for i = (T+1):1:(T+periods);                        
   TempE(i) = (E(T)*(1-ECleanPct(T)))*((1+gXt(T))^(i));       %Assume balanced growth path after period T
   TempE(i) = TempE(i)+ELand(i);
end
    Zt = ones(T+periods,1); 
    Zt(1) = phi23*X0+phi33*Zt0;
    Xt = ones(T+periods,1); 
    Xt(1) = phi12*S_t0+phi22*X0+phi32*Zt0;
    St = ones(T+periods,1); 
    St(1) = ((E0*10)+ELand0)+phi11*S_t0+phi21*X0;
    for i = 2:1:T+periods;
        Zt(i) = phi23*Xt(i-1)+phi33*Zt(i-1);
        Xt(i) = phi12*St(i-1)+phi22*Xt(i-1)+phi32*Zt(i-1);
        St(i) = TempE(i-1)+phi11*St(i-1)+phi21*Xt(i-1);
    end
    Qt = ones(T+periods,1); 
    Qt(1) = Qt0*(1-ksi4)+ksi4*TC0;
    Ft = ones(T+periods,1); 
    Ft(1) = (eta*((log((((St(1)+St(2))/2)+0.000001)/Sbar))/log(2)))+Fx(1); 
    TC = ones(T+periods,1);
    TC(1) = TC0+ksi1*(Ft(1)-ksi2*TC0-ksi3*(TC0-Qt0));
    for i = 2:1:T-1+periods;
        Qt(i) = Qt(i-1)*(1-ksi4)+ksi4*TC(i-1);
        Ft(i) = (eta*((log((((St(i)+St(i+1))/2)+0.000001)/Sbar))/log(2)))+Fx(i);
        TC(i) = TC(i-1)+ksi1*(Ft(i)-ksi2*TC(i-1)-ksi3*(TC(i-1)-Qt(i-1)));
    end
    m = (T-1)+periods;
        Ft(m+1) = (eta*(log(((St(m+1)+0.000001)/Sbar))/log(2)))+Fx((m+1));
        TC(m+1) = TC(m)+ksi1*(Ft(m+1)-ksi2*TC(m)-ksi3*(TC(m)-Qt(m)));
        Qt(m+1) = Qt(m)*(1-ksi4)+ksi4*TC(m);   
 
%%% Compute Output %%%
%%%%%%%%%%%%%%%%%%%%%%
Yt = zeros(T,1);
for j=0:1:T-1;   
    Yt(j+1) = (((1+theta1*(TC(1+j))^2)^(-1))*(Z(1+j))*(((x(T+1+j)*x(3*T+1+j)*N(j+1))^(1-alpha-v))*(((x(2*T+1+j)))^(v))*((N(1+j)*10000*x(4*T+1+j))^alpha))); %Billions of 2005 PPP Int. Dollars per decade
end
 
%%% Compute Factor Prices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MPE = zeros(T,1);   %Marginal Product of Energy per metric ton carbon-equivalent
for i = 1:1:T;
    MPE(i) = (v*Yt(i))/(E(i));                         
end
MPL = zeros(T,1);   %Marginal Product of Labor per worker-decade
for i = 1:1:T;
    MPL(i) = ((1-alpha-v)*Yt(i))/(N(i)*pi1_l(i)*L(i));  
end
MPK = ones(T,1);    %Marginal Product of Capital
for i = 1:1:T;
  MPK(i) = (alpha*Yt(i))/(K1t(i)*10000*N(i));
end
 
%%% Taxes and Wedges %%%
%%%%%%%%%%%%%%%%%%%%%%%%
LaborTax = zeros(T,1);
MRStime = ones(T,1);
CapitalTax = ones(T,1);
CapitalTax(1) = tao_k_0;
CarbonTax_Wedge = zeros(T,1);
Uct = zeros(T,1);
Ult = zeros(T,1);
for i = 1:1:T;
    CarbonTax_Wedge(i) = (MPE(i)-((MPL(i)*(1-pi1_l(i))*N(i)*L(i))/(alphaE*E(i))));
     Uct(i) = (x(i)^(-sigma))*((1-phi_labor*x(T+i))^(gamma_labor*(1-sigma)));
     Ult(i) = (x(i)^(1-sigma))*(gamma_labor)*(-1)*(phi_labor)*((1-phi_labor*x(T+i))^(gamma_labor*(1-sigma)-1));
    LaborTax(i) = 1+((Ult(i)/Uct(i))/((MPL(i)/10000)));      
end
for i = 2:1:T;
  MRStime(i) = Uct(i-1)/(beta*Uct(i));
  CapitalTax(i) = 1 - ((MRStime(i)-1)/((MPK(i))-(1-(1-delta)^10)));
end
 MAC = zeros(T,1);
 at = zeros(T,1);
 b0t = zeros(T,1);
 b1t = zeros(T,1);
 denom = zeros(T,1);
 Eclean = zeros(T,1);
 for i = 1:1:T;
     Eclean(i) = x(5*T+i)*E(i);
     at(i) = a1+a4*log(i);
     b0t(i) = a2+a3*log(i);
     b1t(i) = b1+b2*log(i);
     denom(i) = 1+at(i)*exp(b0t(i)-b1t(i)*(Eclean(i)^b3));
     MAC(i) = ((gamma*Pc(i)*1000)*((denom(i)^(-2))*b1t(i)*b3*(Eclean(i)^b3)*(denom(i)-1))+(gamma*Pc(i)*1000)*denom(i)^(-1));
 end

%%% Marginal Cost of Public Funds %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   lambda1raw = lambda.ineqnonlin(T+1:2*T);
   lambda1 = lambda.ineqnonlin(T+1:2*T);
   for i = 1:1:T;
       lambda1(i) = lambda1(i)/((beta^(i-1))*N(i));
   end
   cons_eq = zeros(T,1);
   for i = 1:1:T;
       cons_eq(i) = 1/(N(i)*10000);
   end
   MCF_num = zeros(T,1);
   for i = 1:1:T;
       MCF_num(i) = (lambda1(i)/Uct(i))*(1/cons_eq(i));
       
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Table 4 Summary Metrics %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%Average labor tax:
avgLab = mean(LaborTax(2:T-1))

%Average capital tax:
avgk = mean(CapitalTax(2:T-1))

%Average MCF 2025-2255
avgMCF = mean(MCF_num(2:T-1))   

%Optimal Carbon Tax 2015, 2025, 2035:
CarbonTax = MAC(1:3)

%Optimal Carbon Tax Adjustment vis-a-vis First-Best:
load('MAC_LS','MAC_LS')
 ratio = zeros(T,1);
 for i = 1:1:T;
     ratio(i) = MAC(i)/MAC_LS(i);
 end
 OptvsFirstBestAdj = 1-mean(ratio(1:10))

%Maximum Temperature Change:
maxTC = max(TC)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Table 5 Decompositions: Utility Damages, Production Damages, Fiscal Constraint Interactions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: The code references Lagrange multipliers under the assumption that whichever tax is NOT fixed
% is also variable over time. Otherwise, the code would have to be adjusted.

if isempty(tao_k_fix)==0 && tao_l_const==0
    CapPsi = lambda.eqnonlin(1:T-1);
  elseif isempty(tao_l_fix)==0 && tao_k_const==0
    CapLam = lambda.eqnonlin(1:T);
end

     
%%% Marginal Tempereature Change dT/dE(1) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Told = TC;
Eold = TempE;
dTdE = zeros(T+periods,T+periods);
Tnew = zeros(T+periods,T+periods);
pert = 0.0001;  
for n = 1:1:T+periods;
    TempE = Eold;
    TempE(n) = TempE(n)+pert;
    Zt = ones(T+periods,1); 
    Zt(1) = phi23*X0+phi33*Zt0;
    Xt = ones(T+periods,1); 
    Xt(1) = phi12*S_t0+phi22*X0+phi32*Zt0;
    St = ones(T+periods,1); 
    St(1) = ((E0*10)+ELand0)+phi11*S_t0+phi21*X0;
    for i = 2:1:T+periods;
        Zt(i) = phi23*Xt(i-1)+phi33*Zt(i-1);
        Xt(i) = phi12*St(i-1)+phi22*Xt(i-1)+phi32*Zt(i-1);
        St(i) = TempE(i-1)+phi11*St(i-1)+phi21*Xt(i-1);
    end
    Qt = ones(T+periods,1); 
    Qt(1) = Qt0*(1-ksi4)+ksi4*TC0;
    Ft = ones(T+periods,1); 
    Ft(1) = (eta*((log((((St(1)+St(2))/2)+0.000001)/Sbar))/log(2)))+Fx(1); 
    TC = ones(T+periods,1);
    TC(1) = TC0+ksi1*(Ft(1)-ksi2*TC0-ksi3*(TC0-Qt0));
    for i = 2:1:T-1+periods;
        Qt(i) = Qt(i-1)*(1-ksi4)+ksi4*TC(i-1);
        Ft(i) = (eta*((log((((St(i)+St(i+1))/2)+0.000001)/Sbar))/log(2)))+Fx(i);
        TC(i) = TC(i-1)+ksi1*(Ft(i)-ksi2*TC(i-1)-ksi3*(TC(i-1)-Qt(i-1)));
    end
    m = (T-1)+periods;
        Ft(m+1) = (eta*(log(((St(m+1)+0.000001)/Sbar))/log(2)))+Fx((m+1));
        TC(m+1) = TC(m)+ksi1*(Ft(m+1)-ksi2*TC(m)-ksi3*(TC(m)-Qt(m)));
        Qt(m+1) = Qt(m)*(1-ksi4)+ksi4*TC(m);   
        Tnew((1:T+periods),n) = TC;
end
for j = 1:1:T+periods;
    for m = 1:1:T+periods;
    dTdE(m,j) = Tnew(m,j)-Told(m);
    end
end
dTdE = dTdE*(1/pert);
dTdE(1);

%Marginal Utility Damages%
%%%%%%%%%%%%%%%%%%%%%%%%%%
TC = Told;
MED_U = zeros(T,1);
MED_Upub = zeros(T,1);
for k = 1:1:T;
    UTt = zeros(T+periods-k+1,1);
    for i = 1:1:T+periods-k+1;
        UTt(i) = (beta^(i-1))*N(k-1+i)*(-1)*((1+alpha0*(TC(k-1+i)^alpha1))^(((-1)*(1-sigma))-1))*(alpha0*2*TC(k-1+i));
    end
    temp = zeros(T+periods+1-k+1,1);
    for j = 1:1:T+periods-k+1;
        temp(j) = UTt(j)*dTdE(j+k-1,k);
    end
UTBGP = (beta^(T+periods-1-k+1))*N(T+periods)*(1/(1-beta))*UTt(T+periods-k+1);
temp(T+periods+1) = UTBGP*dTdE(T+periods,k);
MED_U(k) = sum(temp);                    %Units: (aggregate utility/billions of people) / 1 billion metric tons of E
MED_U(k) = MED_U(k)*(1/(Uct(k)*N(k)));   %Units:  (D c in $10k/billions) / (1 billion metric tons of E)
MED_U(k) = MED_U(k)*10000*N(k);          %Pigouvian Level, Units: $/mtC
MED_Upub(k) = MED_U(k)*(1/MCF_num(k));   %Adjusted Level, Units: $/mtC
end

%Marginal Output Damages%
%%%%%%%%%%%%%%%%%%%%%%%%%%
pi1_k = x(5*T)/(x(5*T)+x(7*T+1));       %Period T share of capital in final goods production
Kfut = ones(periods,1);                 %Continuation aggregate capital stock, bil. int. 2005 PPP dollars
Kfut(1) = sT*(Yt(T)-Gct(T)+(1-Delta)*(N(T)*10000*(x(5*T)+x(7*T+1))));
L(T:1:T+periods) = L(T);
Yfut = zeros(periods,1);                %Continuation output, bil. int. 2005 PPP dollars
for i = 1:1:(periods-1);
  Yfut(i) = ((((1+theta1*(TC(T+i))^2)^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^i))^(1-alpha-v))*((E(T)*((1+gXt(T))^(i)))^(v))*((pi1_k*Kfut(i))^alpha))));
  Kfut(i+1) = sT*(Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i));
  C(T+i) = ((Yfut(i)-Gct(T+i)+(1-Delta)*Kfut(i))*(1-sT))/(N(T)*10000);
end
Yfut(periods) =  (((1+theta1*(TC(T+periods)^2))^(-1))*Z(T)*(((x(4*T)*x(2*T)*N(T)*((1+gXt(T))^periods))^(1-alpha-v))*(((E(T))*((1+gXt(T))^periods))^(v))*((pi1_k*Kfut(periods))^alpha)));
C(T+periods) = (Yfut(periods)-Gct(T+periods)-(Delta+gXt(T))*Kfut(periods))/(N(T)*10000);
for i = 1:1:periods;
    Yt(T+i) = Yfut(i);
    Uct(T+i) = (C(T+i)^(-sigma))*((1-phi_labor*L(T+i))^(gamma_labor*(1-sigma)));
end
YTt = zeros(1,T+periods);
for i =1:1:T+periods;
    YTt(i) = (-1)*((1+theta1*(TC(i))^2)^(-1))*(theta1*2*TC(i))*Yt(i);
end
MDLCon = zeros(T,1);    %Fiscal Constraint Value
MED_Y = zeros(T,1);
MED_Ypub = zeros(T,1);
tao_Pigou = zeros(T,1);
tao_M = zeros(T,1);
tao_M1 = zeros(T,1);
for k = 1:1:T;
    tempcon = zeros(T,1);
    temp = zeros(T+periods+1,1);
    temppub1 = zeros(T,1);
    temppub2 = zeros(periods+1,1);
    for j = 1:1:T+periods;
            temp(j) = (beta^(j-1-k+1))*Uct(j)*YTt(j)*dTdE(j,k);        
    end
    for j = 1:1:T;
        temppub1(j) = lambda1raw(j)*YTt(j)*dTdE(j,k);
        if isempty(tao_l_fix)==0 && tao_k_const==0
            tempcon(j) = CapLam(j)*(1-alpha-v)*(YTt(j)*(1/(N(j)*pi1_l(j)*L(j))))*(1/1000000000)*(100000)*(1-tao_l_fix)*dTdE(j,k);
        end
    end
      for j = 2:1:T;
       if isempty(tao_k_fix)==0 && tao_l_const==0
         tempcon(j) = CapPsi(j-1)*(1/beta)*(alpha)*(YTt(j)*(1/(N(j)*K1t(j))))*(1/1000000000)*(100000)*(1-tao_k_fix)*dTdE(j,k);
       end
      end
    for j = T+1:1:T+periods;
        temppub2(j) = temp(j);
    end
MED_Y(k) = sum(temp);
MED_Y(k) = MED_Y(k)*(1/(Uct(k)));
MED_Ypub(k) = (sum(temppub1)*(1/lambda1raw(k)));
 MDLCon(k) = sum(tempcon)*(1/lambda1raw(k));
MED_Ypub(k) = MED_Ypub(k)+(sum(temppub2)*(1/Uct(k)));
YTBGP = (beta^(T+periods-1-k+1))*((Uct(T+periods)*N(T))/(Uct(k)*N(k)))*(YTt(T+periods))*dTdE(T+periods,k)*(1/(1-beta*(1+gXt(T))^(1-sigma)));
MED_Y(k) = MED_Y(k)+YTBGP;
MED_Ypub(k) = MED_Ypub(k)+YTBGP;
tao_Pigou(k) = (MED_Y(k)+MED_U(k))*(-1);
tao_M1(k) = (MED_Ypub(k)+MED_Upub(k))*(-1);
 tao_M(k) = (MED_Ypub(k)+MED_Upub(k)+MDLCon(k))*(-1);
end

%%% Table 5 Results %%%
%%%%%%%%%%%%%%%%%%%%%%%
if (isempty(tao_k_fix)==0 && tao_l_const==0) || (isempty(tao_l_fix)==0 && tao_k_const==0)
MED_U_Pigou = MED_U(2)      %Pigouvian Utility Damages
MED_U_Opt = MED_Upub(2)     %Adjusted Internalization of Utility Damages
MED_Y_Pigou = MED_Y(2)      %Pigouvian Output Damages
MED_Y_Opt = MED_Ypub(2)     %Adjusted Internalization of Output Damages
FCI = MDLCon(2)             %Fiscal Constraint Interaction Value
avgMCF21stC = mean(MCF_num(1:10))   
end

%%%NOTE: Save Results if Desired%%%
%Naming convention for saved results referenced below:
%Optimal allocation:
% save('ScenarioName','x')
%Other output used for figures, e.g.:
% CarbonTax_ScenarioName = MAC;
% save('CarbonTax_ScenarioName','CarbonTax_ScenarioName')
% tao_Pigou_ScenarioName = tao_Pigou;
% save('tao_Pigou_ScenarioName','tao_Pigou_ScenarioName')
% etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Section 5. Welfare Calculations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcomp = x;

%Note 1a: The parameter values as set in Section 1 must match the ones
%of the scenarios one wishes to compare! For example, if one wishes to compare
%welfare of two model runs where Sigma=2 (Table 7), one needs to run
%the parameter code in Section 1 with Sigma=2. Given the interdependencies
%of the parameters (e.g., labor preferences parameter depends on Sigma),
%all of the parameters must be set simultaneously to ensure consistency.

%Note 1b: The exception to this is for the welfare calculation comparing
%the first-best (lump-sum taxation, Scenario 6) allocation to the second-best
%(optimized distortionary, Scenario 5) in Table 4, which uses the
%parameters from the second-best setting and compare the allocations.
if distortionary==0
    warning('Parameters assume first-best setting, welfare comparison to 1a inconsistent!')
end 
%To replicate this calculation, evaluate the parameters in second-best
%("distortionary=1") and then load the lump-sum scenario results here:
% load('LS_Benchmark','x')     %Table 4 Scenario 6
%  xcomp = x;

%Note 2: The setup below assumes that the model run in question yields higher
%welfare (lower program objective function value f) than the benchmark BAU 
%comparison scenarios 1a. In order to calculate negative welfare effects,
%the code would need to be adjusted (with negative consumption iterations
%and the relevant inequalities reversed).

%Lump-Sum Change in 2015 Consumption:%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load ONE appropriate benchmark comparison scenario:
 load('MAC_0_EWedge_0_taoK3457_taoLconst','x')                   %Table 4 Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK4085_taoLconst_Frisch2','x')     %Table 7, Frisch=2,  Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK394_taoLconst_sig2','x')         %Table 7, Sigma=2,   Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK332_taoLconst_sig11','x')        %Table 7, Sigma=1.1, Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK296_taoconst_g09','x')           %Table 7, G-10%,     Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK38_taoLconst_g11','x')           %Table 7, G+10%,     Scenario 1a
x_wrongE = x;
x_wrongE_vary = x;
[f] = COMET_Objective(xcomp,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip);
target_val = f;
pc = 0.0001;
diff = 1;
while diff > 0.0000001
    pc = pc+0.0001;
    x_wrongE(1) = x_wrongE_vary(1)+pc;
   [f] = COMET_Objective(x_wrongE,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip);
    diff = f-target_val;
end
pc = pc*10000;  %Dollars
WelfareAmt = pc*N(1)    %$2005 Billions:

%Permanent Percentage Change in Consumption%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multip = 0;
%Load appropriate benchmark comparison scenario:
load('MAC_0_EWedge_0_taoK3457_taoLconst','x')                   %Table 4 Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK4085_taoLconst_Frisch2','x')     %Table 7, Frisch=2,  Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK394_taoLconst_sig2','x')         %Table 7, Sigma=2,   Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK332_taoLconst_sig11','x')        %Table 7, Sigma=1.1, Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK296_taoconst_g09','x')           %Table 7, G-10%,     Scenario 1a
% load('MAC_0_EWedge_0_BAU_taoK38_taoLconst_g11','x')           %Table 7, G+10%,     Scenario 1a
x_wrongE = x;
[f] = COMET_Objective(xcomp,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip);
target_val = f;
pct = 0.00000001;
diff = 1;
while diff > 0.00001
    pct = pct+0.000001;
    multip = pct;
   [f] = COMET_Objective(x_wrongE,T,periods,N,K0,A_E,alphaE,S_t0,theta1,Z,alpha,v,Gct,Pc,Delta,phi_labor,gamma_labor,alpha0,alpha1,tao_l_fix,tao_l_const,tao_k_fix,tao_k_const,T_tao_E_fix,tao_E_fix,no_interm_tax,phi23,X0,phi33,Zt0,phi12,phi22,phi32,phi11,phi21,Qt0,ksi4,TC0,eta,Sbar,Fx,ksi1,ksi2,ksi3,E0,ELand0,ELand,gXt,beta,sigma,trans_share,tao_k_0,delta,K1,G,tao_l_0,distortionary,gamma,a1,a2,a3,a4,b1,b2,multip);
	diff = f-target_val;
end
%Percent:
WelfarePct = pct*100
multip = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Section 6: Saved Results and Figures     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Optimal Allocations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Benchmark Calibration, Table 4:
%  load('MAC_0_EWedge_0_taoKconst_taoL3825','x')                 %Table 4 Scenario 1b 
%  load('Opt_MAC_0_EWedge_0','x')                                %Table 4 Scenario 2 
%  load('MAC_LS_EWedge_LS_taoKconst_taoL3825','x')               %Table 4 Scenario 3a 
%  load('BAU_taoKconst_taoL3825','x')                            %Table 4 Scenario 4a 
%  load('MAC_LS_EWedge_LS_taoK3457_taoLconst','x')               %Table 4 Scenario 3b 
%  load('BAU_taoK3457_taoLconst','x')                            %Table 4 Scenario 4b 
%  load('Opt_Benchmark','x')                                     %Table 4 Scenario 5 
%  load('LS_Benchmark','x')                                      %Table 4 Scenario 6

%%%Benchmark Calibration, Table 5:
%  load('BAU_taoK3457_taoLflex','x')                             %Table 5 Scenario 4b'
%  load('BAU_taoKflex_taoL3825','x')                             %Table 5 Scenario 4a'
 
%%%Alternative Ouput Damage Shares, Figure 4:
%  load('Opt_allU','x')                                          %0%   Output Damages, Scen. 5,  Figure 4
%  load('Opt_25pY','x')                                          %25%  Output Damages, Scen. 5,  Figure 4
%  load('Opt_50pY','x')                                          %50%  Output Damages, Scen. 5,  Figure 4
%  load('Opt_allY','x')                                          %100% Output Damages, Scen. 5,  Figure 4
%  load('BAU_taoKconst_taoL3825_allU','x')                       %0%   Output Damages, Scen. 4a, Figure 4
%  load('BAU_taoKconst_taoL3825_25pY','x')                       %25%  Output Damages, Scen. 4a, Figure 4
%  load('BAU_taoKconst_taoL3825_50pY','x')                       %50%  Output Damages, Scen. 4a, Figure 4
%  load('BAU_taoKconst_taoL3825_allY','x')                       %100% Output Damages, Scen. 4a, Figure 4
%  load('BAU_taoK3457_taoLconst_allU','x')                       %0%   Output Damages, Scen. 4b, Figure 4
%  load('BAU_taoK3457_taoLconst_25pY','x')                       %25%  Output Damages, Scen. 4b, Figure 4
%  load('BAU_taoK3457_taoLconst_50pY','x')                       %50%  Output Damages, Scen. 4b, Figure 4
%  load('BAU_taoK3457_taoLconst_allY','x')                       %100% Output Damages, Scen. 4b, Figure 4
 
%%%Frisch elasticity=2, Tables 6 & 7%:
%  load('MAC_LS_EWedge_LS_BAU_taoK4085_taoLconst_Frisch2','x')  %Table 7, Frisch=2, Scenario 3b
%  load('BAU_taoK4085_taoLconst_Frisch2','x')                   %Table 7, Frisch=2, Scenario 4b
%  load('MAC_LS_EWedge_LS_BAU_taoKconst_taoL39_Frisch2','x')    %Table 7, Frisch=2, Scenario 3a
%  load('BAU_taoKconst_taoL39_Frisch2','x')                     %Table 7, Frisch=2, Scenario 4a
%  load('Opt_Frisch2','x')                                      %Table 7, Frisch=2, Scenario 5
%  load('LS_Frisch2','x')                                       %Table 7, Frisch=2, Scenario 6
 
%%%Sigma=2, Tables 6 & 7%:
%  load('MAC_LS_EWedge_LS_BAU_taoK394_taoLconst_Sig2','x')      %Table 7, Sigma=2, Scenario 3b
%  load('BAU_taoK394_taoLconst_Sig2','x')                       %Table 7, Sigma=2, Scenario 4b
%  load('MAC_LS_EWedge_LS_BAU_taoKconst_taoL4075_Sig2','x')     %Table 7, Sigma=2, Scenario 3a
%  load('BAU_taoKconst_taoL4075_sig2','x')                      %Table 7, Sigma=2, Scenario 4a
%  load('Opt_sig2','x')                                         %Table 7, Sigma=2, Scenario 5
%  load('LS_sig2','x')                                          %Table 7, Sigma=2, Scenario 6

%%%Sigma=1.1, Tables 6 & 77:
%  load('MAC_LS_EWedge_LS_BAU_taoK332_taoLconst_sig11','x')     %Table 7, Sigma=1.1, Scenario 3b
%  load('BAU_taoK332_taoLconst_sig11','x')                      %Table 7, Sigma=1.1, Scenario 4b
%  load('MAC_LS_EWedge_LS_BAU_taoKconst_taoL35_sig11','x')      %Table 7, Sigma=1.1, Scenario 3a
%  load('BAU_taoKconst_taoL35_sig11','x')                       %Table 7, Sigma=1.1, Scenario 4a
%  load('Opt_sig11','x')                                        %Table 7, Sigma=1.1, Scenario 5

%%%G=G(.9), Tables 6 & 7:
%  load('MAC_LS_EWedge_LS_BAU_taoK296_taoLconst_g09','x')       %Table 7, G-10%, Scenario 3b
%  load('BAU_taoK296_taoconst_g09','x')                         %Table 7, G-10%, Scenario 4b
%  load('MAC_LS_EWedge_LS_BAU_taoKconst_taoL32_g09','x')        %Table 7, G-10%2, Scenario 3a
%  load('BAU_taoKconst_taoL32_g09','x')                         %Table 7, G-10%, Scenario 4a
%  load('Opt_g09','x')                                          %Table 7, G-10%, Scenario 5
%  load('LS_g09','x')                                           %Table 7, G-10%, Scenario 6
 
%%%G=G(1.1), Tables 6 & 7:
%  load('MAC_LS_EWedge_LS_BAU_taoK38_taoLconst_g11','x')        %Table 7, G+10%, Scenario 3b
%  load('BAU_taoK38_taoLconst_g11','x')                         %Table 7, G+10%, Scenario 4b
%  load('MAC_LS_EWedge_LS_BAU_taoKconst_taoL45_g11','x')        %Table 7, G+10%, Scenario 3a
%  load('BAU_taoKconst_taoL45_g11','x')                         %Table 7, G+10%, Scenario 4a
%  load('Opt_g11','x')                                          %Table 7, G+10%, Scenario 5
%  load('LS_g11','x')                                           %Table 7, G+10%, Scenario 6

%%%Benchmark Calibration, Further Scenarios, Online Appendix Table A5:
%   load('BAU_taoKflex_taoL3825','x')                           %Online Appendix Table A5, Scenario 4a'
%   load('BAU_taoKconst_taoL3825_noTaoInt','x')                 %Online Appendix Table A5, Scenario 4a''
%   load('BAU_taoK3457_taoLflex','x')                           %Online Appendix Table A5, Scenario 4b'
%   load('BAU_taoK3457_taoLconst_noTaoInt','x')                 %Online Appendix Table A5, Scenario 4b''
 
%%%Benchmark Calibration, Further Scenarios, Online Appendix Table A6:
%  load('LS_abt12','x')                                         %Online Appendix Table A6, +20% Abatement Costs, Scenario 6
%  load('Opt_abt12','x')                                        %Online Appendix Table A6, +20% Abatement Costs, Scenario 5
%  load('BAU_taoKconst_taoL3775_abt12','x')                     %Online Appendix Table A6, +20% Abatement Costs, Scenario 4a
%  load('BAU_taoK3457_taoLconst_abt12','x')                     %Online Appendix Table A6, +20% Abatement Costs, Scenario 4b
%  load('LS_abt08','x')                                         %Online Appendix Table A6, -20% Abatement Costs, Scenario 6
%  load('Opt_abt08','x')                                        %Online Appendix Table A6, -20% Abatement Costs, Scenario 5
%  load('BAU_taoKconst_taoL3775_abt08','x')                     %Online Appendix Table A6, -20% Abatement Costs, Scenario 4a
%  load('BAU_taoK3457_taoLconst_abt08','x')                     %Online Appendix Table A6, -20% Abatement Costs, Scenario 4b

%%%No-Recalibration Scenario, Online Appendix Figures A2-A4
%  load('LS_NoRecalib','x')                                     %Online Appendix Figures A2,A3,A4


%%% Figures %%%
%%%%%%%%%%%%%%%

%%% Figure 1: Optimal Carbon Taxes Across Fiscal Scenarios %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MAC_LS','MAC_LS')
load('MAC_Opt','MAC_Opt')
load('MAC_taoK3457_taoLconst','MAC_taoK3457_taoLconst')
load('MAC_taoKconst_taoL3825','MAC_taoKconst_taoL3825')
num = 11;
plot(y(1:num),MAC_LS(1:num),'k*-',y(1:num),MAC_Opt(1:num),'rd-',y(1:num),MAC_taoK3457_taoLconst(1:num),':pm',y(1:num),MAC_taoKconst_taoL3825(1:num),'bO-')
h1leg = legend('First-Best (Scen. 6):     Lump-Sum Taxes,       MCF=1.0','Optimized (Scen. 5):  \tau_K variable, \tau_L variable, MCF~1.1','BAU (Scen. 4b):        \tau_K =34.6%,  \tau_L variable, MCF~1.1','BAU (Scen. 4a):        \tau_K variable, \tau_L = 38.3%, MCF~1.4','Location','NorthWest');
xlabel('Year','FontSize',12)
ylabel('Carbon Tax ($/mtC)','FontSize',12)
title('Optimal Carbon Tax Paths Across Fiscal Scenarios','FontSize',13)


%%% Figure 2: Optimal vs. Pigouvian Carbon Tax Levels %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('tao_Pigou_LS','tao_Pigou_LS')
load('tao_Pigou_Opt','tao_Pigou_Opt')
load('tao_Pigou_taoK3457_taoLconst','tao_Pigou_taoK3457_taoLconst')
load('tao_Pigou_taoKconst_taoL3825','tao_Pigou_taoKconst_taoL3825')
num = 10;
plot(y(1:num),tao_Pigou_LS(1:num),'k*:',y(1:num),MAC_LS(1:num),'k*-',y(1:num),tao_Pigou_taoKconst_taoL3825(1:num),'bO:',y(1:num),MAC_taoKconst_taoL3825(1:num),'o-b')
h1leg = legend('\tau{Pigou} in First-Best (Scen. 6, MCF=1.0)','\tau_E*       in First-Best (Scen. 6, MCF=1.0)','\tau{Pigou} with Dist. Taxes (Scen. 4a, MCF~1.4)','\tau_E*       with Dist. Taxes (Scen. 4a, MCF~1.4)','Location','Northwest');
xlabel('Year','FontSize',12)
ylabel('Carbon Price ($/mtC)','FontSize',12)
title('Optimal vs. Pigouvian Carbon Tax Levels','FontSize',13)


%%% Figure 3: Optimal vs. Pigouvian Carbon Tax Ratios %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i = 1:1:T;
     ratio0(i) = MAC_LS(i)/tao_Pigou_LS(i);
     ratio1(i) = MAC_Opt(i)/tao_Pigou_Opt(i);
     ratio3(i) = MAC_taoKconst_taoL3825(i)/tao_Pigou_taoKconst_taoL3825(i);
     ratio4(i) = MAC_taoK3457_taoLconst(i)/tao_Pigou_taoK3457_taoLconst(i);
 end

num = 10;
plot(y(1:num),ratio0(1:num),'k*-',y(1:num),ratio1(1:num),'rd-',y(1:num),ratio4(1:num),'-pm',y(1:num),ratio3(1:num),'ob-')
h1leg = legend('First-Best (Scen. 6):    Lump-Sum Taxes,        MCF=1.00','Optimized (Scen. 5):  \tau_K variable, \tau_L variable, MCF~1.1','BAU (Scen. 4b):        \tau_K =34.6%,  \tau_L variable, MCF~1.1','BAU (Scen. 4a):        \tau_K variable, \tau_L = 38.3%, MCF~1.4','Location','SouthEast');
xlabel('Year','FontSize',12)
ylabel('Optimal vs. Pigouvian Tax Ratio','FontSize',12)
title('Optimal vs. Pigouvian Carbon Tax Ratios','FontSize',13)

%%% Figure 4: Optimal vs. Pigouvian Carbon Tax across Ouput Share Values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MAC_Opt_allY','MAC_Opt_allY')
load('MAC_Opt_50pY','MAC_Opt_50pY')
load('MAC_Opt_25pY','MAC_Opt_25pY')
load('MAC_Opt_allU','MAC_Opt_allU')
load('tao_Pigou_Opt_allY','tao_Pigou_Opt_allY')
load('tao_Pigou_Opt_50pY','tao_Pigou_Opt_50pY')
load('tao_Pigou_Opt_25pY','tao_Pigou_Opt_25pY')
load('tao_Pigou_Opt_allU','tao_Pigou_Opt_allU')
load('MAC_taoKconst_taoL3825_allY','MAC_taoKconst_taoL3825_allY')
load('MAC_taoKconst_taoL3825_50pY','MAC_taoKconst_taoL3825_50pY')
load('MAC_taoKconst_taoL3825_25pY','MAC_taoKconst_taoL3825_25pY')
load('MAC_taoKconst_taoL3825_allU','MAC_taoKconst_taoL3825_allU')
load('tao_Pigou_taoKconst_taoL3825_allY','tao_Pigou_taoKconst_taoL3825_allY')
load('tao_Pigou_taoKconst_taoL3825_50pY','tao_Pigou_taoKconst_taoL3825_50pY')
load('tao_Pigou_taoKconst_taoL3825_25pY','tao_Pigou_taoKconst_taoL3825_25pY')
load('tao_Pigou_taoKconst_taoL3825_allU','tao_Pigou_taoKconst_taoL3825_allU')
load('MAC_taoK3457_taoLconst_allY','MAC_taoK3457_taoLconst_allY')
load('MAC_taoK3457_taoLconst_50pY','MAC_taoK3457_taoLconst_50pY')
load('MAC_taoK3457_taoLconst_25pY','MAC_taoK3457_taoLconst_25pY')
load('MAC_taoK3457_taoLconst_allU','MAC_taoK3457_taoLconst_allU')
load('tao_Pigou_taoK3457_taoLconst_allY','tao_Pigou_taoK3457_taoLconst_allY')
load('tao_Pigou_taoK3457_taoLconst_50pY','tao_Pigou_taoK3457_taoLconst_50pY')
load('tao_Pigou_taoK3457_taoLconst_25pY','tao_Pigou_taoK3457_taoLconst_25pY')
load('tao_Pigou_taoK3457_taoLconst_allU','tao_Pigou_taoK3457_taoLconst_allU')

yr = 2;
OptVec = zeros(5,1);
OptVec(1) = (MAC_Opt_allY(yr)/tao_Pigou_Opt_allY(yr));
OptVec(2) = (MAC_Opt(yr)/tao_Pigou_Opt(yr));
OptVec(3) = (MAC_Opt_50pY(yr)/tao_Pigou_Opt_50pY(yr));
OptVec(4) = (MAC_Opt_25pY(yr)/tao_Pigou_Opt_25pY(yr));
OptVec(5) = (MAC_Opt_allU(yr)/tao_Pigou_Opt_allU(yr));
taoK3457_taoLconstVec = zeros(5,1);
taoK3457_taoLconstVec(1) = (MAC_taoK3457_taoLconst_allY(yr)/tao_Pigou_taoK3457_taoLconst_allY(yr));
taoK3457_taoLconstVec(2) = (MAC_taoK3457_taoLconst(yr)/tao_Pigou_taoK3457_taoLconst(yr));
taoK3457_taoLconstVec(3) = (MAC_taoK3457_taoLconst_50pY(yr)/tao_Pigou_taoK3457_taoLconst_50pY(yr));
taoK3457_taoLconstVec(4) = (MAC_taoK3457_taoLconst_25pY(yr)/tao_Pigou_taoK3457_taoLconst_25pY(yr));
taoK3457_taoLconstVec(5) = (MAC_taoK3457_taoLconst_allU(yr)/tao_Pigou_taoK3457_taoLconst_allU(yr));
taoKconst_taoL3825Vec = zeros(5,1);
taoKconst_taoL3825Vec(1) = (MAC_taoKconst_taoL3825_allY(yr)/tao_Pigou_taoKconst_taoL3825_allY(yr));
taoKconst_taoL3825Vec(2) = (MAC_taoKconst_taoL3825(yr)/tao_Pigou_taoKconst_taoL3825(yr));
taoKconst_taoL3825Vec(3) = (MAC_taoKconst_taoL3825_50pY(yr)/tao_Pigou_taoKconst_taoL3825_50pY(yr));
taoKconst_taoL3825Vec(4) = (MAC_taoKconst_taoL3825_25pY(yr)/tao_Pigou_taoKconst_taoL3825_25pY(yr));
taoKconst_taoL3825Vec(5) = (MAC_taoKconst_taoL3825_allU(yr)/tao_Pigou_taoKconst_taoL3825_allU(yr));
udamvec = zeros(5,1);
udamvec(1) = 100;
udamvec(2) = 75;
udamvec(3) = 50;
udamvec(4) = 25;
udamvec(5) = 0;

plot(udamvec,OptVec,'-k*',udamvec,taoK3457_taoLconstVec,'-ob',udamvec,taoKconst_taoL3825Vec,'-rd')
h1leg = legend('Optimized (Scenario 5): \tau_K variable, \tau_L variable, MCF~1.1','BAU (Scenario 4b):        \tau_K =34.6%,  \tau_L variable, MCF~1.1','BAU (Scenario 4a):        \tau_K variable, \tau_L = 38.3%, MCF~1.4','Location','Best');
xlabel('Production Damages Share (%)','FontSize',12)
ylabel('Optimal vs. Pigouvian Tax Ratio','FontSize',12)
title('Optimal vs. Pigouvian Carbon Tax Ratios (2025)','FontSize',12)


%%% Figure A1: Optimal Capital Income Tax Path %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('CapitalTax_Opt','CapitalTax_Opt')
num = 15;
plot(y(1:num),CapitalTax_Opt(1:num),'-*k')
h1leg = legend('Capital Income Tax (Scen. 5)','Location','East');
xlabel('Year','FontSize',12)
ylabel('Capital Income Tax','FontSize',12)
title('Optimal Capital Income Tax Path','FontSize',13)

%%% Figure A2: Re-Calibration and Labor Supply %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Lt_LS_noRecalib','Lt_LS_noRecalib')
load('Lt_LS','Lt_LS')
load('Lt_Opt','Lt_Opt')
num = 10;
plot(y(1:num),Lt_LS(1:num),'k*-',y(1:num),Lt_LS_noRecalib(1:num),':pb',y(1:num),Lt_Opt(1:num),'rd-')
h1leg = legend('First-Best (Scen. 6), Benchmark','First-Best (Scen. 6), Not Recalibrated','Optimized Dist. (Scen. 5), Benchmark','Location','Best')
xlabel('Year','FontSize',12)
ylabel('Equilibrium Labor Supply','FontSize',12)
title('Re-Calibration and Labor Supply','FontSize',13)

%%% Figure A3: Re-Calibration and Output %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Yt_LS_noRecalib','Yt_LS_noRecalib')
load('Yt_LS','Yt_LS')
load('Yt_Opt','Yt_Opt')
num = 10;
plot(y(1:num),Yt_LS(1:num),'k*-',y(1:num),Yt_LS_noRecalib(1:num),':pb',y(1:num),Yt_Opt(1:num),'rd-')
h1leg = legend('First-Best (Scen. 6), Benchmark','First-Best (Scen. 6), Not Recalibrated','Optimized Dist. (Scen. 5), Benchmark','Location','Best')
xlabel('Year','FontSize',12)
ylabel('Output ($bil./decade)','FontSize',12)
title('Re-Calibration and Output','FontSize',13)

%%% Figure A4: Re-Calibration and Optimal Carbon Taxes%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MAC_LS_noRecalib','MAC_LS_noRecalib')
load('MAC_LS','MAC_LS')
load('MAC_Opt','MAC_Opt')
num = 10;
plot(y(1:num),MAC_LS(1:num),'k*-',y(1:num),MAC_LS_noRecalib(1:num),':pb',y(1:num),MAC_Opt(1:num),'rd-')
h1leg = legend('First-Best (Scen. 6), Benchmark','First-Best (Scen. 6), Not Recalibrated','Optimized Dist. (Scen. 5), Benchmark','Location','Best')
xlabel('Year','FontSize',12)
ylabel('Carbon Tax ($/mtC)','FontSize',12)
title('Re-Calibration and Optimal Carbon Taxes','FontSize',13)
