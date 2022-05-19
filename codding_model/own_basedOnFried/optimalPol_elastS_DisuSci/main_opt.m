%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%COMET Model M-File%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Lint Barrage, Brown University
%%Version: August 2018
clear
%cd C:\Users\lintb\Desktop\COMET\ReStud
cd '/home/sonja/Documents/projects/subjective_BN/codding_model/own_basedOnFried/optimalPol_elastS_DisuSci'
%%
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
            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 1: Select Scenario        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 12;  % Direct optimization period time horizon: 2020-2080
         % one period = 5 years

lengthh = 5; % number of zears per period         
indic.util =0; % ==0 log utilit, otherwise as in Boppart
indic.target =0; % ==1 if uses emission target
indic.spillovers =1; % ==1 then there are positive spillover effects of scientists within sectors! 
indic.taus =0; % ==1 if taus is present in ramsey problem
indic.noskill =0; % == 1 if no skill calibration of model
indic.notaul=0;
indic.sep =0;% ==1 if uses models with separate markets for scientists
indic.dim=1; %==1 if uses diminishing returns to science
% savedind=indic;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 2: Parameters        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 1) sets direct parameters, 
%          2) calibrates model to indirect params.
if isfile(sprintf('params_spillovers%d.mat', indic.spillovers))
    fprintf('loading parameter values')
    load(sprintf('params_spillovers%d.mat', indic.spillovers),...
        'params', 'Sparams', 'polCALIB', 'init201014', 'init201519', 'list', 'symms', 'Ems', 'Sall', 'x0LF', 'MOM', 'indexx')
else
    fprintf('calibrating model')
    [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx]=get_params_Base( T, indic, lengthh);
    save(sprintf('params_spillovers%d', indic.spillovers))
end

%% 
%%%%% LIST for separate markets
    symms.sepchoice=symms.choice(list.choice~='S'&list.choice~='ws'& list.choice~='gammas');
    syms wsg wsn wsf gammasg gammasf gammasn real
    symms.sepchoice=[symms.sepchoice wsg wsn wsf gammasg gammasf gammasn];
    list.sepchoice=string(symms.sepchoice);
    
    symms.sepallvars=symms.allvars(list.allvars~='ws');
    symms.sepallvars=[symms.sepallvars wsg wsn wsf gammasg gammasf gammasn]; 
    list.sepallvars=string(symms.sepallvars);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 3: BAU Simulation        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this section I simulate the economy starting from 2015-2019
% order of variables in LF_SIM as in list.allvars
 params(list.params=='etaa')=1.1;
 Sparams.etaa=1.1;
for i=[0]
    indic.noskill=i;
if ~isfile(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')))
    [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, polCALIB, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
    helper.LF_SIM=LF_SIM;
%    helper=load(sprintf('LF_BAU_spillovers%d.mat', indic.spillovers));
    [LF_BAU]=solve_LF_VECT(T, list, polCALIB, params,symms, init201519, helper, indic);
    save(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')), 'LF_BAU', 'Sparams')
    clearvars LF_SIM pol FVAL
else
     fprintf('LF_BAU no skill %d exists', indic.noskill);
end
end
%% Competitive equilibrium with policy optimal without spillovers
% DOES NOT SOLVE WITH ETAA ==1
% for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
taus=0;
tauf=0;
taul=0;
lambdaa=1; % balances budget with tauf= taul=0
pol=eval(symms.pol);

  %  if indic.noskill==0
  for i=[0]
      indic.noskill=i;
  if ~isfile(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')))
%       indic.noskill=1;
        [LF_SIM, polLF, FVAL] =solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init201014, indexx, indic, Sall);
        helper.LF_SIM=LF_SIM;
        [LF_SIM]=solve_LF_VECT(T, list, pol, params,symms, init201519, helper, indic);
        save(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')),'LF_SIM', 'Sparams');
        clearvars LF_SIM helper
        else
     fprintf('LF_FB spillover %d no skill %d exists',indic.spillovers, indic.noskill)
  end
  end
%     end
%     if pol~=polLF
%         error('LF not solved under fb policy');
%     end

%%% Check swf value in LF
disc=repmat(Sparams.betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;
hhel= load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')));
sswfbau=vec_discount*hhel.LF_BAU( :, list.allvars=='SWF');
hhblf = load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d', indic.spillovers, indic.noskill, indic.sep));
sswf=vec_discount*hhblf.LF_SIM( :, list.allvars=='SWF');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 4: Sociel Planner allocation                             %%%
% Timing: starting from 2020-2025                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indic.noskill=0;
        for i =0
            indic.noskill=i;
            if ~isfile(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_gammac.mat', indic.spillovers, indic.noskill, indic.sep))
                indic.target=1;
                fprintf('solving Social planner solution with target, noskill%d', indic.noskill);
                SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems);
            else
              fprintf('Social planner solution with target, noskill%d exists', indic.noskill);
            end
            if ~isfile(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_gammac.mat', indic.spillovers, indic.noskill, indic.sep))
                indic.target=0;
                fprintf('solving Social planner solution without target, noskill%d', indic.noskill);
                SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems);
            else
              fprintf('Social planner solution without target, noskill%d exists', indic.noskill);
            end
        end
           
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 5: Solve for Optimal Allocation       %%%
% Timing: starting from 2020-2025 the gov. chooses      %%
% the optimal allocation                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indic.taus  = 0; % with ==0 no taus possible!
indic.notaul=0; % ==0 if labour income tax is available
indic.sep =0;
for i=1
     indic.noskill=i;
     if isfile(sprintf('OPT_target_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul))
         indic.target=1;
         [symms, list, opt_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems);
     else 
        fprintf('OPT solution with target, noskill%d exists', indic.noskill);
     end
    if isfile(sprintf('OPT_notarget_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul))
        indic.target=0;
        [symms, list, opt_all]= OPT_solve(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems);
    else 
       fprintf('OPT solution without target, noskill%d exists', indic.noskill);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Section 6: PLOTS       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plottsSP(list, symms, T);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find good initial starting point for target opt %%
%%% as a function of taul
taul= 0.5;
pf= 1;
exx = polExp(pf, params, polCALIB, list, taul, T, Ems, indexx);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Symbolic approach to solve Ramsey problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indic.target=1;
indic.noskill=0;
indic.sep=0;
indic.target=0; 
%- with etaa=1
params(list.params=='etaa')=1;
%1) get objective function 
if indic.sep==1
    [OB_RAM, list, symms, Ftarget]= model_ram_sep( list, params, T, init201519, indic, Ems, symms);
else
    [OB_RAM, list, symms, Ftarget]= model_ram( list, params, T, init201519, indic, Ems, symms);
end
%- x is a symbolic vector of choice variables! 

%2) take derivatives and write resulting equations as function
if indic.target==1
    [indexx, model, list]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_target_1905', list, indic, indexx);
else
    if indic.sep==0
        [indexx, model]=symmodel_eq(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_1905', list, indic, indexx);
    else
        [indexx, model]=symmodel_eq_sep(OB_RAM, symms.optALL, params,  Ftarget, 'Ram_Model_notarget_sep_1905', list, indic, indexx);
    end
end

%3) solve model using fsolve
if indic.sep==1
    RAM = solve_sym_sep(symms, list, Ftarget, indic);
else
    RAM = solve_sym(symms, list, Ftarget, indic);
end
%%
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