% to solve easy sp and op problems
addpath('tools')
% to save results of loop
keySet={'<1', 'log', 'Bop'};
valueSet= repmat({struct([])},1,length(keySet));
resultsTHETA=containers.Map(keySet, valueSet);
indic.taxsch=0; %==0 uses nonlinear, no lump sum trans
                %==1 then uses linear tax schedule with lump sum transfers
                %==2 linear tax without lump sum transfers
indic.notaul=0; % relevant for optimal solution

for tth=0:1 
    indic.util=tth;
for bbp=0:indic.util % if indic.util==0 then loop is independent of BOP
    indic.Bop=bbp;

    
%% sp solution 
clear Sp
x0=log([0.4,0.4]);

modFF = @(x)easy_sp(x, params, list, indic, init201519);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','Display', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Display', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution
read_in_params;
clear Sp;
Sp.Lf=exp(sol3(2));
Sp.Lg=exp(sol3(1));
Sp.Af=(1+vf)*init201519(list.init=='Af0');
Sp.Ag=(1+vg)*init201519(list.init=='Ag0');
Sp.h=Sp.Lg+Sp.Lf;
Sp.F=Sp.Af*Sp.Lf;
Sp.G=Sp.Ag*Sp.Lg;
Sp.Y=(Sp.F)^(eppsy)*(Sp.G)^(1-eppsy);
Sp.C=Sp.Y;
Sp.pigou =1-(1-eppsy)/eppsy*Sp.Lf/Sp.Lg;
Sp.s = Sp.Lf/Sp.h; % share of labour in fossil sector Lf/h
Sp.w= (Sp.Af*Sp.s)^eppsy*(Sp.Ag*(1-Sp.s))^(1-eppsy);
Sp.pg= eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-eppsy)/eppsy*Sp.Lf*Sp.Af/Sp.Ag/Sp.Lg)^eppsy;
% %- utility
if indic.util==1
    Sp.Ucon=(Sp.C^(1-thetaa)-1)/(1-thetaa);
else
    Sp.Ucon=log(Sp.C);
end
Sp.Ext = -weightext*(omegaa*Sp.F)^extexpp;
Sp.dEdF = -weightext*extexpp*(omegaa*Sp.F)^extexpp/Sp.F;
Sp.Ulab = -chii*Sp.h^(1+sigmaa)/(1+sigmaa);

Sp.SWF = Sp.Ucon+Sp.Ulab+indic.extern*Sp.Ext;
Sp.Ul = -chii*Sp.h^sigmaa;
% test analytic derivation of Sp.h
% Sp.htest= ((Sp.w^(1-thetaa)+Sp.dEdF*Sp.Af*Sp.s*Sp.h^thetaa)/chii)^(1/(sigmaa+thetaa));
Sp.htest = (Sp.w^(1-thetaa)/chii*(1-eppsy)/(1-Sp.s))^(1/(sigmaa+thetaa));

Sp.Ultest = -Sp.C^(-thetaa)*Sp.Y/Sp.G*(1-eppsy)*Sp.Ag; % marginal disutility of labour 
%% Competitive equilibrium

clear LF
x0=log([Sp.pg,Sp.Lg]);
LF.tauf=Sp.pigou;
LF.taul=1-((1-eppsy)/(1-Sp.s))^thetaa; % this tax implements efficient hours worked given tauf=SP.pigou
taus=0;
tauf=LF.tauf;
taul=LF.taul;
lambdaa=1;
pol=eval(symms.pol);
modFF = @(x)easy_lf(x, params, list, pol,  init201519, indic);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','Display', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Display', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution
read_in_params;
read_in_pol;

LF.Af=(1+vf)*init201519(list.init=='Af0');
LF.Ag=(1+vg)*init201519(list.init=='Ag0');
LF.pg=exp(sol3(1));
LF.Lg=exp(sol3(2));

LF.G = LF.Ag*LF.Lg;
LF.w = LF.pg*LF.Ag;
LF.pf = LF.w/((1-LF.tauf)*LF.Af);
LF.Lf = LF.pg/LF.pf*LF.Ag/LF.Af*eppsy/(1-eppsy)*LF.Lg;
LF.F = LF.Af*LF.Lf;
LF.Y = (LF.F)^(eppsy)*(LF.G)^(1-eppsy);
LF.h = LF.Lf+LF.Lg;
LF.lambdaa = (LF.w*LF.h+LF.tauf*LF.pf*LF.F)/((LF.w*LF.h)^(1-LF.taul));
LF.s =LF.Lf/LF.h;
if indic.taxsch==0
    LF.hsup =  (LF.lambdaa^(1-thetaa)*(1-LF.taul)*LF.w^((1-LF.taul)*(1-thetaa))/chii)^(1/(sigmaa+LF.taul+thetaa*(1-LF.taul)));
    LF.taul =  1-LF.h^(thetaa+sigmaa)*chii*(LF.w+LF.tauf*LF.pf*LF.Af*LF.s)^(thetaa-1);
    LF.C  = LF.lambdaa*(LF.w*LF.h)^(1-LF.taul);

elseif indic.taxsch==1 % linear tax schedule
    LF.hsup = (((LF.w+LF.tauf*LF.pf*LF.F/LF.h)^(-thetaa)*LF.w*(1-LF.taul))/(chii))^(1/(sigmaa+thetaa));
    LF.taul =  1-LF.h^(thetaa+sigmaa)*chii*(LF.w+LF.tauf*LF.pf*LF.Af*LF.s)^(thetaa)/LF.w;
    LF.C= LF.w*LF.h+LF.tauf*LF.pf*LF.F;
    LF.T=LF.tauf*LF.pf*LF.F;
elseif indic.taxsch==2 % linear tax schedule but no transfers
    LF.hsup = ((LF.w)^(1-thetaa)*(1-LF.taul)/chii)^(1/(thetaa+sigmaa));
    LF.Gov = LF.tauf*LF.pf*LF.F;
    LF.C = LF.w*LF.h; % transfers from labour tax
    LF.aggdemand =LF.Gov+LF.C;
end

% %- utility
if indic.util==1
    LF.Ucon=(LF.C^(1-thetaa)-1)/(1-thetaa);
else
    LF.Ucon=log(LF.C);
end
LF.Ext = -weightext*(omegaa*LF.F)^extexpp;
LF.dEdF = LF.Ext*extexpp/(LF.F); % negative!
LF.Ulab = -chii*LF.h^(1+sigmaa)/(1+sigmaa);

LF.SWF = LF.Ucon+LF.Ulab+indic.extern*LF.Ext;

% social cost of carbon as: what would a household be willing to pay
LF.scc = -LF.dEdF*LF.C^thetaa/LF.pf;

% analytic derivation of transfers so that competitive eqbm = First Best
% only if utility = log
if indic.util==0
   LF.Tana=-LF.dEdF*LF.Af*LF.s*LF.h^thetaa/(1+LF.dEdF*LF.Af*LF.s*LF.h^thetaa)*LF.w*LF.h;
   LF.Tana_gov= LF.tauf/(1-LF.tauf)*LF.w*LF.s*LF.h;
   LF.taufcheckk=-LF.dEdF*LF.Af*LF.h/(1-LF.dEdF*LF.Af*LF.h*(1-LF.s));
end

%% optimal policy

clear Opt
if indic.notaul==0
    x0=log([Sp.s,Sp.h]);
else
    x0=log([Sp.s]);
end
% if indic.taxsch<=1
    modFF = @(x)easy_opt(x, params, list,  init201519, indic);
% else
%     modFF = @(x)easy_opt_Gov(x, params, list,  init201519, indic);
% end
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','Display', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Display', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution 
Opt.s = exp(sol3(1));
if indic.notaul==0
    Opt.h = exp(sol3(2));
end
% auxiliary
read_in_params;

Opt.Af=(1+vf)*init201519(list.init=='Af0');
Opt.Ag=(1+vg)*init201519(list.init=='Ag0');

Opt.tauf = 1-((1-eppsy)/eppsy)*Opt.s/(1-Opt.s); % tauf determines s
Opt.pg = eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-Opt.tauf)*Opt.Af/Opt.Ag)^eppsy;
Opt.w = Opt.pg*Opt.Ag;
Opt.pf = Opt.w/((1-Opt.tauf)*Opt.Af); 

if indic.notaul==1
    if indic.taxsch==0 % baseline model
        Opt.h = ((Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
    elseif indic.taxsch==1 % with linear tax and transfers
        Opt.h = ((Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(-thetaa)*Opt.w/chii)^(1/(sigmaa+thetaa));
    elseif indic.taxsch==2 % with linear tax no transfers
        Opt.h = (Opt.w^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
    end
end
% labour market clearing: 
Opt.Lg = (1-Opt.s)*Opt.h; 
Opt.Lf = Opt.s*Opt.h;
% production
Opt.G = Opt.Ag*Opt.Lg;
Opt.F = Opt.Af*Opt.Lf;
% income tax/ gov budget

if indic.taxsch==0 % non linear tax scheme 
    Opt.taul= 1-Opt.h^(thetaa+sigmaa)*chii*(Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(thetaa-1);
    Opt.lambdaa = (Opt.w*Opt.h+Opt.tauf*Opt.pf*Opt.F)/((Opt.w*Opt.h)^(1-Opt.taul));
elseif indic.taxsch==1 % linear tax scheme with transfers
    Opt.taul = 1-Opt.h^(thetaa+sigmaa)*chii*(Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(thetaa)/Opt.w;
    Opt.T = Opt.w*Opt.taul*Opt.h+Opt.pf*Opt.tauf*Opt.F;
elseif indic.taxsch ==2
    Opt.taul= 1-(Opt.h^(thetaa+sigmaa)*chii)/(Opt.w)^(1-thetaa);
    Opt.Gov = Opt.tauf*Opt.pf*Opt.F;
    Opt.T =Opt.w*Opt.h*Opt.taul;
    Opt.Cdem = Opt.w*Opt.h;
end
% good market clearing and final production 
Opt.Y = (Opt.F)^(eppsy)*(Opt.G)^(1-eppsy); 

if indic.taxsch<=1
    Opt.C=Opt.Y;
else
    Opt.C=Opt.Y-Opt.Gov;
end
% %- utility
if indic.util==1
    Opt.Ucon=(Opt.C^(1-thetaa)-1)/(1-thetaa);
else
    Opt.Ucon=log(Opt.C);
end
Opt.dUcondC = Opt.C^(-thetaa);
Opt.Ext = -weightext*(omegaa*Opt.F)^extexpp;
Opt.dEdF = -weightext*extexpp*(omegaa*Opt.F)^extexpp/Opt.F;
Opt.Ulab = -chii*Opt.h^(1+sigmaa)/(1+sigmaa);

Opt.SWF = Opt.Ucon+Opt.Ulab+indic.extern*Opt.Ext;
Opt.scc = -Opt.dEdF*Opt.C^thetaa/Opt.pf;

% check analytic version of taul if taxsch==0
if indic.taxsch==0
    Opt.taulcheck = 1-Opt.w*Opt.h/(Opt.Y);
    Opt.taultest= 1-((1-eppsy)/(1-Opt.s))^thetaa; % expression for taul to have optimal h in taxsche==0
end
if indic.notaul==1
    Optnotaul=Opt;
else
    Opttaul=Opt;
end

%check derivative Gov
    
    % derivatives
    dFdh = Opt.Af*Opt.s;
    dFds = Opt.Af*Opt.h;
    dCdh = (Opt.Af*Opt.s)^(eppsy)*(Opt.Ag*(1-Opt.s))^(1-eppsy);
    MPL  = dCdh;
    dCds = dCdh*Opt.h*(eppsy/Opt.s-(1-eppsy)/(1-Opt.s));
    dwdH = 0;
    dwds = (1-eppsy)*Opt.Af^eppsy*Opt.Ag^(1-eppsy)*eppsy*(Opt.s/(1-Opt.s))^(eppsy-1)/(1-Opt.s)^2;
    if indic.taxsch>1
        dGovdh=Opt.tauf*Opt.pf*dFdh;
        dCdh =dCdh-dGovdh;
        dtaufds = -(1-eppsy)/eppsy*(1/(1-Opt.s)^2);
        dpfds =-Opt.pf*(eppsy-1)/(1-Opt.tauf)*dtaufds;
        dGovds = Opt.pf*Opt.F*dtaufds+Opt.tauf*Opt.F*dpfds+Opt.tauf*Opt.pf*dFds;
        dCds = dCds-dGovds;
    % analytic solution to check
    dGovdscheck= Opt.pf*Opt.F*(-(1-eppsy)/eppsy/(1-Opt.s)^2 *(eppsy/Opt.s)+(eppsy-Opt.s)/((1-Opt.s)*eppsy*Opt.s)); %correct!
    taulcheck= (dGovds*Opt.s/Opt.h-dGovdh)/((1-eppsy)*Opt.Y/Opt.G*(-1)*Opt.Ag); 
    taulcheck2= 1+(Opt.h*dCdh-Opt.s*dCds)/(-Opt.h*Opt.w);
    taulcc=dwds*(Opt.s/Opt.w)-dwdH*(Opt.h/Opt.w); 
    % check tauf
    taufcheckk = Opt.scc+dGovds/(Opt.w*Opt.h)*(1-Opt.tauf); 
    taufcheck = Opt.scc+ Opt.tauf-(1-Opt.tauf)/Opt.w*dwds;
    tauff= 1-Opt.scc*Opt.w/dwds; 
    Uh = -chii*Opt.h^(sigmaa);
    d2YdG2=-(1-eppsy)*Opt.Y/Opt.G^2+(1-eppsy)^2*Opt.Y/Opt.G^2;
    taullcheckk= -Opt.h/Opt.w*Opt.Ag^2*d2YdG2;
    end
    
%% save results
%- create structure
st.Sp=Sp;
st.Opt=Opt;
st.LF=LF;
st.thetaa = thetaa;
%- save to map
if indic.util==0
    resultsTHETA('log')= st;
elseif indic.util==1
    if indic.Bop==0
        resultsTHETA('<1')= st;
    else
        resultsTHETA('Bop')= st;

    end
end
end
end

% - create table from results 
kk=keys(resultsTHETA);
Table=table(keys(resultsTHETA)',zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1));
Table.Properties.VariableNames={'Thetaa','FB hours', 'FB SWF', 'FB s', 'FB MPL','FB Pigou', 'Only tauf=pigou hours' , 'Only tauf=pigou SWF' , 'Only tauf=pigou s' , 'Only tauf=pigou scc', 'Only tauf=pigou wage', ...
                                       'Optimal hours', 'Optimal SWF','Optimal s', 'Optimal wage', 'Optimal taul', 'Optimal tauf', 'Optimal scc'};

%- only hours, and policy
TableH=table(keys(resultsTHETA)',zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1),zeros(length(keys(resultsTHETA)),1));
TableH.Properties.VariableNames={'Thetaa','FB hours', 'FB Pigou', 'CE hours',  'CE scc',  ...
                                       'Opt hours','Opt taul', 'Opt tauf', 'Opt scc'};

for i=1:3
    st=resultsTHETA(string(kk(i)));
    Table(i,2:end)={st.Sp.h, st.Sp.SWF,st.Sp.s, st.Sp.w, st.Sp.pigou, st.LF.h, st.LF.SWF, st.LF.s, st.LF.scc, st.LF.w, ...
                    st.Opt.h, st.Opt.SWF, st.Opt.s, st.Opt.w, Opt.taul, st.Opt.tauf, st.Opt.scc};
    TableH(i,2:end)={st.Sp.h, st.Sp.pigou, st.LF.h, st.LF.scc, ...
                    st.Opt.h, Opt.taul, st.Opt.tauf, st.Opt.scc};
end

table2latex(TableH)