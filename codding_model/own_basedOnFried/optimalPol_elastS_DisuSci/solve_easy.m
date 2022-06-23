% to solve easy sp and op problems

%% sp solution 
%  iin=load('init_techgap.mat');
%  init201519_co = iin.initcount;
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
Sp.pg= eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-eppsy)/eppsy*Sp.Lf*Af/Ag/Sp.Lg)^eppsy;
% %- utility
if indic.util==1
    thetaa=2;
    Sp.Ucon=(Sp.C^(1-thetaa)-1)/(1-thetaa);
else
    Sp.Ucon=log(Sp.C);
end
Sp.Ext = -weightext*(omegaa*Sp.F)^extexpp;
Sp.Ulab = -chii*Sp.h^(1+sigmaa)/(1+sigmaa);

Sp.SWF = Sp.Ucon+Sp.Ulab+indic.extern*Sp.Ext;

% prigou from a households problem:
% Uf = -weightext*extexpp*(omegaa)^extexpp*Sp.F.^(extexpp-1);
% Sp.pimarket= -Uf./Sp.C^(-thetaa)/;
%% Competitive equilibrium
indic.taxsch=1; %==1 then uses linear tax schedule

x0=log([Sp.pg,Sp.Lg]);
tauf=Sp.pigou;
taul=0;
taus=0;
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
clear LF
LF.Af=(1+vf)*init201519(list.init=='Af0');
LF.Ag=(1+vg)*init201519(list.init=='Ag0');
LF.pg=exp(sol3(1));
LF.Lg=exp(sol3(2));

LF.G = LF.Ag*LF.Lg;
LF.w = LF.pg*LF.Ag;
LF.pf = LF.w/((1-tauf)*LF.Af);
LF.Lf = LF.pg/LF.pf*LF.Ag/LF.Af*eppsy/(1-eppsy)*LF.Lg;
LF.F = LF.Af*LF.Lf;
LF.Y = (LF.F)^(eppsy)*(LF.G)^(1-eppsy);
LF.h = LF.Lf+LF.Lg;
LF.lambdaa = (LF.w*LF.h+tauf*LF.pf*LF.F)/((LF.w*LF.h)^(1-taul));
LF.s =LF.Lf/LF.h;
if indic.taxsch==0
    LF.hsup =  (LF.lambdaa^(1-thetaa)*(1-taul)*LF.w^((1-taul)*(1-thetaa))/chii)^(1/(sigmaa+taul+thetaa*(1-taul)));
    LF.taul =  1-LF.h^(thetaa+sigmaa)*chii*(LF.w+LF.tauf*LF.pf*LF.Af*LF.s)^(thetaa-1);

else % linear tax schedule
    LF.hsup = (((LF.w+tauf*LF.pf*LF.F/LF.h)^(-thetaa)*LF.w*(1-taul))/(chii))^(1/(sigmaa+thetaa));
    LF.taul =  1-LF.h^(thetaa+sigmaa)*chii*(LF.w+tauf*LF.pf*LF.Af*LF.s)^(thetaa)/LF.w;
end
LF.C  = LF.lambdaa*(LF.w*LF.h)^(1-taul);
LF.Cd= LF.w*LF.h+tauf*LF.pf*LF.F;

% %- utility
if indic.util==1
    thetaa=2;
    LF.Ucon=(LF.C^(1-thetaa)-1)/(1-thetaa);
else
    LF.Ucon=log(LF.C);
end
LF.Ext = -weightext*(omegaa*LF.F)^extexpp;
LF.Ulab = -chii*LF.h^(1+sigmaa)/(1+sigmaa);

LF.SWF = LF.Ucon+LF.Ulab+indic.extern*LF.Ext;
%% optimal policy
indic.notaul=0;
indic.taxsch=0;
clear Opt
if indic.notaul==0
    x0=log([Sp.s,Sp.h]);
else
    x0=log([Sp.s]);
end
modFF = @(x)easy_opt(x, params, list,  init201519, indic);
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
    if indic.taxsch==0
        Opt.h = ((Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
    else
        % with linear tax: 
        Opt.h = ((Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(-thetaa)*Opt.w/chii)^(1/(sigmaa+thetaa));
    end
end
% labour market clearing: 
Opt.Lg = (1-Opt.s)*Opt.h; 
Opt.Lf = Opt.s*Opt.h;
% production
Opt.G = Opt.Ag*Opt.Lg;
Opt.F = Opt.Af*Opt.Lf;
% income tax/ gov budget

if indic.taxsch==0
    Opt.taul= 1-Opt.h^(thetaa+sigmaa)*chii*(Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(thetaa-1);
    Opt.lambdaa = (Opt.w*Opt.h+Opt.tauf*Opt.pf*Opt.F)/((Opt.w*Opt.h)^(1-Opt.taul));
elseif indic.taxsch==1
    Opt.taul = 1-Opt.h^(thetaa+sigmaa)*chii*(Opt.w+Opt.tauf*Opt.pf*Opt.Af*Opt.s)^(thetaa)/Opt.w;
    Opt.T = Opt.w*Opt.taul*Opt.h+Opt.pf*Opt.tauf*Opt.F;
end
% good market clearing and final production 
Opt.C = (Opt.F)^(eppsy)*(Opt.G)^(1-eppsy); 

% %- utility
if indic.util==1
    thetaa=2;
    Opt.Ucon=(Opt.C^(1-thetaa)-1)/(1-thetaa);
else
    Opt.Ucon=log(Opt.C);
end
Opt.Ext = -weightext*(omegaa*Opt.F)^extexpp;
Opt.Ulab = -chii*Opt.h^(1+sigmaa)/(1+sigmaa);

Opt.SWF = Opt.Ucon+Opt.Ulab+indic.extern*Opt.Ext;
if indic.notaul==1
    Optnotaul=Opt;
else
    Opttaul=Opt;
end