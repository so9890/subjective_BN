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
Sp.Lf=exp(sol3(2));
Sp.Lg=exp(sol3(1));
Sp.Af=(1+vf)*init201519(list.init=='Af0');
Sp.Ag=(1+vg)*init201519(list.init=='Ag0');
Sp.h=Sp.Lg+Sp.Lf;
Sp.F=Sp.Af*Sp.Lf;
Sp.G=Sp.Ag*Sp.Lg;
Sp.Y=(Sp.F)^(eppsy)*(Sp.G)^(1-eppsy);
Sp.C=Sp.Y;

%% Optimal solution 
x0=log([0.4,0.4]);

modFF = @(x)easy_opt(x, params, list, pol,  init201519, indic);
options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt','Display', 'Iter');%, );%, );%,  );
[sol2, fval, exitf] = fsolve(modFF, x0, options);
options = optimoptions('fsolve', 'TolFun', 10e-10, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Display', 'Iter');%, );%, );%,  );
[sol3, fval, exitf] = fsolve(modFF, sol2, options);

% save solution
read_in_params;
read_in_pol;

Opt.Af=(1+vf)*init201519(list.init=='Af0');
Opt.Ag=(1+vg)*init201519(list.init=='Ag0');
Opt.pg=exp(sol3(1));
Opt.Lg=exp(sol3(2));

Opt.G = Opt.Ag*Opt.Lg;
Opt.w = Opt.pg*Opt.Ag;
Opt.pf = Opt.w/((1-tauf)*Opt.Af);
Opt.Lf = Opt.pg/Opt.pf*Opt.Ag/Opt.Af*eppsy/(1-eppsy)*Opt.Lg;
Opt.F = Opt.Af*Opt.Lf;
Opt.Y = (Opt.F)^(eppsy)*(Opt.G)^(1-eppsy);
Opt.h = Opt.Lf+Opt.Lg;
Opt.lambdaa = (Opt.w*Opt.h+tauf*Opt.pf*Opt.F)/((Opt.w*Opt.h)^(1-taul));
Opt.hsup =  (Opt.lambdaa^(1-thetaa)*(1-taul)*Opt.w^((1-taul)*(1-thetaa))/chii)^(1/(sigmaa+taul+thetaa*(1-taul)));
Opt.C  = Opt.lambdaa*(Opt.w*Opt.h)^(1-taul);
Opt.Cd= Opt.w*Opt.h+tauf*Opt.pf*Opt.F;
% test goods market clearing
