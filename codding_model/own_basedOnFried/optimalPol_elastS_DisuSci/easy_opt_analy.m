function [SolAn]=easy_opt_analy(params, list, indic, init201519, pol)
read_in_params;
read_in_pol;
 

% exogenous variables: tauf, taul
% endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

Af=(1+vf)*init201519(list.init=='Af0');
Ag=(1+vg)*init201519(list.init=='Ag0');

% analytic solution: only with thetaa=1
pg = (Af*(1-tauf)/Ag)^eppsy*eppsy^eppsy*(1-eppsy)^(1-eppsy);
pf= Ag/((1-tauf)*Af)*pg;
w= Ag*pg;
% variables in relation to h
Lgh= (1-eppsy)/((1-tauf)*eppsy+1-eppsy);
Lfh=(1-tauf)*eppsy/(1-eppsy)*Lgh;
hlpp=w+tauf*pf*Af*Lfh;
h = (((1-taul)/chii)*hlpp^(1-thetaa))^(1/(thetaa+sigmaa)); % from labour market clearing h supply and h demand 

Lg = Lgh*h;
Lf = Lfh*h;
F = Lf*Af;
G = Lg*Ag;
lambdaa = (w*h+tauf*pf*F)/((w*h)^(1-taul));
C  = lambdaa*(w*h)^(1-taul);

% to test labour market clearing
Ysup =(F)^(eppsy)*(G)^(1-eppsy);
s= Lf/h;
taufopt = 1-(1-eppsy)/eppsy*s/(1-s);
end
%- utility
% if indic.util==1
%     thetaa=2;
%     Ucon=(C^(1-thetaa)-1)/(1-thetaa);
% else
%     Ucon=log(C);
% end
% Uc=C^(-thetaa);
% Uh=-chii*h^sigmaa;
% 
% %-externality
% Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);

