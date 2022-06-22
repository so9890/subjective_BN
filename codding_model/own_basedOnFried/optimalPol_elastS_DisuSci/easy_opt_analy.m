read_in_params;
read_in_pol;
 
% Lg = exp(x(1)); 
% Lf = exp(x(2));
% lambdaa =exp(x(3));

% exogenous variables: tauf, taul
% endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

Af=(1+vf)*init201519(list.init=='Af0');
Ag=(1+vg)*init201519(list.init=='Ag0');

% analytic solution
pg = Af*(1-tauf)/(eppsy*Ag+(1-eppsy)*Af*(1-tauf));
pf= Ag/((1-tauf)*Af)*pg;
w= Ag*pg;
Lg = tauf*((Ag^3*eppsy*Af)...
    /(Ag^(3-eppsy)*eppsy*Af*(eppsy*Ag+(1-eppsy)*Af*(1-tauf))...
        - Ag*(1-tauf)*(Ag^2*eppsy+(1-tauf)*(1-eppsy)*Af^2)));
Lf = pg/pf*Ag/Af*eppsy/(1-eppsy)*Lg;
h = Lf+Lg;
F = Lf*Af;
lambdaa = (w*h+tauf*pf*F)/((w*h)^(1-tauf));
C  = lambdaa*(w*h)^(1-taul);
G = Lg*Ag;
Y = C; 

% to test labour market clearing
hsup  = (lambdaa^(1-thetaa)*(1-taul)*w^((1-taul)*(1-thetaa))/chii)^(1/(sigmaa+taul+thetaa*(1-taul)));
Ysup =(F)^(eppsy)*(G)^(1-eppsy);

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

