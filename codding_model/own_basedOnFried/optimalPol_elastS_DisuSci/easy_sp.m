% simple model to compare social planner and 
% competitive equilibrium

function f=easy_sp(x, params, list, indic, init201519)

read_in_params;
Lg = exp(x(1)); 
Lf = exp(x(2));

% auxiliary stuff
Af=(1+vf)*init201519(list.init=='Af0');
Ag=(1+vg)*init201519(list.init=='Ag0');
h=Lg+Lf;
F=Af*Lf;
G=Ag*Lg;
Y=(F)^(eppsy)*(G)^(1-eppsy);
C=Y;

Uc=C^(-thetaa);
Uh=-chii*h^sigmaa;

%-externality

Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);

%- derivatives production 
dYdLg= F^eppsy*(1-eppsy)*Ag^(1-eppsy)*Lg^(-eppsy);
dYdLf= G^(1-eppsy)*(eppsy)*Af^(eppsy)*Lf^(eppsy-1);
dFdLf = Af;
% model
f(1) = Uc*dYdLg+Uh;
f(2) = Uc*dYdLf+Uh+indic.extern*Uf*dFdLf;
end
