% simple model to compare social planner and 
% competitive equilibrium

function f=easy_opt(x, params, list, indic, init201519, pol)

read_in_params;
read_in_pol;

Lg = exp(x(1)); 
Lf = exp(x(2));
lambdaa =exp(x(3));

% auxiliary stuff
Af=(1+vf)*init201519(list.init=='Af0');
Ag=(1+vg)*init201519(list.init=='Ag0');
h=Lg+Lf;
F=Af*Lf;
G=Ag*Lg;
Y=(F)^(eppsy)*(G)^(1-eppsy);
%- utility
if indic.util==1
    thetaa=2;
    Ucon=(C^(1-thetaa)-1)/(1-thetaa);
else
    Ucon=log(C);
end
Uc=C^(-thetaa);
Uh=-chii*h^sigmaa;

%-externality
Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);

%- derivatives production 
pg =F^eppsy*(1-eppsy)*G^(-eppsy);
pf =G^(1-eppsy)*(eppsy)*F^(eppsy-1);
wg= pg*Ag;
wf= pf*Af*(1-tauf);


% model
f(1) = Uc*w+Uh; %HH optimality
f(2) = C-lambdaa*(w*h)^(1-taul); % budget
% gov budget
f(3) = lambdaa -(w*h+tauf*pf*F)/((w*h)^(1-tauf));
% production 
f(4)=1-(pf^(1-eppsy)+pg^(1-eppsy))^(1/(1-eppsy));
% labour market
f(5)=wg-wf;

end
