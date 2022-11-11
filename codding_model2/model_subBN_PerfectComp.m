%% Model with simple subjective basic needs
% eppsilons =1 => no low skill in sustainable sector
% eppsilonn =0 => no high skill in unsustainable sector


% no innovation 
% no monopolistic competition
% no labour sector choice

% provide all variables as a function of ps, pn
clear, clc

cd('/home/sonja/Documents/projects/Overconsumption/codding_model')
% variables
syms ps pn omegas lambdaa cbarh cbarl An As psii real % all other variables are functions of ps, pn


% vector of variables for which to solve model
symm.variables=[ps]; 
% parameters
symm.params= [ omegas, psii, cbarh, cbarl, lambdaa, As, An];

%% model

pn=1;
wh= 0.25/psii*As*ps^2;
wl= 0.25/psii*An*pn^2;

csh= wh^2/(ps-1)^2*(2*omegas-1) -cbarh/(ps-1)-ps/(ps-1)^2*(1-omegas)+1/(ps-1)^2*(omegas);
csl= wl^2/(ps-1)^2*(2*omegas-1) -cbarl/(ps-1)-ps/(ps-1)^2*(1-omegas)+1/(ps-1)^2*(omegas);
lh= wh/(ps-1)*(2*omegas-1);

c=0;
c=c+1;
f(c)= lh*As*ps-(lambdaa*csh+(1-lambdaa)*csl);

S=solve(f==0,ps, 'ReturnConditions', true, 'Real', true)

syms x z real
eq= An^2 + 16*psii^2*x + 16*cbarl*psii^2*x + As^2*lambdaa*x^4 + 4*As^2*psii*x^3 + 2*An^2*lambdaa*omegas + 16*cbarl*lambdaa*psii^2 + 16*cbarh*lambdaa*psii^2*x + 8*As^2*omegas*psii*x^4 -(An^2*lambdaa + 2*An^2*omegas + 16*cbarl*psii^2 + 16*omegas*psii^2 + 16*omegas*psii^2*x + 4*As^2*psii*x^4 + 16*cbarh*lambdaa*psii^2 + 16*cbarl*lambdaa*psii^2*x + 2*As^2*lambdaa*omegas*x^4 + 8*As^2*omegas*psii*x^3);
S=solve(eq==0,x, 'ReturnConditions', true, 'Real', true)
simplify(S.x)
sol_om_cbar= subs(S.x, [As, An, lambdaa, psii], [1,2,0.5, 0.25]);

d=subs(sol_om_cbar, [omegas, cbarl, cbarh], [0.2, 0.5, 1])

m= z^4 - 2*z^3 - (31*z)/6 - 5/6;
S=solve(m==0,z, 'ReturnConditions', true, 'Real', true)

%% alternative model 
% with U=log(C), C=cn^(1-omegas)*cs^omegas;
clear, clc

cd('/home/sonja/Documents/projects/Overconsumption/codding_model')
% variables
syms cnh csh ps pn omegas lambdaa cbarh cbarl An As psii wh lh real % all other variables are functions of ps, pn


eq(1)= cnh^2 +(cbarh-csh-1/(1-ps)*omegas/csh)*cnh +(1-omegas)*ps/(1-ps);
eq(2) = lh-(wh/csh-cnh-csh-cbarh)*wh;
eq(3) = csh*ps-wh*lh+cnh;

S=solve(eq==0, [cnh, csh, lh]);

% choose positive solution for cnh
innd=eval(subs(S, [csh, omegas, ps, cbarh], [1, 0.5, 1.2, 0.4]))>0;
cnh=S(innd);

eq= cnh^2*