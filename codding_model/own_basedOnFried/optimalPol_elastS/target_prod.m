function f= target_prod(x, pf, Ems, params, list, T)
% function to find output side given emission target
% and fossil price

% params
read_in_params;

% read in vars
pn=exp(x((find(list.targprod=='pn')-1)*T+1:find(list.targprod=='pn')*T));
pg=exp(x((find(list.targprod=='pg')-1)*T+1:find(list.targprod=='pg')*T));

F       = (Ems'+deltaa)./omegaa;
G       = F.*(pf./pg).^(eppse);
E       = F.*(1+(pf./pg).^(eppse-1)).^(eppse/(eppse-1));
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); 
N       = (pee./pn).^eppsy.*(1-deltay)/deltay.*E; % optimality final good
Yprod   = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production
Ydem    = pf.*F+pg.*G+pn.*N; 

% to solve
f(1:T*1)     = Yprod-Ydem; % pg ensures output clearing!
f(T*1+1:T*2) = pn - ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

end