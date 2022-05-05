function LF_SIM=aux_solutionLF_VECT(x, pol, list, symms,varrs, params, T)
read_in_params;
tauf=varrs(list.allvars=='tauf', :)';
taus=varrs(list.allvars=='taus', :)';
taul=varrs(list.allvars=='taul', :)';
hhf    = exp(x((find(list.test=='hhf')-1)*T+1:(find(list.test=='hhf'))*T));
 hhg    = exp(x((find(list.test=='hhg')-1)*T+1:(find(list.test=='hhg'))*T));
 hhn    = exp(x((find(list.test=='hhn')-1)*T+1:(find(list.test=='hhn'))*T));
 hln    = exp(x((find(list.test=='hln')-1)*T+1:(find(list.test=='hln'))*T));
 hlf    = exp(x((find(list.test=='hlf')-1)*T+1:(find(list.test=='hlf'))*T));
 hlg    = exp(x((find(list.test=='hlg')-1)*T+1:(find(list.test=='hlg'))*T));
 C      = exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T));
 F      = exp(x((find(list.test=='F')-1)*T+1:(find(list.test=='F'))*T));
 G      = exp(x((find(list.test=='G')-1)*T+1:(find(list.test=='G'))*T));
 Af     = exp(x((find(list.test=='Af')-1)*T+1:(find(list.test=='Af'))*T));
 Ag     = exp(x((find(list.test=='Ag')-1)*T+1:(find(list.test=='Ag'))*T));
 An     = exp(x((find(list.test=='An')-1)*T+1:(find(list.test=='An'))*T));
 hl     = upbarH./(1+exp(x((find(list.test=='HL')-1)*T+1:find(list.test=='HL')*T)));
 hh     = upbarH./(1+exp(x((find(list.test=='HH')-1)*T+1:(find(list.test=='HH'))*T)));
 sff     = exp(x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T));
 sg     = exp(x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T));
 sn     = exp(x((find(list.test=='sn')-1)*T+1:(find(list.test=='sn'))*T));
 S     = exp(x((find(list.test=='S')-1)*T+1:(find(list.test=='S'))*T));

 gammalh = x((find(list.test=='gammalh')-1)*T+1:(find(list.test=='gammalh'))*T).^2;
 gammall = x((find(list.test=='gammall')-1)*T+1:(find(list.test=='gammall'))*T).^2;
 wh     = exp(x((find(list.test=='wh')-1)*T+1:(find(list.test=='wh'))*T));
 wl     = exp(x((find(list.test=='wl')-1)*T+1:(find(list.test=='wl'))*T));
 ws     = exp(x((find(list.test=='ws')-1)*T+1:(find(list.test=='ws'))*T));
 pg     = exp(x((find(list.test=='pg')-1)*T+1:(find(list.test=='pg'))*T));
 pn     = exp(x((find(list.test=='pn')-1)*T+1:(find(list.test=='pn'))*T));
 pee     = exp(x((find(list.test=='pee')-1)*T+1:(find(list.test=='pee'))*T));
 pf     = exp(x((find(list.test=='pf')-1)*T+1:(find(list.test=='pf'))*T));
 lambdaa= exp(x((find(list.test=='lambdaa')-1)*T+1:(find(list.test=='lambdaa'))*T));

% auxiliary variables 

muu = C.^(-thetaa);
E  = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1));
xn=pn.*alphan.*N;
xg=pg.*alphag.*G;
xf=pf.*(1-pol(list.pol=='tauf')).*alphaf.*F;
Cincome=Y-xn-xf-xg;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +tauf.*pf.*F;
wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

Emnet     = omegaa.*F-deltaa; % net emissions

% utility
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
 SWF = Utilcon-Utillab;

% test market clearing
if max(abs(C-Cincome))>1e-10
    error('market clearing does not hold')
else
    fprintf('goods market cleared!')
end

% test variables read in properly
% xx=eval(symms.choice);
% xx=xx(:);
% guess_trans=trans_guess(indexx('LF'), xx, params, list.params);
% f=laissez_faire_nows(guess_trans, params, list, pol, laggs);
% 
% if (max(abs(f)))>1e-8
%     fprintf('f only solved at less than 1e-8')
% else
%     fprintf('saved variables are correct!')
% end

% save stuff
LF_SIM= eval(symms.allvars);
end
 