function LF_SIM=aux_solutionLF_VECT_sep(x, list, symms,varrs, params, T , indic)

read_in_params;

if indic.limit_LF==0
    tauf=varrs(list.allvars=='tauf', :)'; 
else
    tauf=(x((find(list.test=='tauf')-1)*T+1:(find(list.test=='tauf'))*T));
end
taus=varrs(list.allvars=='taus', :)';

gammalh = x((find(list.test=='gammalh')-1)*T+1:(find(list.test=='gammalh'))*T).^2;

if indic.noskill==0
 hhf    = exp(x((find(list.test=='hhf')-1)*T+1:(find(list.test=='hhf'))*T));
 hhg    = exp(x((find(list.test=='hhg')-1)*T+1:(find(list.test=='hhg'))*T));
 hhn    = exp(x((find(list.test=='hhn')-1)*T+1:(find(list.test=='hhn'))*T));
 hln    = exp(x((find(list.test=='hln')-1)*T+1:(find(list.test=='hln'))*T));
 hlf    = exp(x((find(list.test=='hlf')-1)*T+1:(find(list.test=='hlf'))*T));
 hlg    = exp(x((find(list.test=='hlg')-1)*T+1:(find(list.test=='hlg'))*T));
 hl     = upbarH./(1+exp(x((find(list.test=='HL')-1)*T+1:find(list.test=='HL')*T)));
 hh     = upbarH./(1+exp(x((find(list.test=='HH')-1)*T+1:(find(list.test=='HH'))*T)));
 gammall = x((find(list.test=='gammall')-1)*T+1:(find(list.test=='gammall'))*T).^2; 
 wh     = exp(x((find(list.test=='wh')-1)*T+1:(find(list.test=='wh'))*T));
 wl     = exp(x((find(list.test=='wl')-1)*T+1:(find(list.test=='wl'))*T));

else
    
 Ln    = exp(x((find(list.test=='Ln')-1)*T+1:(find(list.test=='Ln'))*T));
 Lg    = exp(x((find(list.test=='Lg')-1)*T+1:(find(list.test=='Lg'))*T));
 Lf    = exp(x((find(list.test=='Lf')-1)*T+1:(find(list.test=='Lf'))*T));
 w     = exp(x((find(list.test=='w')-1)*T+1:(find(list.test=='w'))*T));
 h     = upbarH./(1+exp(x((find(list.test=='H')-1)*T+1:find(list.test=='H')*T)));
 wh=w; wl=w; gammall=gammalh; hh=h; hl=h; hhn=zeros(size(h));hhg=zeros(size(h));
 hhf=zeros(size(h)); hln=zeros(size(h)); hlg=zeros(size(h)); hlf=zeros(size(h)); 

end

C      = exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T));

 F      = exp(x((find(list.test=='F')-1)*T+1:(find(list.test=='F'))*T));
 G      = exp(x((find(list.test=='G')-1)*T+1:(find(list.test=='G'))*T));
 Af     = exp(x((find(list.test=='Af')-1)*T+1:(find(list.test=='Af'))*T));
 Ag     = exp(x((find(list.test=='Ag')-1)*T+1:(find(list.test=='Ag'))*T));
 An     = exp(x((find(list.test=='An')-1)*T+1:(find(list.test=='An'))*T));
  if indic.xgrowth==0 && indic.noskill==1
    sff    = upbarH./(1+exp(x((find(list.test=='sff')-1)*T+1:find(list.test=='sff')*T)));
    sg     = upbarH./(1+exp(x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)));
    sn     = upbarH./(1+exp(x((find(list.test=='sn')-1)*T+1:(find(list.test=='sn'))*T)));
    se     = upbarH./(1+exp(x((find(list.test=='se')-1)*T+1:(find(list.test=='se'))*T)));

 else
    sff     = (x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T)).^2;
    sg     = (x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)).^2;
    sn     = (x((find(list.test=='sn')-1)*T+1:(find(list.test=='sn'))*T)).^2;
    se     = (x((find(list.test=='se')-1)*T+1:(find(list.test=='se'))*T)).^2;

  end
  if indic.sep==0
     S     = upbarS./(1+exp(x((find(list.test=='S')-1)*T+1:(find(list.test=='S'))*T)));
     gammas = (x((find(list.test=='gammas')-1)*T+1:(find(list.test=='gammas'))*T)).^2;
     ws     = (x((find(list.test=='ws')-1)*T+1:(find(list.test=='ws'))*T)).^2;
  elseif indic.sep==2
    ws =wspar*ones(size(C));
    S=sff+sg+sn;
    gammas=zeros(size(S));
  else
      S=(sg+sff+sn);
  end
 gammasg     = (x((find(list.test=='gammasg')-1)*T+1:(find(list.test=='gammasg'))*T)).^2;
 gammasn     = (x((find(list.test=='gammasn')-1)*T+1:(find(list.test=='gammasn'))*T)).^2;
 gammasf     = (x((find(list.test=='gammasf')-1)*T+1:(find(list.test=='gammasf'))*T)).^2;
 wsn     = (x((find(list.test=='wsn')-1)*T+1:(find(list.test=='wsn'))*T)).^2;
 wsg     = (x((find(list.test=='wsg')-1)*T+1:(find(list.test=='wsg'))*T)).^2;
 wsf     = (x((find(list.test=='wsf')-1)*T+1:(find(list.test=='wsf'))*T)).^2;

 pg     = exp(x((find(list.test=='pg')-1)*T+1:(find(list.test=='pg'))*T));
 pn     = exp(x((find(list.test=='pn')-1)*T+1:(find(list.test=='pn'))*T));
 pee     = exp(x((find(list.test=='pee')-1)*T+1:(find(list.test=='pee'))*T));
 pf     = exp(x((find(list.test=='pf')-1)*T+1:(find(list.test=='pf'))*T));
 if indic.notaul ~=6
    lambdaa= (x((find(list.test=='lambdaa')-1)*T+1:(find(list.test=='lambdaa'))*T));
    taul=varrs(list.allvars=='taul', :)';

 else
    taul= (x((find(list.test=='lambdaa')-1)*T+1:(find(list.test=='lambdaa'))*T));   
    lambdaa=varrs(list.allvars=='lambdaa', :)';

 end
% auxiliary variables 
muu      = C.^(-thetaa); % same equation in case thetaa == 1
  
if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf); 
else
    hh=h; hl=h; wh=w; wl=w; % this should suffice to have governmenta budget correct
end


if indic.Sun~=2
    helpS=zeros(size(F)); % helper for gov revenues from taxing scientists
elseif indic.Sun==2
    helpS=(ws.*S-lambdaa.*(ws.*S).^(1-taul)); % version with lump sum transfers to finance subsidies on machine producers
end
if indic.notaul<2 || ...
   indic.notaul == 6 % tauf redistributed via income tax

    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +helpS...
        +tauf.*F;
    Tls =zeros(size(F));    
    GovCon =zeros(size(F));
elseif indic.notaul == 2 ||...
        indic.notaul==3 %2,3,4,5,7
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
           +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
           +helpS;
    GovCon = tauf.*F; % GovCon = env tax consumed by government
    Tls =zeros(size(F)); 
elseif indic.notaul == 4 || indic.notaul ==5
    if indic.noskill==0
        SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +helpS;
    else
        SGov = w.*h-lambdaa.*(w.*h).^(1-taul) +helpS;
    end
    GovCon =zeros(size(F));
    Tls  = tauf.*F;
elseif indic.notaul >= 7 % earmarking
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
         +((1-zh))*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
         +helpS;
    GovCon =zeros(size(F));
    Tls = zeros(size(F));
end
if indic.xgrowth==0
    Tlsall = Tls-ws.*S;
else
    Tlsall=Tls;
end


E  = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
N  =  (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y = (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1));
xn=pn.*alphan.*N;
xg=pg.*(1+taus).*alphag.*G;
xf=pf.*alphaf.*F;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
if indic.sep~=0
    S=sff+sg+sn;
end
wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = (pg.*(1+taus)).^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 

Emnet     = omegaa.*F-deltaa; % net emissions

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end
   
if indic.noskill==0
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+((1-zh)).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
if indic.sep~=3 && indic.sep~=0
    Utilsci = chiis*sn.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sff.^(1+sigmaas)./(1+sigmaas);
elseif indic.sep==3
     Utilsci = chiis*sn.^(1+sigmaas)./(1+sigmaas)+chiis*se.^(1+sigmaas)./(1+sigmaas);
elseif indic.sep==0
    Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
end
 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
Cincome=Y-xn-xf-xg-GovCon-SGov;

diff=C-Cincome;

if max(abs(diff))>1e-7
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
 