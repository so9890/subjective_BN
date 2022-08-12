function LF_SIM=aux_solutionLF_VECT_xgrowth(x, list, symms,varrs, params, T , indic, laggs)
read_in_params;

tauf=varrs(list.allvars=='tauf', :)';
taus=varrs(list.allvars=='taus', :)';
taul=varrs(list.allvars=='taul', :)';

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

%- initial condition

Af_lag=laggs(list.init=='Af0');
An_lag=laggs(list.init=='An0');
Ag_lag=laggs(list.init=='Ag0');

Af=zeros(size(C));
Ag=zeros(size(F));
An=zeros(size(G));

for i=1:T
    An(i)=(1+vn)*An_lag;
    Ag(i)=(1+vg)*Ag_lag;
    Af(i)=(1+vf)*Af_lag;
    %- update laggs
    An_lag=An(i);
    Af_lag=Af(i);
    Ag_lag=Ag(i);
end
sff=zeros(size(gammalh));
sg=zeros(size(gammalh));
sn=zeros(size(gammalh));
S=zeros(size(gammalh));
ws=zeros(size(gammalh));
wsg=ws; wsf=ws; wsn=ws;
gammasg=ws; gammasn=ws; gammasf=ws;
%  Af     = exp(x((find(list.test=='Af')-1)*T+1:(find(list.test=='Af'))*T));
%  Ag     = exp(x((find(list.test=='Ag')-1)*T+1:(find(list.test=='Ag'))*T));
%  An     = exp(x((find(list.test=='An')-1)*T+1:(find(list.test=='An'))*T));
%  sff     = (x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T)).^2;
%  sg     = (x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)).^2;
%  sn     =(x((find(list.test=='sn')-1)*T+1:(find(list.test=='sn'))*T)).^2;
%  S     = (x((find(list.test=='S')-1)*T+1:(find(list.test=='S'))*T)).^2;
%  gammas     = (x((find(list.test=='gammas')-1)*T+1:(find(list.test=='gammas'))*T)).^2;
% 
%  
%  ws     = (x((find(list.test=='ws')-1)*T+1:(find(list.test=='ws'))*T)).^2;
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
xf=pf.*(1-tauf).*alphaf.*F;

A   = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Ln      = hhn.^thetan.*hln.^(1-thetan);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf); 
    
    if indic.notaul<2 % tauf redistributed via income tax
        SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
    else
        SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul));

    end
else
    if indic.notaul<2
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
    else
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul));

    end
end
% gov con
        if 2<=indic.notaul && indic.notaul <4
            GovCon = tauf.*pf.*F;
        else
            GovCon =zeros(size(F));
        end
% lump sum transfers
if indic.notaul >=4
    Tls =tauf.*pf.*F;
else
    Tls =zeros(size(F));
end


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
if indic.noskill==0
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
 Utilsci = zeros(size(gammalh));

 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
Cincome=Y-xn-xf-xg-GovCon;
if max(abs(C-Cincome))>1e-8
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
 