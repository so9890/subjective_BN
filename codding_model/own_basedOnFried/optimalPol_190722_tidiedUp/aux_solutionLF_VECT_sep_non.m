function LF_SIM=aux_solutionLF_VECT_sep_non(x, list, symms,varrs, params, T , indic)

read_in_params;

tauf=varrs(list.allvars=='tauf', :)';
taus=varrs(list.allvars=='taus', :)';
taul=varrs(list.allvars=='taul', :)';

gammalh = x((find(list.test=='gammalh')-1)*T+1:(find(list.test=='gammalh'))*T).^2;

if indic.noskill==0
 hhf    = exp(x((find(list.test=='hhf')-1)*T+1:(find(list.test=='hhf'))*T));
 hhg    = exp(x((find(list.test=='hhg')-1)*T+1:(find(list.test=='hhg'))*T));
 hln    = exp(x((find(list.test=='hln')-1)*T+1:(find(list.test=='hln'))*T));
 hlf    = exp(x((find(list.test=='hlf')-1)*T+1:(find(list.test=='hlf'))*T));
 hlg    = exp(x((find(list.test=='hlg')-1)*T+1:(find(list.test=='hlg'))*T));
 hl     = upbarH./(1+exp(x((find(list.test=='HL')-1)*T+1:find(list.test=='HL')*T)));
 hh     = upbarH./(1+exp(x((find(list.test=='HH')-1)*T+1:(find(list.test=='HH'))*T)));
 gammall = x((find(list.test=='gammall')-1)*T+1:(find(list.test=='gammall'))*T).^2; 
 wh     = exp(x((find(list.test=='wh')-1)*T+1:(find(list.test=='wh'))*T));
 wl     = exp(x((find(list.test=='wl')-1)*T+1:(find(list.test=='wl'))*T));

else
    
 Lg    = exp(x((find(list.test=='Lg')-1)*T+1:(find(list.test=='Lg'))*T));
 Lf    = exp(x((find(list.test=='Lf')-1)*T+1:(find(list.test=='Lf'))*T));
 w     = exp(x((find(list.test=='w')-1)*T+1:(find(list.test=='w'))*T));
 h     = upbarH./(1+exp(x((find(list.test=='H')-1)*T+1:find(list.test=='H')*T)));
 wh=w; wl=w; gammall=gammalh; hh=h; hl=h; hhn=zeros(size(h));hhg=zeros(size(h));
 hhf=zeros(size(h)); hln=zeros(size(h)); hlg=zeros(size(h)); hlf=zeros(size(h)); 

end

if indic.ineq==0
    if indic.BN==0
     C      = exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T));
    else
     C     = B./(1+exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T)));
    end
Ch=C; Cl=C; 
else
    if indic.BN==0
     Ch      = exp(x((find(list.test=='Ch')-1)*T+1:(find(list.test=='Ch'))*T));
     Cl      = exp(x((find(list.test=='Cl')-1)*T+1:(find(list.test=='Cl'))*T));
    else
     Ch     = Bh./(1+exp(x((find(list.test=='Ch')-1)*T+1:(find(list.test=='Ch'))*T)));
     Cl     = Bl./(1+exp(x((find(list.test=='Cl')-1)*T+1:(find(list.test=='Cl'))*T)));
    end
C=Ch;
end
 F      = exp(x((find(list.test=='F')-1)*T+1:(find(list.test=='F'))*T));
 G      = exp(x((find(list.test=='G')-1)*T+1:(find(list.test=='G'))*T));
 Af     = exp(x((find(list.test=='Af')-1)*T+1:(find(list.test=='Af'))*T));
 Ag     = exp(x((find(list.test=='Ag')-1)*T+1:(find(list.test=='Ag'))*T));
 sff    = (x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T)).^2;
 sg     = (x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)).^2;
 gammasg     = (x((find(list.test=='gammasg')-1)*T+1:(find(list.test=='gammasg'))*T)).^2;
 gammasf     = (x((find(list.test=='gammasf')-1)*T+1:(find(list.test=='gammasf'))*T)).^2;
 wsg     = (x((find(list.test=='wsg')-1)*T+1:(find(list.test=='wsg'))*T)).^2;
 wsf     = (x((find(list.test=='wsf')-1)*T+1:(find(list.test=='wsf'))*T)).^2;

 pg     = exp(x((find(list.test=='pg')-1)*T+1:(find(list.test=='pg'))*T));
 pf     = exp(x((find(list.test=='pf')-1)*T+1:(find(list.test=='pf'))*T));
 lambdaa= exp(x((find(list.test=='lambdaa')-1)*T+1:(find(list.test=='lambdaa'))*T));

% auxiliary variables 
if indic.ineq==0
    if indic.BN==0
        muu      = C.^(-thetaa); % same equation in case thetaa == 1
    else
        muu =-(C-B).^(zetaa-1);
    end
muhh=muu;
muul=muu;
else
    if indic.BN==0
        muuh      = Ch.^(-thetaa); % same equation in case thetaa == 1
        muul      = Cl.^(-thetaa); % same equation in case thetaa == 1
    else
        muul =-(Cl-Bl).^(zetaa-1);
        muuh =-(Ch-Bh).^(zetaa-1);
    end
muu=muuh;
end

E  = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1)); 
Y = E; %(deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1));

%- non neutral placeholders
N       = zeros(size(E));
pn      = zeros(size(E));
xn      = zeros(size(E));
sn      = zeros(size(E));
gammasn = zeros(size(E));
Ln      = zeros(size(E));
wln     = zeros(size(E));
hhn     = zeros(size(E));
hln     = zeros(size(E));
wsn     = zeros(size(E));
An      = zeros(size(E));
pee     = zeros(size(E));
%-results continue
xg=pg.*alphag.*G;
xf=pf.*(1-tauf).*alphaf.*F;
Cincome=Y-xf-xg;

A   = (rhof*Af+rhog*Ag)/(rhof+rhog);
S=sff+sg;

if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
else
    SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
end

wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

Emnet     = omegaa.*F-deltaa; % net emissions

% utility
if indic.ineq==0
    if indic.BN==0
        if thetaa~=1
            Utilcon = (C.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = log(C);
        end
    else
        Utilcon=-(C-B).^(zetaa)./(zetaa); 
    end
else
    if indic.BN==0
        if thetaa~=1
            Utilcon = zh.*(Ch.^(1-thetaa))./(1-thetaa)+(1-zh).*(Cl.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = zh.*log(Ch)+(1-zh).*log(Cl);
        end
    else
        Utilcon=zh.*(-(Ch-Bh).^(zetaa)./(zetaa))+(1-zh).*(-(Cl-Bl).^(zetaa)./(zetaa)); 
    end

end


if indic.noskill==0
     Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
else
     Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
end
 Utilsci = chiis*sn.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sff.^(1+sigmaas)./(1+sigmaas);

 SWF = Utilcon-Utillab-Utilsci;

% test market clearing
if indic.ineq==0
    diff=C-Cincome;
else
    diff=zh.*Ch+(1-zh).*Cl-Cincome;
end

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
 