function f=laissez_faireVECT_sep_non(x, params, list, varrs, laggs,T, indic)

% called by script 'test_results.m'
% takes optimal policy results as input

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;

tauf=varrs(list.allvars=='tauf', :)';
taus=varrs(list.allvars=='taus', :)';
taul=varrs(list.allvars=='taul', :)';

% choice variables
%- transform variables directly instead of in code
if indic.noskill==0
 hhf    = exp(x((find(list.test=='hhf')-1)*T+1:(find(list.test=='hhf'))*T));
 hhg    = exp(x((find(list.test=='hhg')-1)*T+1:(find(list.test=='hhg'))*T));
 hlf    = exp(x((find(list.test=='hlf')-1)*T+1:(find(list.test=='hlf'))*T));
 hlg    = exp(x((find(list.test=='hlg')-1)*T+1:(find(list.test=='hlg'))*T));
 gammall = x((find(list.test=='gammall')-1)*T+1:(find(list.test=='gammall'))*T).^2;

 wh     = exp(x((find(list.test=='wh')-1)*T+1:(find(list.test=='wh'))*T));
 wl     = exp(x((find(list.test=='wl')-1)*T+1:(find(list.test=='wl'))*T));
 hl     = upbarH./(1+exp(x((find(list.test=='HL')-1)*T+1:find(list.test=='HL')*T)));
 hh     = upbarH./(1+exp(x((find(list.test=='HH')-1)*T+1:(find(list.test=='HH'))*T)));
else
 w     = exp(x((find(list.test=='w')-1)*T+1:(find(list.test=='w'))*T));
 h     = upbarH./(1+exp(x((find(list.test=='H')-1)*T+1:find(list.test=='H')*T)));
 Lf    = exp(x((find(list.test=='Lf')-1)*T+1:(find(list.test=='Lf'))*T));
 Lg    = exp(x((find(list.test=='Lg')-1)*T+1:(find(list.test=='Lg'))*T));
end

if indic.ineq==0
    if indic.BN==0
     C      = exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T));
    else
     C     = B./(1+exp(x((find(list.test=='C')-1)*T+1:(find(list.test=='C'))*T)));
    end
else
    if indic.BN==0
     Ch      = exp(x((find(list.test=='Ch')-1)*T+1:(find(list.test=='Ch'))*T));
     Cl      = exp(x((find(list.test=='Cl')-1)*T+1:(find(list.test=='Cl'))*T));
    else
     Cl     = Bl./(1+exp(x((find(list.test=='Cl')-1)*T+1:(find(list.test=='Cl'))*T)));
     Ch     = Bh./(1+exp(x((find(list.test=='Ch')-1)*T+1:(find(list.test=='Ch'))*T)));
    end
end
 F      = exp(x((find(list.test=='F')-1)*T+1:(find(list.test=='F'))*T));
 G      = exp(x((find(list.test=='G')-1)*T+1:(find(list.test=='G'))*T));
 Af     = exp(x((find(list.test=='Af')-1)*T+1:(find(list.test=='Af'))*T));
 Ag     = exp(x((find(list.test=='Ag')-1)*T+1:(find(list.test=='Ag'))*T));
 sff     = (x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T)).^2;
 sg     = (x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)).^2;
 gammasf = x((find(list.test=='gammasf')-1)*T+1:(find(list.test=='gammasf'))*T).^2;
 gammasg = x((find(list.test=='gammasg')-1)*T+1:(find(list.test=='gammasg'))*T).^2;

 
 gammalh = x((find(list.test=='gammalh')-1)*T+1:(find(list.test=='gammalh'))*T).^2;
 wsg     = (x((find(list.test=='wsg')-1)*T+1:(find(list.test=='wsg'))*T)).^2;
 wsf     = (x((find(list.test=='wsf')-1)*T+1:(find(list.test=='wsf'))*T)).^2;
 pg     = exp(x((find(list.test=='pg')-1)*T+1:(find(list.test=='pg'))*T));
 pf     = exp(x((find(list.test=='pf')-1)*T+1:(find(list.test=='pf'))*T));
 lambdaa= exp(x((find(list.test=='lambdaa')-1)*T+1:(find(list.test=='lambdaa'))*T));

%% - read in auxiliary equations
%- initial condition
Af_lag=[laggs(list.init=='Af0'); Af(1:T-1)];
An_lag=zeros(size(Af_lag));
Ag_lag=[laggs(list.init=='Ag0'); Ag(1:T-1)];
A_lag   = (rhof*Af_lag+rhog*Ag_lag)/(rhof+rhog);


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

if indic.ineq==0
    if indic.BN==0
        muu      = C.^(-thetaa); % same equation in case thetaa == 1
    else
        muu =-(C-B).^(zetaa-1);
    end
else
    if indic.BN==0
        muuh      = Ch.^(-thetaa); % same equation in case thetaa == 1
        muul      = Cl.^(-thetaa); % same equation in case thetaa == 1

    else
        muuh =-(Ch-Bh).^(zetaa-1);
        muul =-(Cl-Bl).^(zetaa-1);
    end
end
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
        
            % subsidies and profits and wages scientists cancel
Y       =  E; %(deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
% wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
if indic.noskill==0
    if indic.ineq==0
        q=q+1;
        f(q:T)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f((q-1)*T+1:T*q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl
        %3- budget
        q=q+1;
        f((q-1)*T+1:T*q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov-C; %=> determines C
    else  
        q=q+1;
        f((q-1)*T+1:T*q)= chii*hh.^(sigmaa+taul)- ((muuh.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f((q-1)*T+1:T*q)= chii*hl.^(sigmaa+taul) - ((muul.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall.*hl.^taul); %=> determines hl
    
        %3- budget
        q=q+1;
        f((q-1)*T+1:T*q) = lambdaa.*(wl.*hl).^(1-taul)+SGov-Cl; %=> determines C
        q=q+1;
        f((q-1)*T+1:T*q) = lambdaa.*(wh.*hh).^(1-taul)+SGov-Ch;
    end
else
    q=q+1;
    f(q:T)= chii*h.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(w).^(1-taul))-gammalh.*h.^taul); %=> determines hh
   %3- budget
    q=q+1;
    f((q-1)*T+1:T*q) = lambdaa.*(w.*h).^(1-taul)+SGov-C; %=> determines C

end
%4- output fossil
q=q+1;
f((q-1)*T+1:T*q) = ((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%6- output green
q=q+1;
f((q-1)*T+1:T*q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%7- demand green scientists
q=q+1;
f((q-1)*T+1:T*q)= wsf - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
%8
q=q+1;
f((q-1)*T+1:T*q)= wsg - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus).*Ag);

%11
q=q+1;
f((q-1)*T+1:T*q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
%12
q=q+1;
f((q-1)*T+1:T*q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
 
%13- optimality labour input producers
if indic.noskill==0
    %14
    q=q+1;
    f((q-1)*T+1:T*q)= thetag*Lg.*wlg-wh.*hhg;
    
    %16
    q=q+1;
    f((q-1)*T+1:T*q)=(1-thetag)*Lg.*wlg-wl.*hlg;
    %17- demand skill
    q=q+1;
    f((q-1)*T+1:T*q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
            ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
    %18
    q=q+1;
    f((q-1)*T+1:T*q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
            ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
else
    q=q+1; % labour demand fossil
    f((q-1)*T+1:T*q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    q=q+1; % labour demand fossil
    f((q-1)*T+1:T*q) =  w  - wlg; % labour demand fossil

end
% prices and wages
%19- optimality energy producers
q=q+1;
f((q-1)*T+1:T*q) = pf.*F.^(1/eppse)- (G).^(1/eppse).*pg; 



%- definitions prices
%20
q=q+1;
f((q-1)*T+1:T*q) = 1 - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition

%22- market clearing (consumption good=> numeraire)
if indic.noskill==0
    q=q+1;
    f((q-1)*T+1:T*q) = zh.*hh-( hhf+hhg); % high skill market clearing
    %23
    q=q+1;
    f((q-1)*T+1:T*q) = (1-zh).*hl-( hlf+hlg); % low skill market clearing
else
    q=q+1;
    f((q-1)*T+1:T*q) = Lf+Lg-h;  
end

%24
% q=q+1;
% f((q-1)*T+1:T*q) = S-(sn+sff+sg);
% scientists supply
q=q+1;
f((q-1)*T+1:T*q)= (chiis).*sff.^sigmaas-(wsf-gammasf); % scientist hours supply
q=q+1;
f((q-1)*T+1:T*q)= (chiis).*sg.^sigmaas-((wsg-gammasg));

q=q+1;
f((q-1)*T+1:T*q)= gammasf.*(upbarH-sff);
q=q+1;
f((q-1)*T+1:T*q)= gammasg.*(upbarH-sg);

%13- Kuhn Tucker Labour supply
if indic.noskill==0
    %25
    q=q+1;
    f((q-1)*T+1:T*q)= gammalh.*(upbarH-hh);
    %26
    q=q+1;
    f((q-1)*T+1:T*q)= gammall.*(upbarH-hl);
else
    q=q+1;
    f((q-1)*T+1:T*q)= gammalh.*(upbarH-h);
end
q=q+1;
f((q-1)*T+1:T*q)= SGov;

% fprintf('number equations: %d; number variables %d', q, length(list.test));
end