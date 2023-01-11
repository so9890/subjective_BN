function f=laissez_faireVECT_sep_NoRed(x, params, list, varrs, laggs,T, indic, Ems)

% Policy version with GOV=tauf pf F
% can handle version with and without taul

% called by script 'test_results.m'
% takes optimal policy results as input

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;

if indic.limit_LF==0
    tauf=varrs(list.allvars=='tauf', :)'; 
else
    tauf=(x((find(list.test=='tauf')-1)*T+1:(find(list.test=='tauf'))*T));
end

% taus is exogenous
taus=varrs(list.allvars=='taus', :)';

% choice variables
%- transform variables directly instead of in code
if indic.noskill==0
 hhf    = exp(x((find(list.test=='hhf')-1)*T+1:(find(list.test=='hhf'))*T));
 hhg    = exp(x((find(list.test=='hhg')-1)*T+1:(find(list.test=='hhg'))*T));
 hhn    = exp(x((find(list.test=='hhn')-1)*T+1:(find(list.test=='hhn'))*T));
 hln    = exp(x((find(list.test=='hln')-1)*T+1:(find(list.test=='hln'))*T));
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
 Ln    = exp(x((find(list.test=='Ln')-1)*T+1:(find(list.test=='Ln'))*T));
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
        sff    = (x((find(list.test=='sff')-1)*T+1:(find(list.test=='sff'))*T)).^2;
        sg     = (x((find(list.test=='sg')-1)*T+1:(find(list.test=='sg'))*T)).^2;
        sn     = (x((find(list.test=='sn')-1)*T+1:(find(list.test=='sn'))*T)).^2;
        se     = (x((find(list.test=='se')-1)*T+1:(find(list.test=='se'))*T)).^2;
     end
S    = upbarS./(1+exp(x((find(list.test=='S')-1)*T+1:(find(list.test=='S'))*T)));
 if indic.sep~=2
     if indic.sep~=0
         gammasf = x((find(list.test=='gammasf')-1)*T+1:(find(list.test=='gammasf'))*T).^2;
         gammasg = x((find(list.test=='gammasg')-1)*T+1:(find(list.test=='gammasg'))*T).^2;
         gammasn = x((find(list.test=='gammasn')-1)*T+1:(find(list.test=='gammasn'))*T).^2;
         wsn     = (x((find(list.test=='wsn')-1)*T+1:(find(list.test=='wsn'))*T)).^2;
         wsg     = (x((find(list.test=='wsg')-1)*T+1:(find(list.test=='wsg'))*T)).^2;
         wsf     = (x((find(list.test=='wsf')-1)*T+1:(find(list.test=='wsf'))*T)).^2;
     else
        gammas = x((find(list.test=='gammas')-1)*T+1:(find(list.test=='gammas'))*T).^2;
        ws     = (x((find(list.test=='ws')-1)*T+1:(find(list.test=='ws'))*T)).^2;
        
     end
 else
     ws =wspar;
%      wsg=wsgpar;
%      wsf=wsfpar;
 end
 
 gammalh = x((find(list.test=='gammalh')-1)*T+1:(find(list.test=='gammalh'))*T).^2;
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
%% - read in auxiliary equations
%- initial condition
Af_lag=[laggs(list.init=='Af0'); Af(1:T-1)];
An_lag=[laggs(list.init=='An0'); An(1:T-1)];
Ag_lag=[laggs(list.init=='Ag0'); Ag(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);

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
           +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
           +helpS;
    GovCon = tauf.*F; % GovCon = env tax consumed by government
    Tls =zeros(size(F)); 
elseif indic.notaul == 4 || indic.notaul ==5
    if indic.noskill==0
        SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
        +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
        +helpS;
    else
        SGov = w.*h-lambdaa.*(w.*h).^(1-taul) +helpS;
    end
    GovCon =zeros(size(F));
    Tls  = tauf.*F;
elseif indic.notaul >= 7 % earmarking
    SGov = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
         +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
         +helpS;
    GovCon =zeros(size(F));
    Tls = zeros(size(F));
end
if indic.xgrowth==0
    Tlsall=Tls-(ws.*sff+ws.*sn+ws.*(1-taus).*sg);
else
    Tlsall=Tls;
end

muu      = C.^(-thetaa); % same equation in case thetaa == 1
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       =  (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       =  (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = (pg).^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
% wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg*Ag;

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
if indic.noskill==0

        q=q+1;
        f(q:T)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f((q-1)*T+1:T*q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl
        %3- budget
        q=q+1;
        if indic.Sun~=2 || indic.xgrowth==1
            if indic.xgrowth==0
               f((q-1)*T+1:T*q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+ws.*S+Tlsall-C;
            else
               f((q-1)*T+1:T*q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+Tlsall-C;
            end
        elseif indic.Sun==2
            f((q-1)*T+1:T*q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+lambdaa.*(ws.*S).^(1-taul)+Tlsall-C; %=> determines C
        end
else
        q=q+1;
        f(q:T)= chii*h.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(w).^(1-taul))-gammalh.*h.^taul); %=> determines hh
       %3- budget
        q=q+1;
        if indic.Sun~=2 || indic.xgrowth==1
            if indic.xgrowth==0
                f((q-1)*T+1:T*q) = lambdaa.*(w.*h).^(1-taul)+ws.*S+Tlsall-C; %=> determines C
            else
               f((q-1)*T+1:T*q) = lambdaa.*(w.*h).^(1-taul)+Tlsall-C; %=> determines C
            end
       elseif indic.Sun==2
           f((q-1)*T+1:T*q) = lambdaa.*(w.*h).^(1-taul)+lambdaa.*(ws.*S).^(1-taul)+Tlsall-C; %=> determines C
        end
end
%4- output fossil
q=q+1;
f((q-1)*T+1:T*q) = (alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%5- output neutral
q=q+1;
f((q-1)*T+1:T*q) = N-An.*Ln.*(pn.*alphan).^(alphan./(1-alphan)); 

%6- output green
q=q+1;
f((q-1)*T+1:T*q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%7- demand green scientists

if indic.xgrowth==0
    if indic.sep~=0 && indic.sep~=2
            q=q+1;
            f((q-1)*T+1:T*q)= wsf - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
            %8
            q=q+1;
            f((q-1)*T+1:T*q)= wsg - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus).*Ag);
            %9
            q=q+1;
            f((q-1)*T+1:T*q)= wsn - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);
        if indic.sep==3

            q=q+1;
            f((q-1)*T+1:T*q)= wsg -wsf;

        end
    elseif indic.sep==0 || indic.sep==2
            q=q+1;
            f((q-1)*T+1:T*q)= ws - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 
            %8
            q=q+1;
            f((q-1)*T+1:T*q)= ws - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus).*Ag);
            %9
            q=q+1;
            f((q-1)*T+1:T*q)= ws - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag)./(rhon^etaa.*An);
    end
end
%10- LOM technology
q=q+1;
f((q-1)*T+1:T*q) = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
%11
q=q+1;
f((q-1)*T+1:T*q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
%12
q=q+1;
f((q-1)*T+1:T*q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
 
%13- optimality labour input producers
if indic.noskill==0
    q=q+1;
    f((q-1)*T+1:T*q) = thetan*Ln.*wln-wh.*hhn;
    %14
    q=q+1;
    f((q-1)*T+1:T*q)= thetag*Lg.*wlg-wh.*hhg;
    %15
    q=q+1;
    f((q-1)*T+1:T*q)=(1-thetan)*Ln.*wln-wl.*hln;
    %16
    q=q+1;
    f((q-1)*T+1:T*q)=(1-thetag)*Lg.*wlg-wl.*hlg;
    %17- demand skill
    q=q+1;
    f((q-1)*T+1:T*q) = wh - thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
            (pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
    %18
    q=q+1;
    f((q-1)*T+1:T*q) = wl-(1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf)*alphaf^(alphaf/(1-alphaf)).*...
            (pf).^(1/(1-alphaf)).*Af;
else
    q=q+1; % labour demand fossil
    f((q-1)*T+1:T*q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    q=q+1;
    f((q-1)*T+1:T*q) = Ln -pn.*(1-alphan).*N./w; % labour demand neutral 
    q=q+1;
    f((q-1)*T+1:T*q) = Lg - pg.*(1-alphag).*G./w;
end
% prices and wages
%19- optimality energy producers
q=q+1;
f((q-1)*T+1:T*q) = (pf+tauf).*F.^(1/eppse)- (G).^(1/eppse).*pg; 


%- definitions prices
%20
q=q+1;
f((q-1)*T+1:T*q) = pee - ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
f((q-1)*T+1:T*q) =  1-(deltay.*pee.^(1-eppsy)+(1-deltay).*pn.^(1-eppsy)).^(1/(1-eppsy));

%22- market clearing (consumption good=> numeraire)
if indic.noskill==0
    q=q+1;
    f((q-1)*T+1:T*q) = zh.*hh-(hhn + hhf+hhg); % high skill market clearing
    %23
    q=q+1;
    f((q-1)*T+1:T*q) = (1-zh).*hl-(hln + hlf+hlg); % low skill market clearing
else
    q=q+1;
    f((q-1)*T+1:T*q) = Lf+Lg+Ln-h;  
end

%24
% q=q+1;
% f((q-1)*T+1:T*q) = S-(sn+sff+sg);
% scientists supply
if indic.xgrowth==0 && indic.sep~=2
    if indic.sep==1
        q=q+1;
        f((q-1)*T+1:T*q)= (chiis).*sff.^sigmaas-(muu.*wsf-gammasf); % scientist hours supply
        q=q+1;
        f((q-1)*T+1:T*q)= (chiis).*sg.^sigmaas-((muu.*wsg-gammasg));
        
        q=q+1;
        f((q-1)*T+1:T*q)= gammasf.*(sff-upbarH);
        q=q+1;
        f((q-1)*T+1:T*q)= gammasg.*(sg-upbarH);
        
        q=q+1;
        f((q-1)*T+1:T*q)= (chiis).*sn.^sigmaas-((muu.*wsn-gammasn));
        q=q+1;
        f((q-1)*T+1:T*q)= gammasn.*(sn-upbarH);
    elseif indic.sep==3
        q=q+1;
        f((q-1)*T+1:T*q)= (chiis).*se.^sigmaas-((muu.*wsg-gammasg));
        q=q+1;
        f((q-1)*T+1:T*q)= gammasg.*(se-upbarH);  
        q=q+1;
        f((q-1)*T+1:T*q)= sff+sg-se;
        
        q=q+1;
        f((q-1)*T+1:T*q)= (chiis).*sn.^sigmaas-((muu.*wsn-gammasn));
        q=q+1;
        f((q-1)*T+1:T*q)= gammasn.*(sn-upbarH);
    elseif indic.sep==0
        q=q+1;
        f((q-1)*T+1:T*q)= sff+sg+sn-S; % determines wage in neutral sector
        q=q+1;
        if indic.Sun==0
            f((q-1)*T+1:T*q)= (chiis).*S.^sigmaas-((muu.*ws-gammas));
        elseif indic.Sun==1
            f((q-1)*T+1:T*q)= (chiis).*S.^sigmaas-((ws-gammas));
        elseif indic.Sun==2
            f((q-1)*T+1:T*q)= chiis.*S.^(sigmaas+taul)+gammas-muu.*lambdaa.*(1-taul).*ws.^(1-taul);
            
        end
        q=q+1;
        f((q-1)*T+1:T*q)= gammas.*(upbarS-S);
    end
end

%13- Kuhn Tucker Labour supply
if indic.noskill==0
    %25
    q=q+1;
    f((q-1)*T+1:T*q)= gammalh.*(hh-upbarH);
    %26
    q=q+1;
    f((q-1)*T+1:T*q)= gammall.*(hl-upbarH);
else
    q=q+1;
    f((q-1)*T+1:T*q)= gammalh.*(h-upbarH);
end

% income schedule budget clearing
q=q+1;
f((q-1)*T+1:T*q)= SGov-GovRev*Y; % match gov revenues from calibration 

% emisson limit 
if indic.limit_LF==1
    q=q+1;
    f((q-1)*T+1:T*q)=omegaa*F-deltaa-Ems';
end

%fprintf('number equations: %d; number variables %d', q, length(list.test));
end