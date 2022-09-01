function [xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag,sff, sg,  ...
            F, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pg, pf, pee,  ws, wsf, wsg,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom_nn(x, list, params, T, init201519, indic)

read_in_params;

% lambdaa = x((find(list.opt=='lambdaa')-1)*T+1:(find(list.opt=='lambdaa'))*T);
Lf     = x((find(list.opt=='Lf')-1)*T+1:find(list.opt=='Lf')*T);
Lg     = x((find(list.opt=='Lg')-1)*T+1:find(list.opt=='Lg')*T);
C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
h      = x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T);
if indic.xgrowth==0
    sg      = x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T);
    sff      = x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T);
%     gammasg      = x((find(list.opt=='gammasg')-1)*T+1:find(list.opt=='gammasg')*T);
%     gammasf      = x((find(list.opt=='gammasf')-1)*T+1:find(list.opt=='gammasf')*T);
end
%% auxiliary variables
% loop over technology
if indic.xgrowth==0
    Af=zeros(T,1);
    Af_lag=[init201519(list.init=='Af0'); Af(1:T)]; % drop last value later
    Ag=zeros(T,1);
    Ag_lag=[init201519(list.init=='Ag0'); Ag(1:T)];
    A_lag=zeros(T,1);


    for i=1:T
        A_lag(i)   = (rhof*Af_lag(i)+rhog*Ag_lag(i))./(rhof+rhog);

        Af(i)=Af_lag(i).*(1+gammaa*(sff(i)/rhof).^etaa.*(A_lag(i)/Af_lag(i))^phii);
        Ag(i)=Ag_lag(i).*(1+gammaa*(sg(i)/rhog).^etaa.*(A_lag(i)/Ag_lag(i))^phii);

        %-update lags

        Af_lag(i+1)=Af(i);
        Ag_lag(i+1)=Ag(i);

    end

    Af_lag=Af_lag(1:end-1);
    Ag_lag=Ag_lag(1:end-1);
else
    % initial values
    Ag_lag=init201519(list.init=='Ag0');
    Af_lag=init201519(list.init=='Af0');

    Af=zeros(size(F));
    Ag=zeros(size(F));

    for i=1:T
        Ag(i)=(1+vg)*Ag_lag;
        Af(i)=(1+vf)*Af_lag;
        %- update laggs
        Af_lag=Af(i);
        Ag_lag=Ag(i);
    end
    S       = zeros(size(F));
    sff     = zeros(size(F));
    sg      = zeros(size(F));
    A_lag   = (rhof*Af_lag+rhog*Ag_lag)./(rhof+rhog);

end

muu = C.^(-thetaa); % same equation in case thetaa == 1
% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
Y       = E; % (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% wages and policy elements
tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); % production fossil
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    
%- wages scientists  
if indic.xgrowth==0
    wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
    wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus
else
    wsf=zeros(size(sff));
    wsg=zeros(size(sff));
end
%- relevant for code without separate markets
S    = sg+sff;
ws   = chiis*S.^sigmaas; 

% assuming interior solution households:

if indic.notaul==0
    if thetaa~=1
        taul0 = 0.2*ones(size(sg));
        lambdaa0=ones(size(sg));
        ff=@(x)[w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)+tauf.*pf.*F;% balanced budget gov.
                chii*h.^(sigmaa+x(1:T))-(muu.*x(T+1:2*T).*(1-x(1:T)).*(w).^(1-x(1:T)))];
       optionsfs = optimoptions('fsolve', 'TolFun', 10e-12,'Display','none');% 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

        soll = fsolve(ff, [taul0; lambdaa0], optionsfs); 
        taul=soll(1:T);
        lambdaa=soll(T+1:T*2);
    else
        taul    = 1-chii.*h.^(sigmaa+1);
        lambdaa = (w.*h+tauf.*pf.*F)./(w.*h).^(1-taul);  
    end
else
    taul=zeros(size(sg));
    lambdaa=tauf.*pf.*F./(w.*h)+1;
end

xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
        
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhog*Ag)/(rhof+rhog);

gammac =(1+gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii)-1;
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

Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas);
end


 SWF = Utilcon-Utillab-Utilsci;

end