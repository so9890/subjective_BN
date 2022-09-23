function [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, GovCon, Tls, PV,PVSWF, objF]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init201519, indic, MOM)

read_in_params;

% lambdaa = x((find(list.opt=='lambdaa')-1)*T+1:(find(list.opt=='lambdaa'))*T);
Lf     = x((find(list.opt=='Lf')-1)*T+1:find(list.opt=='Lf')*T);
Lg     = x((find(list.opt=='Lg')-1)*T+1:find(list.opt=='Lg')*T);
C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
h      = x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T);
if indic.xgrowth==0
    sn      = x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T);
    sg      = x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T);
    sff      = x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T);
    
 if indic.sep==1
    S    = sn+sg+sff;
 elseif indic.sep==0
    S     = x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T);
 end
else
   sff=sff0*ones(size(F));
    sn=sn0*ones(size(F));
    sg=sg0*ones(size(F));
    S=sn+sff+sg;
end

%% auxiliary variables
% loop over technology
if indic.zero==0
    Af=zeros(T,1);
    Af_lag=[init201519(list.init=='Af0'); Af(1:T)]; % drop last value later
    Ag=zeros(T,1);
    Ag_lag=[init201519(list.init=='Ag0'); Ag(1:T)];
    An=zeros(T,1);
    An_lag=[init201519(list.init=='An0'); An(1:T)];
    A_lag=zeros(T,1);


    for i=1:T
        A_lag(i)   = (rhof*Af_lag(i)+rhon*An_lag(i)+rhog*Ag_lag(i))./(rhof+rhon+rhog);

        Af(i)=Af_lag(i).*(1+gammaa*(sff(i)/rhof).^etaa.*(A_lag(i)/Af_lag(i))^phii);
        Ag(i)=Ag_lag(i).*(1+gammaa*(sg(i)/rhog).^etaa.*(A_lag(i)/Ag_lag(i))^phii);
        An(i)=An_lag(i).*(1+gammaa*(sn(i)/rhon).^etaa.*(A_lag(i)/An_lag(i))^phii);

        %-update lags

        Af_lag(i+1)=Af(i);
        Ag_lag(i+1)=Ag(i);
        An_lag(i+1)=An(i);

    end

    Af_lag=Af_lag(1:end-1);
    An_lag=An_lag(1:end-1);
    Ag_lag=Ag_lag(1:end-1);
else
    % initial values
    An_lag=init201519(list.init=='An0');
    Ag_lag=init201519(list.init=='Ag0');
    Af_lag=init201519(list.init=='Af0');

    Af=zeros(size(F));
    Ag=zeros(size(F));
    An=zeros(size(F));

    for i=1:T
        An(i)=(1+vn)*An_lag;
        Ag(i)=(1+vg)*Ag_lag;
        Af(i)=(1+vf)*Af_lag;
        %- update laggs
        An_lag=An(i);
        Af_lag=Af(i);
        Ag_lag=Ag(i);
    end
    A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)./(rhof+rhon+rhog);

end

muu = C.^(-thetaa); % same equation in case thetaa == 1
% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf); % production fossil

tauf      = (G./F).^(1/eppse).*pg-pf; % optimality energy producers
pee     = ((pf+tauf).^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

Ln      = N./(An.*(pn.*alphan).^(alphan./(1-alphan))); % production neutral

% wages and policy elements
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*(pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    
%- wages scientists  
% if indic.xgrowth==0
    wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*Af_lag)./(Af.*rhof^etaa); 
    wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa); 
    wsg     = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus
% else
%     wsf=zeros(size(sff));
%     wsg=zeros(size(sff));
%     wsn=zeros(size(sff));
% end
%- relevant for code without separate markets
ws   = (chiis*S.^sigmaas)./muu; 

if indic.notaul >=4
    Tls =tauf.*F;
else
    Tls =zeros(size(F));
end

% assuming interior solution households:
if indic.notaul==0 || indic.notaul == 3 || indic.notaul == 4 
    if thetaa~=1
        taul0 = 0.2*ones(size(sn));
        lambdaa0=ones(size(sn));
        if indic.notaul==0
            exxpr = w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)+tauf.*F-GovRev;
        else
            exxpr = w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)-GovRev; % env tax revenues are not redistributed via income tax scheme
        end
        ff=@(x)[exxpr ;% balanced budget gov.
                chii*h.^(sigmaa+x(1:T))-(muu.*x(T+1:2*T).*(1-x(1:T)).*(w).^(1-x(1:T)))];
       optionsfs = optimoptions('fsolve', 'TolFun', 10e-12,'Display','none');% 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

        soll = fsolve(ff, [taul0; lambdaa0], optionsfs); 
        taul=soll(1:T);
        lambdaa=soll(T+1:T*2);
    elseif thetaa ==1
        
        taul    = 1-chii.*h.^(sigmaa+1).*(1+Tls./(w.*h-GovRev));
        if indic.notaul ==0
            lambdaa = (w.*h+tauf.*F-GovRev)./(w.*h).^(1-taul);
        else
            lambdaa = (w.*h-GovRev)./(w.*h).^(1-taul);
        end
   end
    
else % taul cannot be used (indic.notaul == 1, 2, 5)
    taul=zeros(size(sn));
    if indic.notaul==1 
        lambdaa=(tauf.*F-GovRev)./(w.*h)+1;
    elseif indic.notaul == 2 ||  indic.notaul == 5 % env tax revenues not redistributed via income tax scheme; 
        lambdaa=(-GovRev)./(w.*h)+1;
    end
end

xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

if indic.notaul <2  % 
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul)) +tauf.*F;
        GovCon  = zeros(size(F));
else
        SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul));
        if indic.notaul>=4 % lump sum trans
            GovCon = zeros(size(F));
        else
            GovCon = tauf.*F;
        end
end
        
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

% gammac =(1+gammaa.*(sff(T)./rhof).^etaa.*(A(T)./Af(T)).^phii)-1;

% utility
if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end
  

Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
end


 SWF = Utilcon-Utillab-Utilsci;
 %- create discount vector
     disc=repmat(betaa, 1,T);
     expp=0:T-1;
     vec_discount= disc.^expp;
     PVSWF = vec_discount*SWF;
    
    % continuation value
    %- last period growth rate as proxy for future growth rates
    gammay = Y(T)/Y(T-1)-1;
    PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(T);
    PVwork =indic.PVwork*1/(1-betaa)*(Utillab(T)+Utilsci(T)); % this decreases last period work and science 
    PV= betaa^T*(PVconsump-PVwork);

    %Objective function value:
    %!! Dot product!!! so no dot.*
    % f = (-1)*(vec_discount*(Utilcon-Utillab- Utilsci)+PVcontUtil);+
    objF=(vec_discount*(SWF-indic.extern*weightext*(omegaa.*F).^extexpp)+indic.PV*PV);

end