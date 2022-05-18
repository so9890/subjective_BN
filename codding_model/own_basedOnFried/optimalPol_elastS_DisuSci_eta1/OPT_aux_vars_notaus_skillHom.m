function [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init201519, indic)

read_in_params;

% lambdaa = x((find(list.opt=='lambdaa')-1)*T+1:(find(list.opt=='lambdaa'))*T);
Lf     = x((find(list.opt=='Lf')-1)*T+1:find(list.opt=='Lf')*T);
Lg     = x((find(list.opt=='Lg')-1)*T+1:find(list.opt=='Lg')*T);
C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
Af     = x((find(list.opt=='Af')-1)*T+1:find(list.opt=='Af')*T);
Ag     = x((find(list.opt=='Ag')-1)*T+1:find(list.opt=='Ag')*T);
An     = x((find(list.opt=='An')-1)*T+1:find(list.opt=='An')*T);
h      = x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T);  
S      = x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T);
 
% initial values: CALIBRATED dont change
An0=init201519(list.init=='An0');
Ag0=init201519(list.init=='Ag0');
Af0=init201519(list.init=='Af0');

%% auxiliary variables

Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0;Ag(1:T-1)];
An_lag  = [An0;An(1:T-1)];
%A_lag   = [max([Af0,Ag0,An0]);A(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);


muu = C.^(-thetaa); % same equation in case thetaa == 1
% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

Ln      = N./(An.*(pn.*alphan).^(alphan./(1-alphan))); % production neutral

% wages and policy elements
tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); % production fossil
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    
ws      = chiis*S.^sigmaas./muu; 
sff     = ((gammaa*etaa*(A_lag./Af_lag).^phii.*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(ws.*Af.*rhof^etaa)).^(1/(1-etaa));
sg      = ((gammaa*etaa*(A_lag./Ag_lag).^phii.*pg.*G*(1-alphag).*Ag_lag)./(ws.*Ag.*rhog^etaa)).^(1/(1-etaa));
sn      = ((gammaa*etaa*(A_lag./An_lag).^phii.*pn.*N*(1-alphan).*An_lag)./(ws.*An.*rhon^etaa)).^(1/(1-etaa));
% sg      = S -(sff+sn);

%wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

%taus    = 1-wsgtil./ws; % since (wsgtilde/(1-taus)=ws)


% assuming interior solution households
if indic.notaul==0
    taul0 = 0.2*ones(size(sn));
    lambdaa0=ones(size(sn));
    ff=@(x)[w.*h-(w.*h).^(1-x(1:T)).*x(T+1:2*T)+tauf.*pf.*F;% balanced budget gov.
            chii*h.^(sigmaa+x(1:T))-(muu.*x(T+1:2*T).*(1-x(1:T)).*(w).^(1-x(1:T)))];
   optionsfs = optimoptions('fsolve', 'TolFun', 10e-12,'Display','none');% 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

    soll = fsolve(ff, [taul0; lambdaa0], optionsfs); 
    taul=soll(1:T);
    lambdaa=soll(T+1:T*2);
else
    taul=zeros(size(sn));
    lambdaa=tauf.*pf.*F./(w.*h)+1;
end

xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;
        
Emnet     = omegaa*F-deltaa; % net emissions
A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

% wh2 = thetag*(hhg./hlg).^(thetag-1).*wlg;
% utility
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(h.^(1+sigmaa))./(1+sigmaa);
 Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);

 SWF = Utilcon-Utillab-Utilsci;

end