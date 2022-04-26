function [c, ceq] = constraints(x, T, params, init, list, Ems)
% function to read in constraints on government problem
% there is no order of constraints

%% read in stuff
% choice variables
% hhf, hhg, => replace hhn by market clearing
% hlf, hlg  => hln as above
% C, F, G : intermediate outputs
% Af, Ag, An : technology
% hl, hh

 hhf    = x((find(list.opt=='hhf')-1)*T+1:find(list.opt=='hhf')*T);
 hhg    = x((find(list.opt=='hhg')-1)*T+1:(find(list.opt=='hhg'))*T);
 hlf    = x((find(list.opt=='hlf')-1)*T+1:find(list.opt=='hlf')*T);
 hlg    = x((find(list.opt=='hlg')-1)*T+1:find(list.opt=='hlg')*T);
 C      = x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T);
 F      = x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T);
 G      = x((find(list.opt=='G')-1)*T+1:find(list.opt=='G')*T);
 Af     = x((find(list.opt=='Af')-1)*T+1:find(list.opt=='Af')*T);
 Ag     = x((find(list.opt=='Ag')-1)*T+1:find(list.opt=='Ag')*T);
 An     = x((find(list.opt=='An')-1)*T+1:find(list.opt=='An')*T);
 hl     = x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T);
 hh     = x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T);

 % parameters
 read_in_params;


%% auxiliary variables

hhn     = zh*hh-(hhf+hhg);
hln     = (1-zh)*hl-(hlf+hlg); 
hhhl    = hh./hl;
Lg      = hhg.^thetag.*hlg.^(1-thetag);
Ln      = hhn.^thetan.*hln.^(1-thetan);
Lf      = hhf.^thetaf.*hlf.^(1-thetaf);
%A       = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);
Af_lag  = [Af0;Af(1:T-1)]; % shift Af backwards
Ag_lag  = [Ag0;Ag(1:T-1)];
An_lag  = [An0;An(1:T-1)];
%A_lag   = [max([Af0,Ag0,An0]);A(1:T-1)];
A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);


muu      = C.^(-thetaa); % same equation in case thetaa == 1
% Muhh    = -zh*hh.^(-sigmaa);
% Muhl    = -(1-zh)*hl.^(-sigmaa);
% dIdhh   = lambdaa
% dIdhl 
% scientists follow from LOW fossil and green 
sff     = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
sn      =  S-(sff+sg); 

% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
%Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% wages and policy elements

tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); % production fossil
wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf))./(Af./Af_lag*rhof^etaa); 
wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan))./(An./An_lag*rhon^etaa);
wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag))./(Ag./Ag_lag*rhog^etaa);  % to include taus

taus    = 1-wsgtil./wsf;
wsg     = wsgtil./(1-taus); % this ensures wsg=ws

% assuming interior solution households
taul    = (exp(wh./wl)-sigmaa*exp(hhhl))./(exp(hhhl)+exp(wh./wl)); % from equating FOCs wrt skill supply, solve for taul
lambdaa = hl.^(sigmaa+taul)./(muu.*(1-taul).*wl.^(1-taul));     % from FOC on low skill supply (assuming interior solution)

% auxiliary stuff depending on prices
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
            +tauf.*pf.*F;
        % subsidies, profits and wages scientists cancel

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;


xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = zeros(); %  periods and 2 additional ones: 
             %  IMP; skill time endowment

             % the constraints valid for T periods are:

             % 1) budget constraint HH (IMP)
                         

% 4) emission targets coded as linear constraint!
%%% 1. Implementability constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- create discount vector
        % disc=repmat(betaa, 1,T);
        % expp=0:T-1;
        % vec_discount= disc.^expp;
        % 
        % c(2*T+1) =vec_discount*((mu.*C)-...
        %     (((1-zh).*hl.^(sigmaa+1)+zh.*hh.^(sigmaa+1))./(1-taul)+mu.*SGov)); %

% rather without savings technology: one binding constraint per period
c(1:T)= C-zh.*wh.*hh-(1-zh).*wl.*hl-tauf.*pf.*F;
% (muu.*C)-(((1-zh).*hl.^(sigmaa+1)+zh.*hh.^(sigmaa+1))./(1-taul)+muu.*SGov);

%%% 2./3. Time endowment labour %%%
% => can drop kuhn tuckker constraint! in Household problem, the government
% always chooses an interior solution! 
% c(T*1+1:T*2)=hh-upbarH;
% c(T*2+1:T*3)=hl-upbarH;

%%% 4. Emission constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if indic.target==1
% c(T+1:T*2) = omegaa*F-(Ems'+deltaa);
% end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];

 ceq(1:T)       = SGov; % balanced budget
 ceq(T*2+1:T*3) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
 % optimality skills (for fossil used to determined wage rates)
 ceq(T*3+1:T*4) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
 ceq(T*4+1:T*5) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
 ceq(T*5+1:T*6) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
 ceq(T*6+1:T*7) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
 ceq(T*7+1:T*8) = wsf-wsn; % wage scientists neutral
 ceq(T*8+1:T*9) = wsf-wsg; % wage scientists green
 ceq(T*9+1:T*10) = An_lag.*(1+gammaa*(sn./rhon).^etaa.*(A_lag./An_lag).^phii)-An; % LOM neutral technology 
 %ceq(T*10+1:T*11) = C+xf+xn+xg-Y; % final good market clearing Should be
 %superfluous as other markets clear

ceq = ceq';
end