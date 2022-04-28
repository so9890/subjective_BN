% get objective function of ramsey planner
function [OB_RAM, x, list, symms]= OPT_aux_vars( list, params, T, init, indic)

% prepare variables
syms    mu_IMP mu_MarketS mu_wageG mu_wageN mu_NProd mu_OPThhn...
        mu_OPThhg mu_OPThln mu_OPThlg ...
        hhf hhg hlf hlg C F G Af Ag An ...
        hl hh sn sff sg mu_target real
symms.optsym=[ mu_IMP mu_MarketS mu_wageG mu_wageN mu_NProd mu_OPThhn...
        mu_OPThhg mu_OPThln mu_OPThlg hhf hhg hlf hlg C F G Af Ag An hl hh];
if indic.target== 1
    symms.optsym=[symms.optsym mu_target];
end
list.optsym=string(symms.optsym);
    
% create symbolic variables
vecs=sym('a',[T,length([list.optsym])]);
    for s = [list.optsym]  % loop over list entries
        vecs(:,[list.optsym]==s)=sym(sprintf('%s%d',s),  [T,1]); 
    end
% stack into one column vector
x=vecs(:); 

%% read in stuff
% inputs: x is a symbolic vector of choice variables for T periods
%params=symms.params; % by passing a symbolic vector read_in_params gives symbolic pars
read_in_params;

% separate to derive auxiliary variables
 hhf    = x((find(list.optsym=='hhf')-1)*T+1:find(list.optsym=='hhf')*T);
 hhg    = x((find(list.optsym=='hhg')-1)*T+1:(find(list.optsym=='hhg'))*T);
 hlf    = x((find(list.optsym=='hlf')-1)*T+1:find(list.optsym=='hlf')*T);
 hlg    = x((find(list.optsym=='hlg')-1)*T+1:find(list.optsym=='hlg')*T);
 C      = x((find(list.optsym=='C')-1)*T+1:find(list.optsym=='C')*T);
 F      = x((find(list.optsym=='F')-1)*T+1:find(list.optsym=='F')*T);
 G      = x((find(list.optsym=='G')-1)*T+1:find(list.optsym=='G')*T);
 Af     = x((find(list.optsym=='Af')-1)*T+1:find(list.optsym=='Af')*T);
 Ag     = x((find(list.optsym=='Ag')-1)*T+1:find(list.optsym=='Ag')*T);
 An     = x((find(list.optsym=='An')-1)*T+1:find(list.optsym=='An')*T);
 hl     = x((find(list.optsym=='hl')-1)*T+1:find(list.optsym=='hl')*T);
 hh     = x((find(list.optsym=='hh')-1)*T+1:find(list.optsym=='hh')*T);

 % lagrange multiplier
mu_IMP     = x((find(list.optsym=='mu_IMP')-1)*T+1:find(list.optsym=='mu_IMP')*T);
mu_MarketS = x((find(list.optsym=='mu_MarketS')-1)*T+1:find(list.optsym=='mu_MarketS')*T);
mu_NProd   = x((find(list.optsym=='mu_NProd')-1)*T+1:find(list.optsym=='mu_NProd')*T);
mu_OPThhg  = x((find(list.optsym=='mu_OPThhg')-1)*T+1:find(list.optsym=='mu_OPThhg')*T);
mu_OPThhn  = x((find(list.optsym=='mu_OPThhn')-1)*T+1:find(list.optsym=='mu_OPThhn')*T);
mu_OPThlg  = x((find(list.optsym=='mu_OPThlg')-1)*T+1:find(list.optsym=='mu_OPThlg')*T);
mu_OPThln  = x((find(list.optsym=='mu_OPThln')-1)*T+1:find(list.optsym=='mu_OPThln')*T);
mu_wageG   = x((find(list.optsym=='mu_wageG')-1)*T+1:find(list.optsym=='mu_wageG')*T);
mu_wageN   = x((find(list.optsym=='mu_wageN')-1)*T+1:find(list.optsym=='mu_wageN')*T);

% initial values: CALIBRATED dont change
An0=init(list.init=='An0');
Ag0=init(list.init=='Ag0');
Af0=init(list.init=='Af0');

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

% scientists follow from LOW fossil and green 
sff     = ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
sg      = ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
sn      = ((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);

% prices
pg      = (G./(Ag.*Lg)).^((1-alphag)/alphag)./alphag; % from production function green
pf      = (G./F).^(1/eppse).*pg; % optimality energy producers
pee     = (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse));
pn      = ((1-deltay.*pee.^(1-eppsy))./(1-deltay)).^(1/(1-eppsy)); % definition prices and numeraire

% output
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));
N       = (1-deltay)/deltay.*(pee./pn).^(eppsy).*E; % demand N final good producers 
Y       = (deltay^(1/eppsy)*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy)*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % final output production

% wages and policy elements

tauf    = 1-(F./(Af.*Lf)).^((1-alphaf)/alphaf)./(alphaf*pf); % production fossil
wh      = thetaf*(hlf./hhf).^(1-thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af; % from optimality labour input producers fossil, and demand labour fossil
wl      = (1-thetaf)*(hhf./hlf).^(thetaf).*(1-alphaf).*alphaf^(alphaf/(1-alphaf)).*...
        ((1-tauf).*pf).^(1/(1-alphaf)).*Af;
wsf     = (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*F*(1-alphaf).*(1-tauf).*Af_lag)./(Af.*rhof^etaa); 
wsn     = (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N*(1-alphan).*An_lag)./(An.*rhon^etaa);
wsgtil  = (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G*(1-alphag).*Ag_lag)./(Ag.*rhog^etaa);  % to include taus

taus    = 1-wsgtil./wsf;
wsg     = wsgtil./(1-taus); % this ensures wsg=ws

% assuming interior solution households
taul    = (exp(wh./wl)-sigmaa*exp(hhhl))./(exp(hhhl)+exp(wh./wl)); % from equating FOCs wrt skill supply, solve for taul

% lambdaa so that gov budget is balanced
lambdaa = (zh*(wh.*hh)+(1-zh)*(wl.*hl)+tauf.*pf.*F)./...
            (zh*(wh.*hh).^(1-taul)+(1-zh)*(wl.*hl).^(1-taul)); 
        % subsidies, profits and wages scientists cancel

wln     = pn.^(1/(1-alphan)).*(1-alphan).*alphan.^(alphan/(1-alphan)).*An; % price labour input neutral sector
wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
% wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 
% 
% xn      = (alphan*pn).^(1/(1-alphan)).*Ln.*An;
% xf      = (alphaf*pf.*(1-tauf)).^(1/(1-alphaf)).*Lf.*Af;
% xg      = (alphag*pg).^(1/(1-alphag)).*Lg.*Ag;
% 
% SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
%             +(1-zh)*(wl.*hl-lambdaa.*(wl.*hl).^(1-taul))...
%             +tauf.*pf.*F;
%         
% Emnet     = omegaa*F-deltaa; % net emissions
% A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);


%% social welfare

%- create discount vector
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;

%- vector of utilities
if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);

 %% constraints
IMP     = C-zh.*wh.*hh-(1-zh).*wl.*hl-tauf.*pf.*F; % one each period
MarketS = (sff+sg+sn)-S;  % LOM neutral technology 
wageG   = wsf-wsg; % wage scientists green
wageN   = wsf-wsn; % wage scientists neutral
NProd   = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
OPThhn  = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
OPThhg  = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
OPThln  = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
OPThlg  = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low

if indic.target==1
    mu_target
 %% objective function 
OB_RAM=vec_discount*(Utilcon-Utillab...
        - mu_IMP.*IMP- mu_MarketS.*MarketS...
        - mu_wageG.*wageG- mu_wageN.*wageN...
        - mu_NProd.*NProd-mu_OPThhn.*OPThhn...
        - mu_OPThhg.*OPThhg -mu_OPThln.*OPThln...
        - mu_OPThlg.*OPThlg);
end