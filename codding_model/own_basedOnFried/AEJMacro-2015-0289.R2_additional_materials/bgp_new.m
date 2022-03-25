%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: bgp.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes the excess supplies of green energy and fossil energy
% given guesses for pf and pg along the balanced growth path
% Inputs: 
% Guess for pf and pg: pf0, pg0
% Model parameters: alphaf, alphag,alphan, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
%    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, L, S, tau0
% Outputs: 
% Excess supply of green and fossil energy (1X2 vector) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function answ =bgp_new(pf0, pg0, wh0, wl0, zl, zh, zs, a0, alphaf, alphag,alphan, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, tau0, taul0, sigmaa)

% subsidies against monopolistic competition => the equilibrium is efficient
% except for the externality, 
zetaf= (1-alphaf)/alphaf;
zetag= (1-alphag)/alphag;
zetan= (1-alphan)/alphan;

% Note that aggregate skill supply => displays a trend
%S  = hs*zs;
%Hh = hh*zh;
%Hl = hl*zl;

% skill supply ratio given wages from HH optimality
hlhh0 = (wl0/wh0)^((1-taul0)/(sigmaa+taul0));

% using labour firms optimal input shares it follows
hhfhlf0 = thetaf/(1-thetaf)*wl0/wh0;
hhnhln0 = thetan/(1-thetan)*wl0/wh0;
hhghlg0 = thetag/(1-thetag)*wl0/wh0;

% this gives the labour input good normalised

%PRICES
pfStar0 = theta0*pf0; % on bgp this ratio of prices is required to hold!
pfTil0  = ((pf0+ tau0)^(1-epsf)*deltaff^epsf + (pfStar0 + tau0)^(1-epsf)*deltafo^epsf)^(1/(1-epsf)); 
pe0     = (pfTil0^(1-epse)*deltaef^epse + pg0^(1-epse)*deltaeg^epse)^(1/(1-epse));
pn0     = ((1 - pe0^(1-epsy)*deltaye^epsy)/deltayn^epsy)^(1/(1-epsy));

% due to the subsidy, monopolists produce at the efficient price
pfix0 = 1/(alphaf*(1+zetaf));
pgix0 = 1/(alphag*(1+zetag)); 
pnix0 = 1/(alphan*(1+zetan)); 

%TECHONOLOGIES 
% follows from 1) optimality labour producers, substitute price for labour input good, 
%              2) substitute production function by intermediate good,
%              3) in production function of intermediate good suubstitute
%              labour input good by production function and optimal input
%              skills derived in 2)
%              4) => defines productivity as a function of wh0, wl0, pf,
%              pn, pg
%              5) below results from division

afg0 = (1-alphag)/(1-alphaf) *pg0^(1/(1-alphag))/pf0^(1/(1-alphaf))*(alphag/alphaf)^2 ...
        *((1-thetag)/wl0)^(1-thetag)*((thetag)/wh0)^(thetag)/(((1-thetaf)/wl0)^(1-thetaf)*((thetaf)/wh0)^(thetaf)); 
ang0 = (1-alphag)/(1-alphan) *pg0^(1/(1-alphag))/pn0^(1/(1-alphan))*(alphag/alphan)^2 ...
        *((1-thetag)/wl0)^(1-thetag)*((thetag)/wh0)^(thetag)/(((1-thetan)/wl0)^(1-thetan)*((thetan)/wh0)^(thetan));

% technologies in levels!
ag0 = (1 + 1/rhonf + 1/rhong)*a0/(1/rhong + 1/rhonf*afg0 + ang0);
af0 = ag0*afg0; 
an0 = ag0*ang0;  

% use equations for ag0 and explicit one to  and equate to expressions resulting under
% point 4)

% technology in deviation from aggregate growth 
afHat0 = af0/a0; agHat0 = ag0/a0; anHat0 = an0/a0;

% using skill market clearing, production functions of intermediate  I solve for hlg/hh as a
% function of the wage rate, etc. the variable is stable by assumption 

% helpers skill market clearing
%first output ratios
GF = (pfTil0/pg0)^(epse)*(deltaff+(deltafo)*((pf0+tau0)/(pfStar0+tau0)*(deltafo)/deltaff)^((epsf-1)))^(epsf/(epsf-1));
GN = (pn0/pe0*deltaye/deltayn)^(epsy)*((pg0/pfTil0)^(epse-1)+1)^(epse/(epse-1)); 
FN = GN/GF;

LMHHg= (thetaf*(1-alphaf)*pf0*1/GF+thetag*(1-alphag)*pg0+thetan*(1-alphan)*pn0/GN);
LMHHn= (thetaf*(1-alphaf)*pf0*FN+thetag*(1-alphag)*pg0*GN+thetan*(1-alphan)*pn0);
LMHHf= (thetaf*(1-alphaf)*pf0+thetag*(1-alphag)*pg0*GF+thetan*(1-alphan)*pn0/FN);

helpn = (alphan^2*pn0)^(alphan/(1-alphan));
helpg = (pg0*alphag^2)^(alphag/(1-alphag));
helpf = (alphaf^2*pf0)^(alphaf/(1-alphaf)); 

hlfHh0 = wh0/(LMHHf*helpf*afHat0*(thetaf/(1-thetaf)*wl0/wh0)^thetaf); % hlf/(hh*zh);
hlgHh0 = wh0/(LMHHg*helpg*agHat0*(thetag/(1-thetag)*wl0/wh0)^thetag); % hlg/(hh*zh);
hlnHh0 = wh0/(LMHHn*helpn*anHat0*(thetan/(1-thetan)*wl0/wh0)^thetan); % hln/(hh*zh);

% => detrended labour input good
LfHh0 = (hhfhlf0)^thetaf*hlfHh0;
LgHh0 = (hhghlg0)^thetag*hlgHh0;
LnHh0 = (hhnhln0)^thetan*hlnHh0;

% scientists => all scientists are alike
sgf0 = (1/afg0)^(phi/eta)*(rhonf/rhong); % follows from LOM of technology
snf0 = (afg0/ang0)^(phi/eta)*(rhonf);
% detrend the level of scientist=> with falling hours there is a trend in
% labour supply, but assume shares to be constant

sfHat0 = 1/(1 + sgf0 + snf0); % divide by total scientist supply (hs*zs) which moves
sgHat0 = sfHat0*sgf0; 
snHat0 = sfHat0*snf0;
gamma = n0/(snHat0^eta*(a0/an0)^phi); % n is growth rate in non-energy technology: n0=An'/An-1

%Labour input good by sector

%lf/lg follows from 1) dividing production functions of F/G
%                   2) substitute F/G from demand side: using optimality
%                   conditions of energy producers and fossil composite
%                   producers => lf/lg as a function of prices and
%                   corrective taxes!
%   => these are demanded ratios in equilibrium


% another definition of labour input ratios from the demand side => has to
% be equal in equilibrium 
lgf0 = helpf/helpg*afg0*GF;
lgn0 = helpn/helpg*ang0*GN;

% excess labour market demand
d1 = lgf0 - (LgHh0/LfHh0);
d2 = lgn0 - (LgHh0/LgHh0); 

%MACHINES
xfHat0 = (alphaf*pf0/pfix0)^(1/(1-alphaf))*LfHh0*afHat0;
xgHat0= (alphag*pg0/pgix0)^(1/(1-alphag))*LgHh0*agHat0;
xnHat0 = (alphan*pn0/pnix0)^(1/(1-alphan))*LnHh0*anHat0;
 
%SCIENTISTS WAGES
wsfHat0 = eta*gamma*rhonf^(eta)*(1-alphaf)/alphaf*xfHat0/(sfHat0^(1-eta)*afHat0^phi*(1+n0)); %=> wsfHat= wsf/(A*Hh*Hs^(eta-1)/rhon^eta)
wsgHat0 = eta*gamma*rhong^(eta)*(1-alphag)/alphag*xgHat0/(sgHat0^(1-eta)*agHat0^phi*(1+n0)); % this detrends machines and scientist supply
wsnHat0 = eta*gamma*(1-alphan)/alphan*xnHat0/(snHat0^(1-eta)*anHat0^phi*(1+n0)); % division by (1/rhon)^eta => rhon^eta cancels

%SUPPLIES
fsHat0 = (alphaf^2*pf0)^(alphaf/(1-alphaf))*afHat0*LfHh0;
gsHat0 = (alphag^2*pg0)^(alphag/(1-alphag))*agHat0*LgHh0;
nsHat0 = (alphan^2*pn0)^(alphan/(1-alphan))*anHat0*LnHh0; %=> constant: nshat0= ns/(A*Hh )

% DEMANDS=> follow from final good producers
term1 = (pn0/pe0)^epsy*(deltaye/deltayn)^epsy; % optimal input ratio Energy and non-Energy goods in final good production 
term2 = (deltaeg + (pg0/pfTil0)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1)); % production and optimal input energy sector
term3 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*term1; % optimality energy sector times optimal ratio energy/non-energy
term4 = (deltaff + ((pf0 + tau0)/(pfStar0 +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));

pcons = pn0 + pg0*term1/term2 + pf0*term3/(term2*term4);

% assume market for non-energy good clears, what does it take to have the
% others clear?

ndHat0 = 1/pcons*(pf0 *fsHat0 + pg0*gsHat0 + pn0*nsHat0); 
edHat0 = ndHat0*(pn0/pe0)^epsy*(deltaye/deltayn)^epsy;

gdHat0 = edHat0/term2; 
fdHatTil0 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*gdHat0;

fdHat0 = fdHatTil0/term4; 
d3 = (gsHat0-gdHat0)./gsHat0;  %gives pf and pg.  know pfStar0 so  can find pn0.
d4 = (fsHat0 - fdHat0)./fsHat0;
% d4 = hlhh0- wh0*zh*(1/(LMHHg*helpg*(thetag/(1-thetag)*wl0/wh0)^thetag)+ 1/(LMHHn*helpn*(thetan/(1-thetan)*wl0/wh0)^thetan) ...
%         + 1/(LMHHf*helpf*(thetaf/(1-thetaf)*wl0/wh0)^thetaf)); % follows from skill market clearings
% d5 = ag0-1/(ag0^2*pg0^ (1/(1-alphag))*(1-alphag)*((1-thetag)/wl0)^(1-thetag)*(thetag/wh0)^thetag); % follows from aggregate productivity definition 

answ = [d2, d3, d4, d5];


