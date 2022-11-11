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


function answ =bgp_new(pf0, pg0, wh0, wl0, sf0, sg0, sn0, gamma0, zl, zh, S, alphaf, alphag,alphan, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
    epse, epsf, epsy,eta, n0, rhon, rhof, theta0, tau0, taul0, sigmaa, lambda0)

% subsidies against monopolistic competition => the equilibrium is efficient

%MISTAKE IN THE FOLLOWING AS HAVE CODED PXIF WITHOUT TAX
% except for the externality, 
zetaf= (1-alphaf)/alphaf;
zetag= (1-alphag)/alphag;
zetan= (1-alphan)/alphan;

% Note that aggregate skill supply => displays a trend
%S  = hs*zs;
%Hh = hh*zh;
%Hl = hl*zl;

% skill supply ratio given wages from HH optimality
HlHh0 = (wl0/wh0)^((1-taul0)/(sigmaa+taul0))*zh/zl;

% labour firms optimal input shares imply
hhfhlf0 = thetaf/(1-thetaf)*wl0/wh0;
hhnhln0 = thetan/(1-thetan)*wl0/wh0;
hhghlg0 = thetag/(1-thetag)*wl0/wh0;

%PRICES:

pfStar0 = theta0*pf0; % having pfStar as a function of fossil price results from assuming a constant markup by foreign producers
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
af0 = (alphaf^(2*alphaf/(1-alphaf))*pf0^(1/(1-alphaf))*(1-alphaf)*(1-thetaf)^(1-thetaf)*thetaf^thetaf*(wl0/wh0)^thetaf*wl0^(-1))^(-1);
ag0 = (alphag^(2*alphag/(1-alphag))*pg0^(1/(1-alphag))*(1-alphag)*(1-thetag)^(1-thetag)*thetag^thetag*(wl0/wh0)^thetag*wl0^(-1))^(-1);
an0 = (alphan^(2*alphan/(1-alphan))*pn0^(1/(1-alphan))*(1-alphan)*(1-thetan)^(1-thetan)*thetan^thetan*(wl0/wh0)^thetan*wl0^(-1))^(-1);

% afg0 = (1-alphag)/(1-alphaf) *pg0^(1/(1-alphag))/pf0^(1/(1-alphaf))*(alphag/alphaf)^2 ...
%         *((1-thetag)/wl0)^(1-thetag)*((thetag)/wh0)^(thetag)/(((1-thetaf)/wl0)^(1-thetaf)*((thetaf)/wh0)^(thetaf)); 
% ang0 = (1-alphag)/(1-alphan) *pg0^(1/(1-alphag))/pn0^(1/(1-alphan))*(alphag/alphan)^2 ...
%         *((1-thetag)/wl0)^(1-thetag)*((thetag)/wh0)^(thetag)/(((1-thetan)/wl0)^(1-thetan)*((thetan)/wh0)^(thetan));

% technologies in levels!
% a0 = max([af0,ag0, an0]);

% use equations for ag0 and explicit one to  and equate to expressions resulting under
% point 4)

% using skill market clearing, production functions of intermediate  I solve for hlg/hh as a
% function of the wage rate, etc. the variable is stable by assumption 

% helpers skill market clearing
%first output ratios=> NOT STABLE WHEN HAVING TO HAVE F AT A CONSTANT LEVEL
%ON THE BGP
% instead look at expenditure ratios

pgGpfF = pg0/pf0*(pfTil0/pg0)^(epse)*(deltaff+(deltafo)*((pf0+tau0)/(pfStar0+tau0)*(deltafo)/deltaff)^((epsf-1)))^(epsf/(epsf-1));
pgGpnN = pg0/pn0*(pn0/pe0*deltaye/deltayn)^(epsy)*((pg0/pfTil0)^(epse-1)+1)^(epse/(epse-1)); 
pfFpnN = pgGpnN/pgGpfF;

LMHHg= (thetaf*(1-alphaf)/pgGpfF+thetag*(1-alphag)+thetan*(1-alphan)/pgGpnN); % division by pgG
LMHHn= (thetaf*(1-alphaf)*pfFpnN+thetag*(1-alphag)*pgGpnN+thetan*(1-alphan)); % division by pnN
LMHHf= (thetaf*(1-alphaf)+thetag*(1-alphag)*pgGpfF+thetan*(1-alphan)/pfFpnN); % division by pfF

helpn = (alphan^2)^(alphan/(1-alphan));
helpg = (alphag^2)^(alphag/(1-alphag));
helpf = (alphaf^2)^(alphaf/(1-alphaf)); 

hlfHh0 = wh0/(LMHHf*helpf*(thetaf/(1-thetaf)*wl0/wh0)^thetaf*af0*pf0^(1/(1-alphaf))); % hlf/(hh*zh); => constant due to constant labour input 
hlgHh0 = wh0/(LMHHg*helpg*(thetag/(1-thetag)*wl0/wh0)^thetag*ag0*pg0^(1/(1-alphag))); % hlg/(hh*zh);
hlnHh0 = wh0/(LMHHn*helpn*(thetan/(1-thetan)*wl0/wh0)^thetan*an0*pn0^(1/(1-alphan))); % hln/(hh*zh);

% => detrended labour input good
LfHh0 = (hhfhlf0)^thetaf*hlfHh0;
LgHh0 = (hhghlg0)^thetag*hlgHh0;
LnHh0 = (hhnhln0)^thetan*hlnHh0;

% scientists => all scientists are alike; ratio of scientists is stable
% from free movement of scientists => eqbm level of scientists
% substitute LOM technology => scientists as function of known stuff

d4 = (sf0/sn0)^(eta-1)-((1+gamma0*sf0^eta*rhof^(-eta))/(1+gamma0*sn0^eta*rhon^(-eta))*alphan*(1-alphan)/(alphaf*(1-alphaf))*rhofn^eta*pfFpnN);
d5 = (sf0/sg0)^(eta-1)-((1+gamma0*sf0^eta*rhof^(-eta))/(1+gamma0*sg0^eta*rhon^(-eta))*alphag*(1-alphag)/(alphaf*(1-alphaf))*rhofg^eta*pgGpnN);
d6 = sn0 - S/(1 + sf0/sn0 + sg0/sn0); % 
d7 = gamma0 - n0*(rhon/sn0)^eta; % n is growth rate in non-energy technology: n0=An'/An-1

% determine growth in other sectors from here?
% Labour input good by sector

%lf/lg follows from 1) dividing production functions of F/G
%                   2) substitute F/G from demand side: using optimality
%                   conditions of energy producers and fossil composite
%                   producers => lf/lg as a function of prices and
%                   corrective taxes!
%   => these are demanded ratios in equilibrium


% another definition of labour ratios from production (demand)
lgf0= pgGpfF*af0*pf0^(1/(1-alphaf))*alphaf^(2*alphaf/(1-alphaf))/(ag0*pg0^(1/(1-alphag))*alphag^(2*alphag/(1-alphag)));
lgn0= pgGpnN*an0*pn0^(1/(1-alphan))*alphan^(2*alphan/(1-alphan))/(ag0*pg0^(1/(1-alphag))*alphag^(2*alphag/(1-alphag)));

% excess labour market demand => pins down high and low wage
d1 = lgf0 - (LgHh0/LfHh0);
d2 = lgn0 - (LgHh0/LgHh0); 

% %MACHINES
xfHat0 = (alphaf*pf0/pfix0)^(1/(1-alphaf))*LfHh0*af0/wh0; % division by wh*Hh
xgHat0= (alphag*pg0/pgix0)^(1/(1-alphag))*LgHh0*ag0/wh0;
xnHat0 = (alphan*pn0/pnix0)^(1/(1-alphan))*LnHh0*an0/wh0;

%SCIENTISTS WAGES
wsHat0 = eta*gamma0*rhon^(-eta)*(1-alphan)/alphan/(1+n0)*S^(eta-1)*xnHat0; %=> wsfHat= wsf/(A*Hh*Hs^(eta-1)/rhon^eta)

% wsfHat0 = eta*gamma0*rhonf^(eta)*(1-alphaf)/alphaf*xfHat0/(sfHat0^(1-eta)*afHat0^phi*(1+n0)); %=> wsfHat= wsf/(A*Hh*Hs^(eta-1)/rhon^eta)
% wsgHat0 = eta*gamma0*rhong^(eta)*(1-alphag)/alphag*xgHat0/(sgHat0^(1-eta)*agHat0^phi*(1+n0)); % this detrends machines and scientist supply
% wsnHat0 = eta*gamma0*(1-alphan)/alphan*xnHat0/(snHat0^(1-eta)*anHat0^phi*(1+n0)); % division by (1/rhon)^eta => rhon^eta cancels

% AGGREGATE OUTPUT MARKET CLEARING C=Y
income = lambda0*S*wsHat0
C = 

%SUPPLIES in expenditures: output by assumption not constant on BGP! 
pffsHat0 = (alphaf^2)^(alphaf/(1-alphaf))*pf0^(1/(1-alphaf))*af0/wh0*LfHh0; % division by HH*wh and times pf0 => constant
pggsHat0 = (alphag^2)^(alphag/(1-alphag))*pg0^(1/(1-alphag))*ag0/wh0*LgHh0;
pnnsHat0 = (alphan^2)^(alphan/(1-alphan))*pn0^(1/(1-alphan))*an0/wh0*LnHh0; 

% resulting aggregate output
y0
% DEMANDS=> follow from final good producers
% term1 = (pn0/pe0)^epsy*(deltaye/deltayn)^epsy; % optimal input ratio Energy and non-Energy goods in final good production 
% term2 = (deltaeg + (pg0/pfTil0)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1)); % production and optimal input energy sector
% term3 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*term1; % optimality energy sector times optimal ratio energy/non-energy
% term4 = (deltaff + ((pf0 + tau0)/(pfStar0 +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));
% 
% pcons = pn0 + pg0*term1/term2 + pf0*term3/(term2*term4);

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


