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


function answ =bgp(pf0, pg0, wh0, wl0, a0, alphaf, alphag,alphan, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, tau0, taul0, sigmaa)

% subsidies against monopolistic competition => the equilibrium is efficient
% except for the externality, 
zetaf= (1-alphaf)/alphaf;
zetag= (1-alphag)/alphag;
zetan= (1-alphan)/alphan;

% skill supply ratio given wages from HH optimality
hlhh = (wl/wh)^((1-taul0)/(sigmaa+taul0));

%PRICES
pfStar0 = theta0*pf0; % on bgp this ratio of prices is required to hold!
pfTil0 = ((pf0+ tau0)^(1-epsf)*deltaff^epsf + (pfStar0 + tau0)^(1-epsf)*deltafo^epsf)^(1/(1-epsf)); 
pe0 = (pfTil0^(1-epse)*deltaef^epse + pg0^(1-epse)*deltaeg^epse)^(1/(1-epse));
pn0 = ((1 - pe0^(1-epsy)*deltaye^epsy)/deltayn^epsy)^(1/(1-epsy));

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

% use equations for af0, ag0, an0 and equate to expressions resulting under
% point 4)

% technology in deviation from aggregate growth 
afHat0 = af0/a0; agHat0 = ag0/a0; anHat0 = an0/a0;


IN THE FOLLOWING REPLACE DETRENDED GROWTH RATES

sgf0 = (ag0/af0)^(phi/eta)*(rhonf/rhong); % follows from LOM of technology
snf0 = (an0/af0)^(phi/eta)*(rhonf); 
sf0 = hs0*zs/(1 + sgf0 + snf0);
sg0 = sf0*sgf0; 
sn0 = sf0*snf0;
gamma = n0/(sn0^eta*(a0/an0)^phi); % n is growth rate in non-energy technology: n0=An'/An-1

%WORKERS
%lf/lg follows from 1) dividing production functions of F/G
%                   2) substitute F/G from demand side: using optimality
%                   conditions of energy producers and fossil composite
%                   producers => lf/lg as a function of prices and
%                   corrective taxes!
%   => these are demanded ratios in equilibrium

%first output ratios
GF = (pfTil0/pg0)^(epse)*(deltaff+(deltafo)*((pf0+tau0)/(pfStar0+tau0)*(deltafo)/deltaff)^((epsf-1)))^(epsf/(epsf-1));
GN =   (pn0/pe0*deltaye/deltayn)^(epsy)*((pg0/pfTil0)^(epse-1)+1)^(epse/(epse-1)); 
% helper
helpn = (alphan^2*pn0)^(alphan/(1-alphan))*an0;
helpg = (pg0*alphag^2)^(alphag/(1-alphag))*ag0;
helpf = (alphaf^2*pf0)^(alphaf/(1-alphaf))*an0; 

lgf0 = helpf/helpg *GF;
lgn0 = helpn/helpg*GN;

    
lg0 = L/(1 + 1/lgn0 + 1/lgf0); 
lf0 = lg0/lgf0; 
ln0 = lg0/lgn0; 

% from here on there seems to be a detrending
%MACHINES
xfHat0 = (alphaf*pf0/pfix0)^(1/(1-alphaf))*lf0*afHat0;
xgHat0= (alphag*pg0/pgix0)^(1/(1-alphag))*lg0*agHat0;
xnHat0 = (alphan*pn0/pnix0)^(1/(1-alphan))*ln0*anHat0;
 
%SCIENTISTS WAGES
wsfHat0 = eta*gamma*rhonf^(eta)*(1-alphaf)/alphaf*xfHat0/(sf0^(1-eta)*afHat0^phi*(1+n0));
wsgHat0 = eta*gamma*rhong^(eta)*(1-alphag)/alphag*xgHat0/(sg0^(1-eta)*agHat0^phi*(1+n0));
wsnHat0 = eta*gamma*(1-alphan)/alphan*xnHat0/(sn0^(1-eta)*anHat0^phi*(1+n0));

%SUPPLIES
fsHat0 = xfHat0^(alphaf)*afHat0^(1-alphaf)*lf0^(1-alphaf);
gsHat0 = xgHat0^(alphag)*agHat0^(1-alphag)*lg0^(1-alphag);
nsHat0 = xnHat0^(alphan)*anHat0^(1-alphan)*ln0^(1-alphan);

% DEMANDS
term1 = (pn0/pe0)^epsy*(deltaye/deltayn)^epsy; 
term2 = (deltaeg + (pg0/pfTil0)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1));
term3  = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*term1;
term4= (deltaff + ((pf0 + tau0)/(pfStar0 +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));

% helpers skill market clearing
FN=GN/GF;
LMHHg= (thetaf*(1-alphaf)*pf0*1/GF+thetag*(1-alphag)*pg0+thetan*(1-alphan)*pn0/GN);
LMHHn= (thetaf*(1-alphaf)*pf0*FN+thetag*(1-alphag)*pg0*GN+thetan*(1-alphan)*pn0);
LMHHf= (thetaf*(1-alphaf)*pf0+thetag*(1-alphag)*pg0*GF+thetan*(1-alphan)*pn0/FN);

pcons = pn0 + pg0*term1/term2 + pf0*term3/(term2*term4);

ndHat0 = 1/pcons*(pf0 *fsHat0 + pg0*gsHat0 + pn0*nsHat0); 
edHat0 = ndHat0*(pn0/pe0)^epsy*(deltaye/deltayn)^epsy;

gdHat0 = edHat0/term2; 
fdHatTil0 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*gdHat0;

fdHat0 = fdHatTil0/term4; 
d2 = (gsHat0-gdHat0)./gsHat0;  %gives pf and pg.  know pfStar0 so  can find pn0.
d3 = (fsHat0 - fdHat0)./fsHat0;
d4 = hlhh- wh0*(1/(LMHHg*helpg*(thetag/(1-thetag)*wl0/wh0)^thetag)+ 1/(LMHHn*helpn*(thetan/(1-thetan)*wl0/wh0)^thetan) ...
        + 1/(LMHHf*helpf*(thetaf/(1-thetaf)*wl0/wh0)^thetaf)); % follows from skill market clearings
answ = [d2, d3];


