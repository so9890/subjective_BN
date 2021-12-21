%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: bgp.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes the excess supplies of green energy and fossil energy
% given guesses for pf and pg along the balanced growth path
% Inputs: 
% Guess for pf and pg: pf0, pg0
% Model parameters: alphaf, alphag,alpham, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
%    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, L, S, tau0
% Outputs: 
% Excess supply of green and fossil energy (1X2 vector) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function answ =bgp(pf0, pg0, a0, alphaf, alphag,alpham, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, L, S, tau0)

%PRICES
pfStar0 = theta0*pf0;
pfTil0 = ((pf0+ tau0)^(1-epsf)*deltaff^epsf + (pfStar0 + tau0)^(1-epsf)*deltafo^epsf)^(1/(1-epsf)); 
pe0 = (pfTil0^(1-epse)*deltaef^epse + pg0^(1-epse)*deltaeg^epse)^(1/(1-epse));
pm0 = ((1 - pe0^(1-epsy)*deltaye^epsy)/deltayn^epsy)^(1/(1-epsy));

pfix0 = 1/(alphaf);
pgix0 = 1/(alphag); 
pmix0 = 1/(alpham); 

%TECHONOLOGIES 
%Labor market clearing implies:
afg0 = (1-alphag)/(1-alphaf) *pg0^(1/(1-alphag))/pf0^(1/(1-alphaf))*(alphag/pgix0)^(alphag/(1-alphag))*(pfix0/alphaf)^(alphaf/(1-alphaf)); 
amg0 = (1-alphag)/(1-alpham) *pg0^(1/(1-alphag))/pm0^(1/(1-alpham))*(alphag/pgix0)^(alphag/(1-alphag))*(pmix0/alpham)^(alpham/(1-alpham)); 

ag0 = (1 + 1/rhonf + 1/rhong)*a0/(1/rhong + 1/rhonf*afg0 + amg0);
af0 = ag0*afg0; 
am0 = ag0*amg0;  

afHat0 = af0/a0; agHat0 = ag0/a0; amHat0 = am0/a0;

sgf0 = (ag0/af0)^(phi/eta)*(rhonf/rhong); 
smf0 = (am0/af0)^(phi/eta)*(rhonf); 
sf0 = S/(1 + sgf0 + smf0);
sg0 = sf0*sgf0; 
sm0 = sf0*smf0;
gamma = n0/(sm0^eta*(a0/am0)^phi); 

%WORKERS
%Scientist market clearing implies (wsg = wsm)
lgm0 = (sg0^(1-eta)/sm0^(1-eta))*pm0^(1/(1-alpham))/pg0^(1/(1-alphag))*(am0/ag0)^(1-phi)*(1-alpham)/(1-alphag)*eta/eta...
    *pgix0^(alphag/(1-alphag))/pmix0^(alpham/(1-alpham))*alpham^(1/(1-alpham))/alphag^(1/(1-alphag))*(1/rhong)^(eta);

lgf0 = (sg0^(1-eta)/sf0^(1-eta))*pf0^(1/(1-alphaf))/pg0^(1/(1-alphag))*(af0/ag0)^(1-phi)*(1-alphaf)/(1-alphag)*eta/eta...
    *pgix0^(alphag/(1-alphag))/pfix0^(alphaf/(1-alphaf))*alphaf^(1/(1-alphaf))/alphag^(1/(1-alphag))*(rhonf/rhong)^(eta); 

lg0 = L/(1 + 1/lgm0 + 1/lgf0); 
lf0 = lg0/lgf0; 
lm0 = lg0/lgm0; 

%MACHINES
xfHat0 = (alphaf*pf0/pfix0)^(1/(1-alphaf))*lf0*afHat0;
xgHat0= (alphag*pg0/pgix0)^(1/(1-alphag))*lg0*agHat0;
xmHat0 = (alpham*pm0/pmix0)^(1/(1-alpham))*lm0*amHat0;
 
%SCIENTISTS WAGES
wsfHat0 = eta*gamma*rhonf^(eta)*(1-alphaf)/alphaf*xfHat0/(sf0^(1-eta)*afHat0^phi*(1+n0));
wsgHat0 = eta*gamma*rhong^(eta)*(1-alphag)/alphag*xgHat0/(sg0^(1-eta)*agHat0^phi*(1+n0));
wsmHat0 = eta*gamma*(1-alpham)/alpham*xmHat0/(sm0^(1-eta)*amHat0^phi*(1+n0));

%SUPPLIES
fsHat0 = xfHat0^(alphaf)*afHat0^(1-alphaf)*lf0^(1-alphaf);
gsHat0 = xgHat0^(alphag)*agHat0^(1-alphag)*lg0^(1-alphag);
msHat0 = xmHat0^(alpham)*amHat0^(1-alpham)*lm0^(1-alpham);

% DEMANDS
term1 = (pm0/pe0)^epsy*(deltaye/deltayn)^epsy; 
term2 = (deltaeg + (pg0/pfTil0)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1));
term3  = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*term1;
term4= (deltaff + ((pf0 + tau0)/(pfStar0 +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));

pcons = pm0 + pg0*term1/term2 + pf0*term3/(term2*term4);

mdHat0 = 1/pcons*(pf0 *fsHat0 + pg0*gsHat0 + pm0*msHat0); 
edHat0 = mdHat0*(pm0/pe0)^epsy*(deltaye/deltayn)^epsy;

gdHat0 = edHat0/term2; 
fdHatTil0 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*gdHat0;

fdHat0 = fdHatTil0/term4; 
d2 = (gsHat0-gdHat0)./gsHat0;  %gives pf and pg.  know pfStar0 so  can find pm0.
d3 = (fsHat0 - fdHat0)./fsHat0;

answ = [d2, d3];


