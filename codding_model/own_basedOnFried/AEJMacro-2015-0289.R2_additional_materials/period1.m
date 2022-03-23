%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: period1.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes FOC and market clearing in the shock period endog-model.
% Inputs: 
% Guesses for pf, pg, lf, lg, sf sg: pf1, pg1, lf1, lg1, sf1, sg1
% Model parameters and values: alphaf,alphag, alpham, deltaef, deltaeg, deltaff, 
%     deltafo, deltaye, deltayn, epse, epsf, epsy, pfStar1, tau, tauStar, L 
% Outputs: 
% First order and market clearing conditions (should all equal zero in eq.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function answ = period1(lf, lg, pf, pg, sf1, sg1, af0, ag0, am0,a0, pfStar, xfHat0, xgHat0, xmHat0,  ...
    alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn, epse, epsf, epsy, eta,...
    gamma, nu, n0,phi,rhonf, rhong, S,tau0, L)

%PRICES
pfTil = ((pf+ tau0)^(1-epsf)*deltaff^epsf + (pfStar + tau0)^(1-epsf)*deltafo^epsf)^(1/(1-epsf)); 
pe = (pfTil^(1-epse)*deltaef^epse + pg^(1-epse)*deltaeg^epse)^(1/(1-epse));
pm = ((1 - pe^(1-epsy)*deltaye^epsy)/deltayn^epsy)^(1/(1-epsy));

% M sector
lm = L - lf - lg; 

%TECHNOLOGY
sm1 = S - sf1 - sg1;
af1 = af0*(1+ gamma*(rhonf)^(eta)*(sf1).^eta.*(a0./(af0)).^phi); 
ag1 = ag0*(1+ gamma*(rhong)^(eta)*(sg1).^eta.*(a0./(ag0)).^phi); 
am1 = am0*(1+ gamma*sm1.^eta.*(a0./am0).^phi); 
a1  = (1/rhonf*af1 + 1/rhong*ag1 + am1)/(1+1/rhonf + 1/rhong);

afHat1 = af1./a1;
agHat1 = ag1./a1;
amHat1 = am1./a1;

afHat0 = af0./a0;
agHat0 = ag0./a0;
amHat0 = am0./a0;

%MACHINES FIXED AT THE BGP VALUES WITH NO CHANGE IN A
xfHat1 = xfHat0*a0*(1+n0)/a1; 
xgHat1 = xgHat0*a0*(1+n0)/a1; 
xmHat1 = xmHat0*a0*(1+n0)/a1;

%GOODS MARKET CLEARING
fsHat = nu*lf.^(1-alphaf).*xfHat1.^(alphaf).*afHat1.^(1-alphaf);
gsHat = lg.^(1-alphag).*xgHat1.^(alphag).*agHat1.^(1-alphag);
msHat = lm.^(1-alpham).*xmHat1.^(alpham).*amHat1.^(1-alpham);

term1 = (pm/pe)^epsy*(deltaye/deltayn)^epsy; 
term2 = (deltaeg + (pg/pfTil)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1));
term3  = (pg/pfTil)^epse*(deltaef/deltaeg)^epse*term1;
term4= (deltaff + ((pf + tau0)/(pfStar +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));

pcons = pm + pg*term1/term2 + pf*term3/(term2*term4);

mdHat = 1/pcons*(pf*fsHat + pg*gsHat + pm*msHat); 
edHat = mdHat*(pm/pe)^epsy*(deltaye/deltayn)^epsy;

gdHat = edHat/term2;
fdHatTil = (pg/pfTil)^epse*(deltaef/deltaeg)^epse*gdHat;
fdHat = fdHatTil/term4; 
fdHatStar = ((pf + tau0)/(pfStar+tau0))^epsf*(deltafo/deltaff)^epsf*fdHat;

fdHat = fdHatTil/term4;

%excess supply of g and f
d2 = (gsHat-gdHat);  %gives pf and pg.  know pfStar0 so can find pm0.
d3 = (fsHat - fdHat);

%LABOR                      
wlfHat = nu*(1-alphaf).*pf.*lf.^(-alphaf).*xfHat1.^alphaf.*afHat1.^(1-alphaf);
wlgHat = (1-alphag).*pg.*lg.^(-alphag).*xgHat1.^alphag.*agHat1.^(1-alphag);
wlmHat = (1-alpham).*pm.*lm.^(-alpham).*xmHat1.^alpham.*amHat1.^(1-alpham);

%Workers' wages should be equated across sectors
d4 = (wlfHat - wlgHat); 
d5 = (wlfHat - wlmHat);


%Scientist market clearing
pfix = nu*alphaf.*pf.*(afHat1.*lf./xfHat1).^(1-alphaf);
pgix = alphag.*pg.*(agHat1.*lg./xgHat1).^(1-alphag);
pmix = alpham.*pm.*(amHat1.*lm./xmHat1).^(1-alpham);

%SCIENTIST MARKET CLEARING
wsfHat = eta*gamma*rhonf^(eta)*(1-alphaf)*xfHat1.*pfix./(sf1.^(1-eta).*afHat0.^(phi-1).*afHat1.*a1./a0);
wsgHat = eta*gamma*rhong^(eta)*(1-alphag)*xgHat1.*pgix./(sg1.^(1-eta).*agHat0.^(phi-1).*agHat1.*a1./a0);
wsmHat = eta*gamma*(1-alpham)*xmHat1.*pmix./(sm1.^(1-eta).*amHat0.^(phi-1).*amHat1.*a1./a0);

%Scientists' wages should be equated across sectors
d6 = wsfHat - wsgHat; 
d7 = wsfHat - wsmHat;

answ =[d2, d3, d4, d5, d6, d7];



