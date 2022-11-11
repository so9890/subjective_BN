%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: allinoneSubs.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes FOC and market clearing in the endog-model.
% Inputs: 
% Guesses for pf, pg, lf, lg, xfHat, xgHat, xmHat, sf and sg: pf1, pg1, 
%    lf1, lg1,xfHat1, xmHat1, xgHat1, sf1, sg1
% Model parameters and values: alphaf,alphag, alpham, deltaef, deltaeg, deltaff, 
%     deltafo, deltaye, deltayn, epse, epsf, epsy, pfStar1, tau, tauStar, L, v 
% Outputs: 
% First order and market clearing conditions (should all equal zero in eq.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function answ = allinoneSubs(lf1,lg1,pf1,pg1,sf1,sg1,...
    xfHat1, xgHat1, xmHat1, af0, ag0, am0,a0, ...
    alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn,...
    epse, epsf, epsy, eta, gamma, pfStar1, phi,rhonf, rhong, tau, tauStar, L, S, v)

%TECHNOLOGY
sm = S - sf1 - sg1;
af1 = af0*(1+ gamma*rhonf^(eta)*(sf1).^eta.*(a0./(af0)).^phi);
ag1 = ag0*(1+ gamma*rhong^(eta)*(sg1).^eta.*(a0./(ag0)).^phi); 
am1 = am0*(1+ gamma*sm.^eta.*(a0./am0).^phi); 
a1  = (1/rhonf*af1 + 1/rhong*ag1 + am1)/(1/rhonf + 1/rhong + 1);

afHat1 = af1./a1;
agHat1 = ag1./a1;
amHat1 = am1./a1;

afHat0 = af0./a0;
agHat0 = ag0./a0;
amHat0 = am0./a0;

%PRICES
pfTil1 = ((pf1 + tau).^(1-epsf)*deltaff^epsf + (pfStar1  + tauStar).^(1-epsf)*deltafo^epsf).^(1/(1-epsf)); 
pe1 = (pfTil1.^(1-epse)*deltaef^epse + pg1.^(1-epse)*deltaeg^epse).^(1/(1-epse));
pm1 = ((1 - pe1.^(1-epsy)*deltaye^epsy)/deltayn^epsy).^(1/(1-epsy));

%LABOR IN SECTOR M
lm1 = L - lf1 - lg1;
 
%GOODS MARKET CLEARING
fsHat1 = lf1.^(1-alphaf).*xfHat1.^(alphaf).*afHat1.^(1-alphaf);
gsHat1 = lg1.^(1-alphag).*xgHat1.^(alphag).*agHat1.^(1-alphag);
msHat1 = lm1.^(1-alpham).*xmHat1.^(alpham).*amHat1.^(1-alpham);

term1 = (pm1./pe1).^epsy*(deltaye/deltayn)^epsy; 
term2 = (deltaeg + (pg1./pfTil1).^(epse-1)*deltaeg^(1-epse)*deltaef^epse).^(epse/(epse-1));
term3  = (pg1./pfTil1).^epse*(deltaef/deltaeg).^epse.*term1;
term4= (deltaff + ((pf1 + tau)./(pfStar1 +  tauStar)).^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf).^(epsf/(epsf-1));

pcons = pm1 + pg1.*term1./term2 + pf1.*term3./(term2.*term4);
mdHat1 = 1./pcons.*(pf1.*fsHat1 + pg1.*gsHat1 + pm1.*msHat1);
edHat1 = mdHat1.*(pm1./pe1).^epsy*(deltaye/deltayn)^epsy;
gdHat1 = edHat1./term2; 
fdHatTil1 = (pg1./pfTil1).^epse*(deltaef/deltaeg).^epse.*gdHat1;
fdHat1 = fdHatTil1./term4;

%Excess supplies of f and g
d1 = (gsHat1 - gdHat1);
d2 = (fsHat1 - fdHat1); 

%LABOR MARKET CLEARING
wlfHat1 = (1-alphaf).*pf1.*lf1.^(-alphaf).*xfHat1.^alphaf.*afHat1.^(1-alphaf);
wlgHat1 = (1-alphag).*pg1.*lg1.^(-alphag).*xgHat1.^alphag.*agHat1.^(1-alphag);
wlmHat1 = (1-alpham).*pm1.*lm1.^(-alpham).*xmHat1.^alpham.*amHat1.^(1-alpham);

%Workers' wages should be equated across sectors
d3= (wlfHat1 - wlgHat1); 
d4 = ((wlfHat1 - wlmHat1)); 

%FOC FOR THE NUMBER OF MACHINES
pfix1 = alphaf.*pf1.*(afHat1.*lf1./xfHat1).^(1-alphaf);
pgix1 = alphag.*pg1.*(agHat1.*lg1./xgHat1).^(1-alphag);
pmix1 = alpham.*pm1.*(amHat1.*lm1./xmHat1).^(1-alpham);

d5 = alphaf*pfix1 -1;
d6 = alphag*pgix1 -1;
d7 = alpham*pmix1 -1;

%SCIENTIST MARKET CLEARING
wsfHat1 = eta*gamma*rhonf^(eta)*(1-alphaf)*xfHat1.*pfix1./(sf1.^(1-eta).*afHat0.^(phi-1).*afHat1.*a1./a0);
wsgHat1 = eta*gamma*rhong^(eta)*(1-alphag)*xgHat1.*pgix1./(sg1.^(1-eta).*agHat0.^(phi-1).*agHat1.*a1./a0)*(1/(1-v)); 
wsmHat1 = eta*gamma*(1-alpham)*xmHat1.*pmix1./(sm.^(1-eta).*amHat0.^(phi-1).*amHat1.*a1./a0); 

%Scientists' wages should be equated across sectors
d8 = (wsfHat1 - wsmHat1); 
d9 = (wsgHat1 - wsfHat1); 

answ =[d1, d2, d3, d4, d5,d6, d7, d8, d9];




