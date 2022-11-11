%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: simulation.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Simulate the economy under the following four scenarios:
% (1) Baseline with no carbon tax; (2) Carbon tax with endogenous
% innovation; (3) carbon tax with exogenous innovation; (4) Green research
% subsidy
% Note: "hat" indicatees a detrended quantity
% Calls programs momModel3.m and cev.m
% Saves results in baseResults.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear;
home =0;

if home ==1
    bpath = '/Users/steph/Dropbox/Research/Price_shocks/Matlab_files';
    f_path = '/Users/steph/Dropbox/Research/Price_shocks/Presentations/ASU2017';
else
    bpath = '/Users/sfried/Dropbox/Research/Price_shocks/Matlab_files';
    f_path = '/Users/sfried/Dropbox/Research/Price_shocks/Presentations/ASU2017';
    
end

cd(strcat(bpath, '/Uncertainty/AEJFiles'))
load(strcat(bpath, '/Uncertainty/AEJFiles/','base'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iniatlize quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(strcat(bpath, '/Uncertainty/AEJFiles/','tax')); %Carbon tax necessary to achieve target
tauVec = [0, tauI, tauNI]; %tax that achieves a 30% reduction in emissions; keep tax the same to calculate elasticity for different cases.
Tsim = 100; %number of simulated periods 
T = Tsim;
nExp = 3; %number pf experiments

pp = [alphag, eta, deltaye, epsf, deltaff, nu1];
theta1 = 2.1579; % ratio of pfStar/pf (equals average value f from 2001-2010 )
tauMat = zeros(T, nExp);
tauStarMat = zeros(T, nExp);
vMat = zeros(T, nExp);  %green research subsidy (set to zero)

tauStarVec = 1.045328173*tauVec; %conversion factor is metric tons co2 per mbtu FdStar/metric tons co2 per MBTU Fd

for i = 1:nExp
    tauMat(2:end, i) = tauVec(i)*ones(Tsim-1, 1);
    tauStarMat(2:end, i) = tauStarVec(i)*ones(Tsim-1,1);
end

a = zeros(Tsim, nExp);
af = zeros(Tsim, nExp);
ag = zeros(Tsim, nExp);
am = zeros(Tsim, nExp);
agg = zeros(Tsim, nExp);
cHat = zeros(Tsim, nExp);
edHat = zeros(Tsim, nExp);
fdHat = zeros(Tsim, nExp);
fdHatStar = zeros(Tsim, nExp);
fsHat = zeros(Tsim, nExp);
gdpHat= zeros(Tsim, nExp);
gdHat = zeros(Tsim, nExp);
gsHat = zeros(Tsim, nExp);
lf = zeros(Tsim, nExp);
lg = zeros(Tsim, nExp);
mdHat = zeros(Tsim, nExp);
msHat = zeros(Tsim, nExp);
pf = zeros(Tsim, nExp);
pfStar = zeros(Tsim, nExp);
pfTil = zeros(Tsim, nExp);
pg = zeros(Tsim, nExp);
pm= zeros(Tsim, nExp);
sf = zeros(Tsim, nExp);
sg = zeros(Tsim, nExp);
wlfHat = zeros(Tsim, nExp);
wlgHat= zeros(Tsim, nExp);
wlmHat = zeros(Tsim, nExp);
wsfHat = zeros(Tsim, nExp);
wsgHat = zeros(Tsim, nExp);
wsmHat = zeros(Tsim, nExp);
xfHat= zeros(Tsim, nExp);
xgHat= zeros(Tsim, nExp);
xmHat= zeros(Tsim, nExp);
yHat= zeros(Tsim, nExp);
pe= zeros(Tsim, nExp);
profitsFxHat= zeros(Tsim, nExp);
profitsFHat= zeros(Tsim, nExp);
profitsGxHat= zeros(Tsim, nExp);
profitsMxHat= zeros(Tsim, nExp);
pfix= zeros(Tsim, nExp);
pgix= zeros(Tsim, nExp);
pmix= zeros(Tsim, nExp);
share= zeros(Tsim, nExp);
imports= zeros(Tsim, nExp);
solve= zeros(Tsim, nExp);
noInnov = 0;
c1 = 1; c2 = Tsim;

for i = 1:nExp %[baseline, Innovation, No-innovation]
    i
    if i ==3
        noInnov=1;
        c1 = 1;
        c2 =2;
    end
    
    [dist, solve(:, i), guessBgp1,  momM1, gamma, a(c1:c2, i), af(c1:c2, i), ag(c1:c2, i), am(c1:c2, i),  agg(c1:c2, i),...
        cHat(c1:c2, i), edHat(c1:c2, i), fdHat(c1:c2, i),...
        fdHatStar(c1:c2, i), fsHat(c1:c2, i), gdpHat(c1:c2, i), gdHat(c1:c2, i), gsHat(c1:c2, i), lf(c1:c2, i),...
        lg(c1:c2, i), mdHat(c1:c2, i), msHat(c1:c2, i),pe(c1:c2, i), pf(c1:c2, i), pfStar(c1:c2, i), pfTil(c1:c2, i), pg(c1:c2, i),...
        pm(c1:c2, i), pfix(c1:c2, i), pgix(c1:c2, i), pmix(c1:c2, i),...
        profitsFHat(c1:c2, i), profitsFxHat(c1:c2, i), profitsGxHat(c1:c2, i), profitsMxHat(c1:c2, i),...
        sf(c1:c2, i), sg(c1:c2, i), wlfHat(c1:c2, i), wlgHat(c1:c2, i), wlmHat(c1:c2, i), wsfHat(c1:c2, i), wsgHat(c1:c2, i),...
        wsmHat(c1:c2, i), xfHat(c1:c2, i), xgHat(c1:c2, i), xmHat(c1:c2, i), yHat(c1:c2, i), share(c1:c2, i), imports(c1:c2, i)] ...
        = momModel3(pp,a0,alphaf,alpham, deltaef, epse, epsf, epsy, L,noInnov,phi ,rhonf, rhong,...
        S,tauMat(:,i), tauStarMat(:,i), theta1, n0, guess0, guessBgp, options,p0,Tsim, vMat(:,i));
    
end

sm = S - sf - sg;
lm = L- lf-lg;


%ADJUST TECHNOLOGIES (momModel3 only computes BGP values once)
for i = 1:Tsim -1
    a(i+1,1) = (1+n0)*a(i,1);
    af(i+1,1) = (1+n0)*af(i,1);
    ag(i+1,1) = (1+n0)*ag(i,1);
    am(i+1,1) = (1+n0)*am(i,1);
end


%ADJUST NOINNOVATION CASE (momModel3 only computes no-innovation case once)
if nExp > 2
    c3 = c2+1;
    c4 = Tsim;
    cHat(c3:c4,3) = cHat(2, 3)*ones(size(cHat(c3:c4))); edHat(c3:c4,3) = edHat(2, 3)*ones(size(edHat(c3:c4)));
    fdHat(c3:c4,3) = fdHat(2, 3)*ones(size(fdHat(c3:c4))); fdHatStar(c3:c4,3) = fdHatStar(2, 3)*ones(size(cHat(c3:c4)));
    fsHat(c3:c4,3) = fsHat(2, 3)*ones(size(cHat(c3:c4))); gsHat(c3:c4,3) = gsHat(2,3)*ones(size(cHat(c3:c4)));
    gdpHat(c3:c4,3) = gdpHat(2,3)*ones(size(cHat(c3:c4))); msHat(c3:c4,3) = msHat(2,3)*ones(size(cHat(c3:c4)));
    mdHat(c3:c4,3) = mdHat(2,3)*ones(size(cHat(c3:c4))); gdHat(c3:c4,3) = gdHat(2,3)*ones(size(cHat(c3:c4)));
    pm(c3:c4,3) = pm(2,3)*ones(size(cHat(c3:c4))); pg(c3:c4,3) = pg(2,3)*ones(size(cHat(c3:c4)));
    pf(c3:c4,3) = pf(2,3)*ones(size(cHat(c3:c4))); xfHat(c3:c4,3) = xfHat(2,3)*ones(size(cHat(c3:c4)));
    xgHat(c3:c4,3) = xgHat(2,3)*ones(size(cHat(c3:c4))); xmHat(c3:c4,3) = xmHat(2,3)*ones(size(cHat(c3:c4)));
    wlfHat(c3:c4,3) = wlfHat(2,3)*ones(size(cHat(c3:c4))); wlgHat(c3:c4,3) = wlgHat(2,3)*ones(size(cHat(c3:c4)));
    wlmHat(c3:c4,3) = wlmHat(2,3)*ones(size(cHat(c3:c4))); yHat(c3:c4,3) = yHat(2,3)*ones(size(cHat(c3:c4)));
    share(c3:c4,3) = share(2,3)*ones(size(cHat(c3:c4))); agg(c3:c4,3) = agg(2,3)*ones(size(cHat(c3:c4)));
    pe(c3:c4,3) = pe(2,3)*ones(size(cHat(c3:c4)));
    
    a(c2:c4, 3)= a(c2:c4, 1); af(c2:c4, 3)= af(c2:c4, 1); ag(c2:c4, 3)= ag(c2:c4, 1); am(c2:c4, 3)= am(c2:c4, 1);
    sf(c2:c4, 3)= sf(c2:c4, 1); sg(c2:c4, 3)= sg(c2:c4, 1); sm(c2:c4, 3)= sm(c2:c4, 1);
end


%Calculate aggregate values
c = cHat.*a; ed = edHat.*a; fd = fdHat.*a; fdStar = fdHatStar.*a; fs = fsHat.*a; gdp  = gdpHat.*a; gd = gdHat.*a;
gs = gsHat.*a; md = mdHat.*a; ms = msHat.*a; wlf = wlfHat.*a; wlg = wlgHat.*a; wlm = wlmHat.*a; wsf = wsfHat.*a;
wsg = wsgHat.*a; wsm = wsmHat.*a; xf = xfHat.*a; xg = xgHat.*a; xm = xmHat.*a;  y= yHat.*a;

% Calculate emissions
%units:
ind2015=2;%first period of the tax
ind2030 = 5; %20-year mark (start of the 2030-2034 time period)

d1 = 282.030793*1e9/fd(ind2015);  %MBTU fd 2010/ fd2010 maps fd in model to fd in data
d3 = 0.072556484/(3.667*1e9);  %giga  tons Carbon per MBTU fd emissions
d4 =0.07454/(3.667*1e9); %giga tons carbon per MBTU fdStar emissions fdStar per btu
d5 = 6.187196613;%emissions world per emissions US in 2010

cfd = d1*d3*d5;
cfdStar = d1*d4*d5;
emissions = cfd*fd + cfdStar*fdStar;

% PERCENT CHANGES
%values of variables in 2030-2034 time period
var2030 = [sf(ind2030, :); sg(ind2030, :); sm(ind2030, :);...
     ag(ind2030,:)./af(ind2030, :); ag(ind2030,:)./am(ind2030, :);...
   pg(ind2030,:)./(pf(ind2030,:));  pg(ind2030,:)./(pm(ind2030,:)); pe(ind2030,:)./(pm(ind2030,:));
     gd(ind2030, :)./fd(ind2030, :); ed(ind2030, :)./md(ind2030, :); wsf(ind2030, :)./wlf(ind2030, :);emissions(ind2030, :)];

 indlr = Tsim;

 
 varlr = [sf(indlr, :); sg(indlr, :); sm(indlr, :);...
     ag(indlr,:)./af(indlr, :); ag(indlr,:)./am(indlr, :);...
   pg(indlr,:)./(pf(indlr,:));  pg(indlr,:)./(pm(indlr,:)); pe(indlr,:)./(pm(indlr,:));
     gd(indlr, :)./fd(indlr, :); ed(indlr, :)./md(indlr, :); wsf(indlr, :)./wlf(indlr, :); emissions(indlr, :)];
 
%Percent change from baseline and tax in 2030-2034 time period with innovataion
change2 =(var2030(:, 2) - var2030(:, 1))./var2030(:, 1)*100;
%Percent change from baseline and tax in 2030-2034 time period without innovataion
if nExp>2
    change3= (var2030(:, 3) - var2030(:, 1))./var2030(:, 1)*100;
else change3 =0;
end
 
 
%Percent change from baseline and tax in long run with innovataion
changelr =  (varlr(:, 2) - varlr(:, 1))./varlr(:, 1)*100;
 
% Calculate CEV
theta =2; %CRRA parameter
discount = .985; %rate of time preference = 1.5%
discountHat = discount*(1.01)^(1-theta); %detrend discount by population
discountHat5 = discountHat^5; %adjust for the 5 year time period
guess = 1;
cevNoDamage = zeros(nExp,1);

cBase2 = c(ind2015+1:indlr, 1);
for i = 2:nExp
    cTax2 = c(ind2015+1:indlr, i);
    dum = @(lambda) cev(lambda, cBase2(1:end), cTax2(1:end), discountHat5, theta);
    [cevNoDamages, fval, exit] = fzero(dum, guess);
    cevNoDamage(i) = (cevNoDamages-1)*100;
end

cevNoClimateValue = cevNoDamage/100*53671; %Convert to (2012) consumption

% ELASTICITY OF GREEN IDEAS WRT COMPOSITE FOSSIL ENERGY PRICE
deltaSg = (sg(:,2) - sg(:,1))./sg(:,1);
greenIdeas = gamma*(sg.*rhong).^eta.*a.^phi.*ag.^(1-phi);
deltaIdeas = (greenIdeas(:,2) - greenIdeas(:,1))./greenIdeas(:,1);
deltaPfTil = (pfTil(:,2) - pfTil(:,1))./pfTil(:,1);
pfT = pf + tauMat;
deltaPfT = (pfT(:,2) - pfT(:,1))./pfT(:,1);
elasIdeas = deltaIdeas./deltaPfTil;
elasIdeas(2);

deltaEmissionsI = (emissions(:,2) - emissions(:,1))./emissions(:,1)*100;
deltaEmissionsNI = (emissions(:,3) - emissions(:,1))./emissions(:,1)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save baseResults

%% Sonja plots

% figure(1)
% plot(1:Tsim, pm(:,1))
