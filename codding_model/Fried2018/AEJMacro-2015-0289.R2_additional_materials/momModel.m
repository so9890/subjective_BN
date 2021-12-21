%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: momModel.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Reports model values of moments used for calibration and other
% equilibrium quantities.
% Inputs: 
% Parameters calibrated with a method of moments: pp
% Directly calibrated parameters: alphaf,alpham, deltaef, epse, epsf, epsy,
% L, phi ,rhonf, rhong, S
% Percent change in the foriegn oil price: deltaPfStarD
% Empirical ratio between foreign oil price and domestic fossil energy
% price on BGP: theta0
% BGP growth rate: n0
% Empirical values of the moments: momD
% Initial guesses: guess0, guessBGP
% fsolve options: options
% Number of time periods: Ttarg
% Maximum value of the simplex from minDist: yHigh
% Aggregate technology on the BGP: a0
% Outputs: 
% Sum of the squared distance between momD and momM: dist
% Fval value in BGP and shock period. Should be zero: solve (2X1)
% Initial guess for BGP: guessBGP
% model values of the moments: momM
% Scientist efficieny: gamma
% Values in the BGP and shock period. "Hat" indicates the de-trended
% quantities. 2X1 vectors. First entry is the BGP and second entry is the
% shock period.
% Note: "0" at the end of a variable indicates a BGP quantity while "1" at
% the end of a variable indicates a shock period quantity 
% Calls programs bgp.m and period1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, solve, guessBgp,  momM, gamma, a, af, ag, am,  agg, cHat, edHat, fdHat,...
    fdHatStar, fsHat, gdpHat, gdHat, gsHat, lf,lg, mdHat, msHat,pe, pf, pfStar, pfTil, pg,pm, ...
    pfix, pgix, pmix, profitsFHat, profitsFxHat, profitsGxHat, profitsMxHat, ...
    sf, sg, wlfHat, wlgHat, wlmHat, wsfHat, wsgHat, wsmHat, ...
    xfHat, xgHat, xmHat, yHat, shareCon, imports] = momModel(pp,a0,alphaf,alpham, deltaef, epse, epsf, epsy, L, phi ,rhonf, rhong,...
    S, theta0, deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  yHigh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatlize quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Carbon taxes are zero in the BGP and the shock period
tau = zeros(Ttarg, 1); 
tauStar = zeros(Ttarg,1);
tau0 = tau(1);
tauStar0 = tauStar(1);


tolInc = 1e-3;
tolIncBGP = 1e-7;

testing =0;
if testing == 1
    pp = [alphag, eta, deltaye, epsf, deltaff, nu1]; %don't include z1 or z10 pVec does not include p0 or p10
end

alphag = pp(1);  eta = pp(2); deltaye = pp(3); 
epsf = pp(4); deltaff = pp(5); nu1 = pp(6);


deltayn = 1-deltaye;
deltaeg=1-deltaef;
deltafo = 1-deltaff;

nuVec = [1, nu1]';
solve = zeros(Ttarg+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Pf and Pg along BGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repeat =0; error =0; %error turns on if something bad happens

options.TypicalX = guessBgp';

dum0 = @(vec)bgp(vec(1), vec(2), a0, alphaf, alphag,alpham, deltaef, deltaeg, deltaff, deltafo,deltaye, deltayn,...
    epse, epsf, epsy,eta, n0, phi, rhonf, rhong, theta0, L, S, tau0);


[val0, fval, exit] = fsolve(dum0,guessBgp, options);

try
    [val0, fval, exit] = fsolve(dum0,guessBgp, options);
    solve0 = max(abs(fval));
    
catch
    solve0 = 1;
    fval=1;
end

if solve0 > 1e-12 || isreal(fval) ==0
    repeat =1;
    guessBgpNew = guess0;
end

inc= 2;
while repeat ==1
    %%
    incVec = [0:1/inc:1];
    i =1;
    while i <=  length(incVec);
        %%
        incVec(i);
        pNew = incVec(i)*pp + (1-incVec(i))*p0;
        
        alphag1 = pNew(1);  eta1 = pNew(2); deltaye1 = pNew(3); 
        epsf1 = pNew(4); deltaff1 = pNew(5);   
        
        deltafo1 = 1-deltaff1;
        deltayn1 = 1-deltaye1;
        
        dum0 = @(vec) bgp(vec(1), vec(2), a0,alphaf,alphag1,alpham, deltaef, deltaeg, deltaff1, deltafo1, deltaye1, deltayn1,...
            epse, epsf1, epsy, eta1, n0, phi, rhonf, rhong, theta0, L, S, tau0);
        
        [val0, fval, exit] = fsolve(dum0,guessBgpNew, options);
        
        solve0 = max(abs(fval));
        
        if solve0 < 1e-12 && incVec(i) >1 && isreal(fval)==1;
            repeat =0;
            i = i+1;
        elseif solve0 < 1e-12 && isreal(fval) ==1
            i = i+1; % continue the loop
            guessBgpNew = val0;
        else
            inc = 2*inc; %go back one and try again
            incVec1 = [incVec(i-1):1/inc:1];
            incVec =[incVec(1:i-2), incVec1];
            i = i-1;
            if 1/inc< tolIncBGP
                disp('Steady state did not solve')
                repeat =0;
                error=1;
                i = length(incVec)+1;
                pp';
            end
        end
    end
end

solve(1) = solve0;

%Store solution as new initial guess
if error ==0
    guessBgp = val0;
end

pf0 = val0(1); 
pg0 = val0(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BGP quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PRICES
pfStar0 = theta0*pf0;
pfTil0 = ((pf0+ tau0)^(1-epsf)*deltaff^epsf + (pfStar0 + tau0)^(1-epsf)*deltafo^epsf)^(1/(1-epsf)); 
pe0 = (pfTil0^(1-epse)*deltaef^epse + pg0^(1-epse)*deltaeg^epse)^(1/(1-epse));
pm0 = ((1 - pe0^(1-epsy)*deltaye^epsy)/deltayn^epsy)^(1/(1-epsy));

pfix0 = 1/(alphaf);
pgix0 = 1/(alphag); 
pmix0 = 1/(alpham); 

%TECHNOLOGIES AND SCIENTISTS
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
    *pgix0^(alphag/(1-alphag))/pmix0^(alpham/(1-alpham)) *alpham^(1/(1-alpham))/alphag^(1/(1-alphag))*(1/rhong)^(eta);

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

%DEMANDS
term1 = (pm0/pe0)^epsy*(deltaye/deltayn)^epsy;
term2 = (deltaeg + (pg0/pfTil0)^(epse-1)*deltaeg^(1-epse)*deltaef^epse)^(epse/(epse-1));

term3  = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*term1;

term4= (deltaff + ((pf0 + tau0)/(pfStar0 +  tau0))^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf)^(epsf/(epsf-1));
pcons = pm0 + pg0*term1/term2 + pf0*term3/(term2*term4);

mdHat0 = 1/pcons*(pf0*fsHat0 + pg0*gsHat0 + pm0*msHat0);
edHat0 = mdHat0*(pm0/pe0)^epsy*(deltaye/deltayn)^epsy;
gdHat0 = edHat0/term2;
fdHatTil0 = (pg0/pfTil0)^epse*(deltaef/deltaeg)^epse*gdHat0;

fdHat0 = fdHatTil0/term4;
fdHatStar0 = ((pf0 + tau0)/(pfStar0+tau0))^epsf*(deltafo/deltaff)^epsf*fdHat0;

%WORKERS WAGES
wlfHat0 = (1-alphaf).*pf0.*lf0.^(-alphaf).*xfHat0.^alphaf.*afHat0.^(1-alphaf);
wlgHat0 = (1-alphag).*pg0.*lg0.^(-alphag).*xgHat0.^alphag.*agHat0.^(1-alphag);
wlmHat0 = (1-alpham).*pm0.*lm0.^(-alpham).*xmHat0.^alpham.*amHat0.^(1-alpham);

fdHatTil0 = (deltaff*fdHat0^((epsf-1)/epsf) + deltafo*fdHatStar0^((epsf-1)/epsf))^(epsf/(epsf-1));
edHat0 = (deltaeg*gdHat0^((epse-1)/epse) + deltaef*fdHatTil0^((epse-1)/epse))^(epse/(epse-1));
yHat0 = (deltaye*edHat0^((epsy-1)/epsy) + deltayn*mdHat0^((epsy-1)/epsy))^(epsy/(epsy-1));

%PROFITS
%machine producer profits should be positive 
profitsFx0 = pfix0.*xfHat0 - xfHat0 - wsfHat0*sf0;
profitsGx0 = pgix0*xgHat0 - xgHat0 - wsgHat0*sg0;
profitsMx0 = pmix0*xmHat0 - xmHat0 - wsmHat0*sm0;

%all other profits should be zero
profitsF0 = pf0*fsHat0 - wlfHat0*lf0 - pfix0*xfHat0 ;
profitsG0 = pg0*gsHat0 - wlgHat0*lg0 - pgix0*xgHat0;
profitsM0 = pm0*msHat0 - wlmHat0*lm0 - pmix0*xmHat0;
profitsE0 = pe0*edHat0 - pfTil0*fdHatTil0 - pg0*gdHat0;
profitsFtil0 = pfTil0*fdHatTil0 - pf0*fdHat0 - pfStar0*fdHatStar0;
profitsY0 = yHat0 - pe0*edHat0 - pm0*mdHat0;

%consumption
cHat0 = profitsFx0 + profitsGx0 + profitsMx0 + wlfHat0*L + wsfHat0*S;
%machines
xHat0 = xfHat0 + xgHat0 + xmHat0;
gdpHat0 =yHat0 - pfStar0*fdHatStar0 - xHat0;

%Check that agg resource constraint is zero
agg0 = gdpHat0  - cHat0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatilize variables for shock period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First period is BGP, second period is shock period

am = ones(Ttarg, 1);
ag = ones(Ttarg, 1);
af = ones(Ttarg, 1);
a = ones(Ttarg, 1);
sm = ones(Ttarg, 1);
sg = zeros(Ttarg, 1);
sf = ones(Ttarg, 1);

%Aggregates
pg = zeros(Ttarg,1);
pf = zeros(Ttarg,1);
lf = zeros(Ttarg,1);
lg = zeros(Ttarg,1);

xfHat = zeros(Ttarg,1);
xgHat=zeros(Ttarg,1);
xmHat = zeros(Ttarg,1);

%Initial Conditions
a(1) =a0;
af(1) = afHat0*a0;
ag(1) =agHat0*a0;
am(1) = amHat0*a0;

xfHat(1) = xfHat0;
xgHat(1)  = xgHat0;
xmHat(1) = xmHat0;
sf(1) = sf0;
sg(1) = sg0;
sm(1) = S - sf(1)- sg(1);
lf(1) = lf0;
lg(1) = lg0;
pf(1) =pf0;
pg(1) = pg0;

%CALCULATE PFSTAR FROM THE DATA
pfStar = (1+deltaPfStarD(1:Ttarg))*pfStar0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for pf, pg, lf, and lg in shock period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error==0
    repeat =0;
    guess1 = [lf0, lg0, pf0, pg0, sf0, sg0];
    options.TypicalX = guess1';
  
        dum1 = @(vec) period1(vec(1), vec(2), vec(3), vec(4), vec(5), vec(6), af0, ag0, am0,a0, pfStar(2), xfHat0, xgHat0, xmHat0,  ...
        alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn, epse, epsf, epsy, eta,...
        gamma, nu1, n0,phi,rhonf, rhong, S,tau0, L);
    
    
    [val1, fval1, exitflag1] = fsolve(dum1, guess1, options);
    solve1 = max(abs(fval1));

    if solve1 > 1e-12 || isreal(fval1) ==0
        repeat =1;
        guess1New = guess1;
    end
    
    inc = 2;
    incVec = [0:1/inc:1];
    i =1;
    while repeat ==1 && error ==0
        
        while i <=  length(incVec);
            incVec(i);
            pfStar2 = incVec(i)*pfStar(2) + (1-incVec(i))*pfStar(1);
            nu12 = incVec(i)*nu1 + (1-incVec(i))*1;
            
            dum1 = @(vec) period1(vec(1), vec(2), vec(3), vec(4), vec(5), vec(6), af0, ag0, am0,a0, pfStar2,...
                xfHat0, xgHat0, xmHat0,  ...
                alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn, epse, epsf, epsy, eta,...
                gamma, nu12, n0,phi,rhonf, rhong, S,tau0, L);
            
          
            [val1, fval1, exitflag1] = fsolve(dum1, guess1New, options);
            
            solve1 = max(abs(fval1));
            
            if solve1 < 1e-12 && incVec(i) ==1 && isreal(fval1) ==1
                repeat =0;
                i = i+1;
            elseif solve1 < 1e-12 && isreal(fval1) ==1
                i = i+1; % continue the loop
                guess1New = val1;
            else
                inc = 2*inc; %go back one and try again
                incVec1 = [incVec(i-1):1/inc:1];
                incVec =[incVec(1:i-2), incVec1];
                i = i-1;
                if 1/inc< tolInc
                    disp('disaster in period 1!')
                    repeat =0;
                    i = length(incVec)+1;
                    error=1;
                end
            end
            
        end
        
        
    end
   
    
    lf(2) = val1(1); lg(2) = val1(2); pf(2) = val1(3); pg(2) = val1(4); sf(2) = val1(5); sg(2) = val1(6); 
    sm(2) = S - sf(2) - sg(2);
         
    lm = L - lf - lg;
    
    %Update Technology;
    af(2) = af(1)*(1 + gamma*rhonf^(eta)*(sf(2))^eta*(a(1)/(af(1)))^phi);
    ag(2) = ag(1)*(1+ gamma*rhong^(eta)*(sg(2))^eta*(a(1)/(ag(1)))^phi);
    am(2) = am(1)*(1+ gamma*sm(2)^eta*(a(1)/am(1))^phi);
    a(2) = (1/rhonf*af(2) + 1/rhong*ag(2) + am(2))/(1 + 1/rhonf + 1/rhong);
    afHat = af./a; agHat = ag./a; amHat = am./a;   

    solve(2) = solve1;
    
    xfHat(2) = xfHat0*a0*(1+n0)/a(2); 
    xgHat(2) = xgHat0*a0*(1+n0)/a(2); 
    xmHat(2) = xmHat0*a0*(1+n0)/a(2);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Shock period (and BGP) quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %PRICES
    pfTil = ((pf+ tau(1:Ttarg)).^(1-epsf)*deltaff^epsf + (pfStar(1:Ttarg) + tauStar(1:Ttarg)).^(1-epsf)*deltafo^epsf).^(1/(1-epsf));
    pe = (pfTil.^(1-epse)*deltaef^epse + pg.^(1-epse)*deltaeg^epse).^(1/(1-epse));
    pm = ((1 - pe.^(1-epsy)*deltaye^epsy)/deltayn^epsy).^(1/(1-epsy));
    
    pfix = nuVec.*alphaf.*pf.*(afHat.*lf./xfHat).^(1-alphaf);
    pgix = alphag.*pg.*(agHat.*lg./xgHat).^(1-alphag);
    pmix = alpham.*pm.*(amHat.*lm./xmHat).^(1-alpham);
    
    %SCIENTISTS WAGES
    wsfHat = [wsfHat0; eta*gamma*rhonf^(eta)*(1-alphaf)*xfHat(2:end).*pfix(2:end)...
        ./(sf(2:end).^(1-eta).*afHat(1:end-1).^(phi-1).*afHat(2:end).*a(2:end)./a(1:end-1))];
    wsgHat = [wsgHat0; eta*gamma*rhong^(eta)*(1-alphag)*xgHat(2:end).*pgix(2:end)...
        ./(sg(2:end).^(1-eta).*agHat(1:end-1).^(phi-1).*agHat(2:end).*a(2:end)./a(1:end-1))];
    wsmHat = [wsmHat0; eta*gamma*(1-alpham)*xmHat(2:end).*pmix(2:end)...
        ./(sm(2:end).^(1-eta).*amHat(1:end-1).^(phi-1).*amHat(2:end).*a(2:end)./a(1:end-1))];
       
   %SUPPLIES
    fsHat = nuVec.*lf.^(1-alphaf).*xfHat.^(alphaf).*afHat.^(1-alphaf);
    gsHat = lg.^(1-alphag).*xgHat.^(alphag).*agHat.^(1-alphag);
    msHat = lm.^(1-alpham).*xmHat.^(alpham).*amHat.^(1-alpham);
    
    %DEMANDS
    term1 = (pm./pe).^epsy*(deltaye/deltayn)^epsy;
    term2 = (deltaeg + (pg./pfTil).^(epse-1)*deltaeg^(1-epse)*deltaef^epse).^(epse/(epse-1));
    term3  = (pg./pfTil).^epse*(deltaef/deltaeg).^epse.*term1;
    term4= (deltaff + ((pf + tau(1:Ttarg))./(pfStar(1:Ttarg) +  tauStar(1:Ttarg))).^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf).^(epsf/(epsf-1));
    
    pcons = pm + pg.*term1./term2 + pf.*term3./(term2.*term4);
    
    mdHat = 1./pcons.*(pf.*fsHat + pg.*gsHat + pm.*msHat);
    
    edHat = mdHat.*(pm./pe).^epsy*(deltaye/deltayn)^epsy;
    
    gdHat = edHat./term2;
    fdHatTil = (pg./pfTil).^epse*(deltaef/deltaeg).^epse.*gdHat;
    fdHat = fdHatTil./term4;
    fdHatStar = ((pf + tau(1:Ttarg))./(pfStar(1:Ttarg)+tauStar(1:Ttarg))).^epsf.*(deltafo/deltaff).^epsf.*fdHat;
    esHat = (deltaeg*gdHat.^((epse -1)/epse) + deltaef*fdHatTil.^((epse-1)/epse)).^(epse/(epse-1));

    yHat = (deltaye*esHat.^((epsy-1)/epsy) + deltayn*mdHat.^((epsy-1)/epsy)).^(epsy/(epsy-1));

    
    %WORKERS WAGES
    wlfHat = nuVec.*(1-alphaf).*pf.*lf.^(-alphaf).*xfHat.^alphaf.*afHat.^(1-alphaf);
    wlgHat = (1-alphag).*pg.*lg.^(-alphag).*xgHat.^alphag.*agHat.^(1-alphag);
    wlmHat = (1-alpham).*pm.*lm.^(-alpham).*xmHat.^alpham.*amHat.^(1-alpham);

    
    %PROFITS
    %Machine producer profits should be positive
    profitsFxHat = pfix.*xfHat - xfHat - wsfHat.*sf;
    profitsGxHat = pgix.*xgHat - xgHat - wsgHat.*sg;
    profitsMxHat = pmix.*xmHat - xmHat - wsmHat.*sm;
    
    %All other profits should be zero
    profitsFHat = pf.*fsHat - wlfHat.*lf - pfix.*xfHat;
    profitsGHat = pg.*gsHat - wlgHat.*lg - pgix.*xgHat;
    profitsMHat = pm.*msHat - wlmHat.*lm - pmix.*xmHat;
    profitsEHat = pe.*edHat - pfTil.*fdHatTil- pg.*gdHat;
    profitsFtilHat = pfTil.*fdHatTil - pf.*fdHat - pfStar(1:Ttarg).*fdHatStar;
    
    profitsYHat = yHat - pe.*edHat - pm.*mdHat;
    
    %CHECK AGGREGATE RESOURCE CONSTAINT
    %Consumption
    cHat = profitsFxHat + profitsGxHat + profitsMxHat + wlfHat*L + wsfHat.*sf + wsgHat.*sg + wsmHat.*sm...
        +tau(1:Ttarg).*fdHat + tauStar(1:Ttarg).*fdHatStar;
    
    xHat = xfHat + xgHat+ xmHat;
    gdpHat =yHat - pfStar(1:Ttarg).*fdHatStar - xHat;
    
    %Aggregate resource constraint should be zero
    agg = gdpHat -cHat;
    
    imports = pfStar(1:Ttarg).*fdHatStar./(pf.*fdHat + pfStar(1:Ttarg).*fdHatStar); %value of energy from abroad divided by total energy
    shareCon = (pf.*fdHat+ pfStar(1:Ttarg).*fdHatStar)./gdpHat;
    shareImports = pfStar(1:Ttarg).*fdHatStar./gdpHat;
    shareProd = pf.*fdHat./gdpHat;
    
    wsg = wsgHat.*a;
    wsm = wsmHat.*a;
    wsf = wsfHat.*a;
    
    fracSf = wsf.*sf./(wsf.*sf + wsg.*sg  + wsm.*sm);
    fracSg = wsg.*sg./(wsf.*sf + wsg.*sg  + wsm.*sm);
    
    momM  =[shareImports(1:Ttarg); shareProd(1:Ttarg); fracSf(2:Ttarg);...
        fracSg(2:Ttarg)]*100;
            
if error ==0
    dist = ((momM - momD))'*((momM - momD));
else dist = yHigh +1;
end



 


