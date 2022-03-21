%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: momModel3.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Computes equilibrium quantities in the endogenous and exogenous
% innovation models along the BGP and under two
% different policies (1) a carbon tax and (2) a subsidy to green research.
% Inputs:
% Parameters calibrated with a method of moments: pp
% Directly calibrated parameters: alphaf,alpham, deltaef, epse, epsf, epsy,
% L, phi ,rhonf, rhong, S
% Percent change in the foriegn oil price: deltaPfStarD
% Empirical ratio between foreign oil price and domestic fossil energy
% price on BGP: theta0
% BGP growth rate: n0
% Initial guesses: guess0, guessBGP
% fsolve options: options
% Number of time periods for the simulation: Tsim
% Green research subsidy (set to zero for paper results): v
% Aggregate technology on the BGP: a0

% Outputs:
% Indicator of whether the fsolve was successful: solve
% Initial guess for BGP: guessBGP
% model values of the moments: momM
% Scientist efficieny: gamma
% Equilibrium values from the model: TsimX1 vectors. First entry is the BGP.
% "Hat" indicates the de-trended
% quantities.
% Note: output dist and momM are not meaningful and are set to constants.
% Note: "0" at the end of a variable indicates a BGP quantity while "1" at
% the end of a variable indicates a shock period quantity
% Calls programs bgp.m, allinoneSubs.m and allinoneNoS.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, solve, guessBgp,  momM, gamma, a, af, ag, am,  agg, cHat, edHat, fdHat,...
    fdHatStar, fsHat, gdpHat, gdHat, gsHat, lf,lg, mdHat, msHat,pe, pf, pfStar, pfTil, pg,pm, ...
    pfix, pgix, pmix, profitsFHat, profitsFxHat, profitsGxHat, profitsMxHat, ...
    sf, sg, wlfHat, wlgHat, wlmHat, wsfHat, wsgHat, wsmHat, ...
    xfHat, xgHat, xmHat, yHat, shareCon, imports] = momModel3(pp,a0,alphaf,alpham, deltaef, epse, epsf, epsy, L, noInnov,phi ,rhonf, rhong,...
    S, tau, tauStar, theta0, n0, guess0, guessBgp, options,p0,Tsim, v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatlize quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Carbon taxes are zero in the BGP
tau0 = 0;

tolInc = 1e-3;
tolIncBGP = 1e-7;

alphag = pp(1);  eta = pp(2); deltaye = pp(3);
epsf = pp(4); deltaff = pp(5);

deltayn = 1-deltaye;
deltaeg = 1-deltaef;
deltafo = 1-deltaff;

solve = zeros(Tsim,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Pf and Pg along BGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


repeat =0; error =0; %turns on if something bad happens

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


%% the following looks like the calibration routine
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
        
        
        deltaeg=1-deltaef;
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
                disp('BGP did not solve')
                repeat =0;
                error=1;
                i = length(incVec)+1;
                pp';
            end
        end
    end
end

%Store solution as new initial guess
solve(1) = solve0;
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

% DEMANDS
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
cHat0 = profitsFx0 + profitsGx0 + profitsMx0 + profitsF0 + wlfHat0*L + wsfHat0*S;
%machines
xHat0 = xfHat0 + xgHat0 + xmHat0;
gdpHat0 =yHat0 - pfStar0*fdHatStar0 - xHat0;

%Check that the agg resource constraint is zero.
agg0 = gdpHat0  - cHat0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatilize variables for the policy periods (2-Tsim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if noInnov==1
    Tsim =2; %Only simulate policy for one period in exog. model.
    %(because go immediately to new long-run eq.)
elseif max(abs(tau))==0 && max(abs(v)) ==0 %Don't simulate policy in
    %in baseline (no carbon tax and no subsidy)
    Tsim =1;
end

am = ones(Tsim, 1);
ag = ones(Tsim, 1);
af = ones(Tsim, 1);
a = ones(Tsim, 1);
sm = ones(Tsim, 1);
sg = zeros(Tsim, 1);
sf = ones(Tsim, 1);

%Aggregates
pg = zeros(Tsim,1);
pf = zeros(Tsim,1);
lf = zeros(Tsim,1);
lg = zeros(Tsim,1);

xfHat = zeros(Tsim,1);
xgHat=zeros(Tsim,1);
xmHat = zeros(Tsim,1);

%Period 1 is the BGP
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

%Tech grows at BGP rate in exogenous model
if noInnov ==1
    af(2) = (1+n0)*af(1);
    ag(2) = (1+n0)*ag(1);
    am(2) = (1+n0)*am(1);
    a(2) = (1+n0)*a(1);
end

pfStar = pfStar0*ones(Tsim,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for pf, pg, lf, lg, sf, sg, xfHat, xgHat, xmHat under policy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Endogenous model
t =2;
error =0;
if (max(abs(tau))>0 || max(abs(v)) >0) && noInnov ==0
    guess2 = [lf0, lg0, pf0, pg0,sf0, sg0, xfHat0, xgHat0, xmHat0];
    guess = guess2;
    while t <= Tsim
        
        options.TypicalX = guess2';
        
        repeat =0;
        
        dum = @(vec) allinoneSubs(vec(1), vec(2),vec(3), vec(4), vec(5), vec(6), vec(7), vec(8), vec(9), ....
            af(t-1), ag(t-1), am(t-1), a(t-1), alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn,...
            epse, epsf, epsy, eta, gamma,pfStar0, phi,rhonf, rhong, tau(t),tauStar(t), L, S, v(t));
        
        [val, fval, exitflag] = fsolve(dum, guess, options);
        solve(t) = max(abs(fval));
        
        if solve(t) > 1e-12 || isreal(fval) ==0
            repeat =1;
            sv = [af(t-1), ag(t-1), am(t-1), a(t-1), tau(t), tauStar(t)];
        end
        
        guessNew = guess2; %Should be the correct solution for the BGP
        inc = 2;
        incVec = [0:1/inc:1];
        i =1;
        while repeat ==1 && error ==0
            while i <=  length(incVec);
                incVec(i);
                
                sv1 = incVec(i)*sv+ (1-incVec(i))*[af0, ag0, am0, a0, 0,0];
                
                dum = @(vec) allinoneSubs(vec(1), vec(2),vec(3), vec(4), vec(5), vec(6), vec(7), vec(8), vec(9), ...
                    sv1(1), sv1(2), sv1(3),sv1(4), alphaf, alphag, alpham, ...
                    deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn,...
                    epse,epsf, epsy, eta, gamma, pfStar0, phi, rhonf, rhong,  ...
                    sv1(5),sv1(6), L, S, v(t));
                
                [val, fval, exitflag] = fsolve(dum, guessNew, options);
                
                solve(t) = max(abs(fval));
                solve2 = solve(t);
                
                if solve2 < 1e-12 && incVec(i) == 1 && isreal(fval)==1;
                    repeat =0;
                    i = i+1;
                elseif solve2 < 1e-12 && isreal(fval) ==1;
                    i = i+1; % continue the loop
                    guessNew = val;
                else
                    inc = 2*inc; %go back one and try again with smaller increments
                    if i >1
                        incVec1 = [incVec(i-1):1/inc:1];
                        incVec =[incVec(1:i-2), incVec1];
                        i = i-1;
                    else incVec =[0:1/inc:1];
                    end
                    if 1/inc< tolInc
                        disp('Policy did not solve in endog-model')
                        exitflag
                        repeat =0;
                        error=1;
                        i = length(incVec)+1;
                    end
                end
            end
        end
        
        lf(t)=  val(1);
        lg(t) = val(2);
        pf(t) = val(3);
        pg(t) = val(4);
        sf(t) = val(5);
        sg(t)=  val(6);
        xfHat(t) = val(7);
        xgHat(t) = val(8);
        xmHat(t)= val(9);
        
        sm(t) = S - sf(t) - sg(t);
        af(t) = af(t-1)*(1 + gamma*rhonf^(eta)*sf(t)^eta*(a(t-1)/(af(t-1)))^phi);
        ag(t) = ag(t-1)*(1+ gamma*rhong^(eta)*(sg(t))^eta*(a(t-1)/(ag(t-1)))^phi);
        am(t) = am(t-1)*(1+ gamma*sm(t)^eta*(a(t-1)/am(t-1))^phi);
        a(t) = (1/rhonf*af(t) + 1/rhong*ag(t) + am(t))/(1 + 1/rhonf + 1/rhong);
        
        t = t+1;
    end
end


%Exogneous model
if  noInnov ==1 && max(abs(tau))> 0
    %%
    guess3 = [lf0, lg0, pf0, pg0,xfHat0, xgHat0, xmHat0];
    
    af(t) = (1+n0)*af(t-1);
    ag(t) = (1+n0)*ag(t-1);
    am(t) = (1+n0)*am(t-1);
    a(t) = (1+n0)*a(t-1);
    
    options.TypicalX = guess3';
    
    dum = @(vec) allinoneNoS(vec(1), vec(2),vec(3), vec(4), vec(5), vec(6), vec(7),af(t), ag(t), am(t), a(t), ...
        alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn,...
        epse, epsf, epsy, pfStar0, tau(t), tauStar(t), L);
    
    [val, fval, exitflag] = fsolve(dum, guess3, options);
    solve(t) = max(abs(fval));
    
    if solve(t) > 1e-12 || isreal(fval) ==0
        repeat =1;
        sv = [af(t), ag(t), am(t), a(t), pfStar0,tau(t+1), tauStar(t+1)];
    end
    
    guessNew = guess3; %Should be the correct solution for the BGP
    inc = 2;
    incVec = [0:1/inc:1];
    i =1;
    while repeat ==1 && error ==0
        while i <=  length(incVec);
            incVec(i);
            
            sv1 = incVec(i)*sv+ (1-incVec(i))*[af0, ag0, am0, a0, pfStar0, 0, 0];
            
            dum = @(vec) allinoneNoS(vec(1), vec(2),vec(3), vec(4), vec(5), vec(6), vec(7),sv1(1),sv1(2),sv1(3),sv1(4), ...
                alphaf,alphag, alpham, deltaef, deltaeg, deltaff, deltafo, deltaye, deltayn,...
                epse, epsf, epsy, sv1(5), sv1(6), sv1(7), L);
            
            [val, fval, exitflag] = fsolve(dum, guessNew, options);
            
            solve(t) = max(abs(fval));
            solve2 = solve(t);
            
            if solve2 < 1e-12 && incVec(i) == 1 && isreal(fval)==1;
                repeat =0;
                i = i+1;
            elseif solve2 < 1e-12 && isreal(fval) ==1;
                i = i+1; % continue the loop
                guessNew = val;
            else
                inc = 2*inc; %go back one and try again with smaller increments
                if i >1
                    incVec1 = [incVec(i-1):1/inc:1];
                    incVec =[incVec(1:i-2), incVec1];
                    i = i-1;
                else incVec =[0:1/inc:1];
                end
                if 1/inc< tolInc
                    disp('Policy did not solve in exog-model')
                    repeat =0;
                    error=1;
                    i = length(incVec)+1;
                end
            end
        end
    end
    
    lf(t)= val(1);
    lg(t) = val(2);
    pf(t) =  val(3);
    pg(t) =  val(4);
    xfHat(t) = val(5);
    xgHat(t) = val(6);
    xmHat(t)= val(7);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate equilibrium quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

afHat = af./a; agHat = ag./a; amHat = am./a;
nuVec = ones(Tsim,1);
lm = L - lf - lg;

%PRICES
pfTil = ((pf+ tau(1:Tsim)).^(1-epsf)*deltaff^epsf + (pfStar(1:Tsim) + tauStar(1:Tsim)).^(1-epsf)*deltafo^epsf).^(1/(1-epsf));
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
term4= (deltaff + ((pf + tau(1:Tsim))./(pfStar(1:Tsim) +  tauStar(1:Tsim))).^(epsf-1)*deltaff^(1-epsf)*deltafo^epsf).^(epsf/(epsf-1));

pcons = pm + pg.*term1./term2 + pf.*term3./(term2.*term4);

mdHat = 1./pcons.*(pf.*fsHat + pg.*gsHat + pm.*msHat);

edHat = mdHat.*(pm./pe).^epsy*(deltaye/deltayn)^epsy;

gdHat = edHat./term2;
fdHatTil = (pg./pfTil).^epse*(deltaef/deltaeg).^epse.*gdHat;
fdHat = fdHatTil./term4;
fdHatStar = ((pf + tau(1:Tsim))./(pfStar(1:Tsim)+tauStar(1:Tsim))).^epsf.*(deltafo/deltaff).^epsf.*fdHat;

%WORKERS WAGES
wlfHat = nuVec.*(1-alphaf).*pf.*lf.^(-alphaf).*xfHat.^alphaf.*afHat.^(1-alphaf);
wlgHat = (1-alphag).*pg.*lg.^(-alphag).*xgHat.^alphag.*agHat.^(1-alphag);
wlmHat = (1-alpham).*pm.*lm.^(-alpham).*xmHat.^alpham.*amHat.^(1-alpham);

esHat = (deltaeg*gdHat.^((epse -1)/epse) + deltaef*fdHatTil.^((epse-1)/epse)).^(epse/(epse-1));

yHat = (deltaye*esHat.^((epsy-1)/epsy) + deltayn*mdHat.^((epsy-1)/epsy)).^(epsy/(epsy-1));

%PROFITS
%Machine prpducer profits should be positive
profitsFxHat = pfix.*xfHat - xfHat - wsfHat.*sf;
profitsGxHat = pgix.*xgHat - xgHat - wsgHat.*sg;
profitsMxHat = pmix.*xmHat - xmHat - wsmHat.*sm;

%All other profits should be zero
profitsFHat = pf.*fsHat - wlfHat.*lf - pfix.*xfHat;
profitsGHat = pg.*gsHat - wlgHat.*lg - pgix.*xgHat;
profitsMHat = pm.*msHat - wlmHat.*lm - pmix.*xmHat;
profitsEHat = pe.*edHat - pfTil.*fdHatTil- pg.*gdHat;
profitsFtilHat = pfTil.*fdHatTil - pf.*fdHat - pfStar(1:Tsim).*fdHatStar;

profitsYHat = yHat - pe.*edHat - pm.*mdHat;

%CHECK TOE AGG RESOURCE CONSTRAINT IN THE ENDOG-MODEL
cHat1 = profitsFxHat + profitsGxHat + profitsMxHat + profitsFHat + wlfHat*L + wsfHat.*sf + wsgHat.*sg + wsmHat.*sm...
    +tau(1:Tsim).*fdHat + tauStar(1:Tsim).*fdHatStar; %should be NAN in exog-model

xHat = xfHat + xgHat+ xmHat;
gdpHat =yHat - pfStar(1:Tsim).*fdHatStar - xHat;
agg = gdpHat -cHat1;

cHat = gdpHat;

imports = pfStar(1:Tsim).*fdHatStar./(pf.*fdHat + pfStar(1:Tsim).*fdHatStar); %value of energy from abroad divided by total energy
shareCon = (pf.*fdHat+ pfStar(1:Tsim).*fdHatStar)./gdpHat;


momM =1;
dist =0;


