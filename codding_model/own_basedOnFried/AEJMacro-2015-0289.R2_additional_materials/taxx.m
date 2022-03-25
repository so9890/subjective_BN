%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: tax.m
% By: Stephie Fried
% Date: Spring 2017
% Purpose: Calculate the carbon tax necessary to achieve a 30 percent
% reduction in emissions in 20 years, in the exogenous and endogenous
% innovation models. Uses a bracketing algorithm.
% Calls program momModel3.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
% home =0;
% 
% if home ==1
%     bpath = '/Users/steph/Dropbox/Research/Price_shocks/Matlab_files';
% else
     bpath = '/home/sonja/Documents/projects/Overconsumption/codding_model/own_basedOnFried/AEJMacro-2015-0289.R2_additional_materials';
% end
% 
 cd(bpath)


%Load parameter values
load('base'); 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iniatlize quantities
% notes sonja:
% pp = calibrated variables using MoM approach 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of simulation periods (period 1=BGP, period 2 =2015-2019,
% period 3 = 2020-2024, period 4 = 2025-2029, period 5 = 2030-2034). Introduce tax in
% period 2
Tsim = 5; 
ind2015 =2; ind2030 = 5;

%Number of experiments: 1 = endog Innovation; 2= exogenous innovation
nExp = 2; 
pp = [alphag, eta, deltaye, epsf, deltaff, nu1];  % parameters from MOM
theta1 = 2.1579; %average from 2001-2010
tauMat = zeros(Tsim, nExp);
tauStarMat = zeros(Tsim, nExp);
noInnov=0; %noInnov =1 =>exogenous innovation

emissionsTarg = 30; %30 percent reduction in emissions
v = zeros(Tsim, 1); %subsidy

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bracketing Routine: Computes emissions for different tax values and 
% solves for the tax for which emissions are 30 percent below their
% baseline value in 2030-2034 time period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:2
    
        ub =5;
        lb = .001;
        diff =1;
        tau = (ub+lb)/2;
        
        while max(abs(diff)) > 1e-10
            
            %Adjust tau star by difference in carbon content between F and 
            %Ostar.
            tauStar = 1.045328173*tau;  
            tauMat(2:end, 2) = tau*ones(Tsim-1, 1);
            tauStarMat(2:end, 2) = tauStar*ones(Tsim-1,1);
            
            a = zeros(Tsim, nExp);
            af = zeros(Tsim, nExp); % machine share fossil fuel production
            ag = zeros(Tsim, nExp); % machine share green energy production
            am = zeros(Tsim, nExp); % machine share  production non-energz
            agg = zeros(Tsim, nExp);
            cHat = zeros(Tsim, nExp); % deviation consumption from BGP
            edHat = zeros(Tsim, nExp); % energy output?
            fdHat = zeros(Tsim, nExp); 
            fdHatStar = zeros(Tsim, nExp);
            fsHat = zeros(Tsim, nExp);
            gdpHat= zeros(Tsim, nExp);
            gdHat = zeros(Tsim, nExp);
            gsHat = zeros(Tsim, nExp);
            lf = zeros(Tsim, nExp); % labour fossil
            lg = zeros(Tsim, nExp); % labour green sector
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
            pfix= zeros(Tsim, nExp); % price machines fossil
            pgix= zeros(Tsim, nExp); % price machines green
            pmix= zeros(Tsim, nExp); % price machines non-energz good
            share= zeros(Tsim, nExp);
            imports= zeros(Tsim, nExp);
            solve= zeros(Tsim, nExp);
            noInnov = 0;
            c1 = 1; 
            c2 = Tsim;
            
            for i = 1:2 %i = 1 => BGP, i => tax
                
                if k ==2 && i ==2
                    noInnov =1;
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
                    = momModel3(pp,a0,alphaf,alpham,deltaef, epse, epsf, epsy, L,noInnov,phi ,rhonf, rhong,...
                    S,tauMat(:,i), tauStarMat(:,i), theta1, n0, guess0, guessBgp, options,p0,Tsim, v);
            end

            sm = S - sf - sg;
            
            %Adjust technology on BGP (momModel3 only computes BGP once, not Tsim times)
            for i = 1:Tsim -1
                a(i+1,1) = (1+n0)*a(i,1);
                af(i+1,1) = (1+n0)*af(i,1);
                ag(i+1,1) = (1+n0)*ag(i,1);
                am(i+1,1) = (1+n0)*am(i,1);
            end
            
           
            if k ==2
                 %No transition dynamics with exogenous inovation. Adjust future
                 %de-trended values to equal values in period 2
                c3 = c2+1;
                c4 = Tsim;
                fdHat(c3:c4,2) = fdHat(2, 2)*ones(size(fdHat(c3:c4)));
                fdHatStar(c3:c4,2) = fdHatStar(2, 2)*ones(size(cHat(c3:c4)));
                %technology with exogenous innovation is the 
                %same as on the BGP
                a(c2:c4, 2)= a(c2:c4, 1); 
            end
            
            %unhat stuff
            fd = fdHat.*a; fdStar = fdHatStar.*a; fs = fsHat.*a;
            
            %Convert energy use to emissions
            d1 = 282.030793*1e9/fd(ind2015);  %MBTU fd 2010/ fd2010 maps fd in model to fd in data
            d3 = 0.072556484/(3.667*1e9);  %giga  tons Carbon per MBTU fd emissions
            d4 =0.07454/(3.667*1e9); %giga tons carbon per MBTU fdStar emissions fdStar per btu
            d5 = 6.187196613;%emissions world per emissions US in 2010
            
            cfd = d1*d3*d5;
            cfdStar = d1*d4*d5;
            emissions = cfd*fd + cfdStar*fdStar;
            
            %Calculate the percent reduction in emissions in 2030-2-34 time
            %period
            deltaEmiss = (emissions(ind2030, 1) - emissions(ind2030,2))./emissions(ind2030, 1)*100;
            diff = deltaEmiss - emissionsTarg;
            
            %Bracket the minimum
            if diff <0 %emissions are too big => raise tau
                lb = tau;
                tau = (ub + tau)/2;
            else ub = tau;
                tau = (lb + tau)/2;
            end
            
        end
        
        diff  
  
    if k ==1
        tauI = tau;
    else tauNI  = tau;
    end
    
end 

%Calculate percent difference in carbon tax from endogenous innovation
tauDiff = (tauI - tauNI)./tauNI*100;
td = -tauDiff; %make absolute value
tauI = tauI; 
tauNI = tauNI; 

tauVec = [tauI, tauNI];
tauStarVec = 1.045328173*tauVec; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert tax to (2013) dollars per ton CO2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tauStarData = tauStarVec./pfStar(1,1)*10.1554 *106.5850/100; %2013 dollars per MBTU
%mean(racOilImportsReal from 2001-2010 is 10.1554)
% 106.5850/100 is the GDP delfator in 2013/ gdp deflator in 2009
tax = tauStarData/0.07454; %2013 dollars per ton CO2 
taxI = tax(1); taxNI = tax(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars ={'td', 'tauI', 'tauNI', 'taxI', 'taxNI'};
save('tax', 'td', 'tauI', 'tauNI', 'taxI', 'taxNI'); 

disp('tax data saved');
