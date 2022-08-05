function [deltaa, E_vec, MOM]= calibration_emissions(T, lengthh, MOM)

% 2019 is base year
% function calibrates emissions and targets

addpath('data')

%-- read in data
% % data: us_ghg_inventory PDF
grosemissions1519_equi= [6.6711, 6.5203, 6.4833, 6.6714, 6.5583]; % gros emissions by year in Gt
grosemissions1519_co2 = [5.3718, 5.2480, 5.2078, 5.3755, 5.2558];
averagecontrib_co2 = mean(grosemissions1519_co2./grosemissions1519_equi); % to calculate sink capacity (equalt share as contributed to gros emissions is used from total sinks)
netemissions1519_equi = [5.9073, 5.6775, 5.7172, 5.8700, 5.7691];
%- read in population data 
datapop =readtable('data/population/unpopulation_dataportal_totalPop.csv', 'ReadVariableNames', true);
%-- us population share

datapop=datapop(datapop.sexId==3,:);

datasmall=datapop(:,["locationId" "timeLabel" "value"]);
datasmall=unstack(datasmall, 'value', 'locationId');
datasmall.popshare=datasmall.x840./datasmall.x900;

%- calibrate sink capacity in US as average over initial period expressed
%  in model periods, i.e. times 5
deltaa         = mean(grosemissions1519_equi-netemissions1519_equi)*averagecontrib_co2*5;  % sinks in the model given in Gt here from equivalents

%--- emissions and target emissions
% emissionsUS2019     = grosemissions1519(end); 
netemissionsUS2019  = grosemissions1519_co2(end)-deltaa/5; % annual 
netemissionsGlobal2019 = 0.75*59; % from IPCC report table a page 6 SPM


%-- emission limits based on equital per capity distribution of limiting emissions
globaltarget_annual30s   = netemissionsGlobal2019*0.5; % annual 
popshare_30to50          = datasmall(2035<=datasmall.timeLabel & datasmall.timeLabel<2050,:);
USnetlimit_annual30s_dynamic = globaltarget_annual30s.*popshare_30to50.popshare;
USnetlimit_modelperiods  = [sum(USnetlimit_annual30s_dynamic(1:5)), sum(USnetlimit_annual30s_dynamic(6:10)), sum(USnetlimit_annual30s_dynamic(11:15))]; %, sum(USnetlimit_annual30s_dynamic(16:20))];

%-- emissin limits
%new report 2022: reducing emissions by 50\% compared to 2019
% nettargetUS_annual30s    = netemissionsUS2019*0.5; 
% nettargetUS_period30s    = nettargetUS_annual30s*5;
globalcarbonBudget20_50  = 510; 
if globalcarbonBudget20_50- globaltarget_annual30s*5*3<globalcarbonBudget20_50/2
    USnetlimit_modelperiods = USnetlimit_modelperiods*globaltarget_annual30s*5*3/globalcarbonBudget20_50;
    globalcarbonBudget20_30  = globalcarbonBudget20_50/2; % from budget deduct what is consumed in 30s
else
    globalcarbonBudget20_30  = globalcarbonBudget20_50- globaltarget_annual30s*5*3;
end
popshare20_30 = mean(datasmall.popshare(datasmall.timeLabel>=2020 & datasmall.timeLabel<2035 )); 

MOM.US_Budget20_30           = globalcarbonBudget20_30*popshare20_30; % see whether code solves when can budget is for 2 periods

% new share of US gros emissions when target is reached
% shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector      
% starting from 2030s (3rd period in model)
length_firstBu = (T-(length(USnetlimit_modelperiods)+(T*lengthh-(50-20))/lengthh));
E_vec          = [repmat(MOM.US_Budget20_30/length_firstBu, 1, length_firstBu ), USnetlimit_modelperiods, zeros(1,(T*lengthh-(50-20))/lengthh)]; 
                % from t=1 to t=29 E=30; from t=30 (2050) to t=T E=0
%- check if emission budget is satisfied
% (sum(E_vec(1:3)./popshare_30to50.popshare)+MOM.US_Budget20_30/mean(datasmall.popshare(datasmall.timeLabel>=2020 & datasmall.timeLabel<2035)))

%- save gros emissions as average of baseline period expressed in model periods to
%   calibrate omega
MOM.grosemissionsUS2019 = mean(grosemissions1519_co2)*5; 
end