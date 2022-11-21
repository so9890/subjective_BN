function [deltaa, E_vec, MOM, StatsEms]= calibration_emissions(T, lengthh, MOM)

% 2019 is base year
% function calibrates emissions and targets

addpath('data')

%-- read in data
% % data: us_ghg_inventory PDF
grosemissions1519_equi= [6.6711, 6.5203, 6.4833, 6.6714, 6.5583]; % gros emissions by year in giga tons (in report given in  million metric tons, therefore divided by 1,000)
grosemissions1519_co2 = [5.3718, 5.2480, 5.2078, 5.3755, 5.2558];
averagecontrib_co2 = mean(grosemissions1519_co2./grosemissions1519_equi); % to calculate sink capacity (equalt share as contributed to gros emissions is used from total sinks)
netemissions1519_equi = [5.9073, 5.6775, 5.7172, 5.8700, 5.7691];
%- for comparison to politcal target (Biden)
grosCO2_2005 = 6.1345;
grostotal_ems2005 = 7.423;
shareCO2 =grosCO2_2005/grostotal_ems2005;
net_total2005 = 6.635;
sink_2005 =grostotal_ems2005-net_total2005;
sinkCO2_2005 = sink_2005*shareCO2;
netCO2_2005 = grosCO2_2005-sinkCO2_2005;
%- read in population data 
datapop =readtable('data/population/unpopulation_dataportal_totalPop.csv', 'ReadVariableNames', true);

%-- us population share

datapop=datapop(datapop.sexId==3,:);
datasmall=datapop(:,["locationId" "timeLabel" "value"]);
datasmall=unstack(datasmall, 'value', 'locationId');
datasmall.popshare=datasmall.x840./datasmall.x900; % get population share of US

%- calibrate sink capacity in US as average over initial period expressed
%  in model periods, i.e. times 5
deltaa         = sum(grosemissions1519_equi-netemissions1519_equi)*averagecontrib_co2;  % sinks in the model given in Gt here from equivalents

%--- emissions and target emissions
% emissionsUS2019     = grosemissions1519(end); 
netemissionsUS2019  = grosemissions1519_co2(end)-deltaa/5; % annual 
StatsEms.netems_sum1519  =sum(grosemissions1519_co2)-deltaa;
% netemlimit35_equalratio = netemissionsUS2019*0.5*5
netemissionsGlobal2019 = 0.75*59; % from IPCC report table a page 6 SPM
StatsEms.contribUS2019=netemissionsUS2019/netemissionsGlobal2019; % => in 2019 the US contributed 10.44% to global net emissions; 
               % => this is twice as much as its population share! which was equal to 4.3 percent in 2019
               % => even without a new limit the US would have to reduce
               % emissions if wanted to have equitable capital emissions

%-- emission limits based on equital per capity distribution of limiting emissions
globaltarget_annual30s   = netemissionsGlobal2019*0.5; % annual 
popshare_35to50          = datasmall(2035<=datasmall.timeLabel & datasmall.timeLabel<2050,:);
USnetlimit_annual30s_dynamic = globaltarget_annual30s.*popshare_35to50.popshare; % each country can contribute as much as its population share to global target
USnetlimit_modelperiods  = [sum(USnetlimit_annual30s_dynamic(1:5)), sum(USnetlimit_annual30s_dynamic(6:10)), sum(USnetlimit_annual30s_dynamic(11:15))]; %, sum(USnetlimit_annual30s_dynamic(16:20))];

%-- emission limits
% new report 2022: reducing emissions by 50\% compared to 2019
% nettargetUS_annual30s    = netemissionsUS2019*0.5; 
% nettargetUS_period30s    = nettargetUS_annual30s*5;
globalcarbonBudget20_50  = 510; 
if globalcarbonBudget20_50- globaltarget_annual30s*5*3<globalcarbonBudget20_50/2 % what remains for 20-35 period is smaller than half
    % then use equal share of remaining carbon budget in 2020-2035 and
    % 2035-2050 period
%     USnetlimit_modelperiods = USnetlimit_modelperiods*globaltarget_annual30s*5*3/globalcarbonBudget20_50;
    globaltarget_annual30s   = globalcarbonBudget20_50/(2*15); % division by 15 to get annual limit
    USnetlimit_annual30s_dynamic = globaltarget_annual30s.*popshare_35to50.popshare;
    USnetlimit_modelperiods  = [sum(USnetlimit_annual30s_dynamic(1:5)), sum(USnetlimit_annual30s_dynamic(6:10)), sum(USnetlimit_annual30s_dynamic(11:15))]; %, sum(USnetlimit_annual30s_dynamic(16:20))];

    globalcarbonBudget20_35  = globalcarbonBudget20_50/2; % from budget deduct what is consumed in 30s
else
    globalcarbonBudget20_35  = globalcarbonBudget20_50- globaltarget_annual30s*5*3;
end
global_budget_peryear20_30   = globalcarbonBudget20_35/15;
popshare_20to35          = datasmall(2020<=datasmall.timeLabel & datasmall.timeLabel<2035,:);
USnetlimit_annual2035_dynamic = global_budget_peryear20_30.*popshare_20to35.popshare;
USnetlimit_modelperiods2035  = [sum(USnetlimit_annual2035_dynamic(1:5)), sum(USnetlimit_annual2035_dynamic(6:10)),...
    sum(USnetlimit_annual2035_dynamic(11:15))]; %, sum(USnetlimit_annual30s_dynamic(16:20))];

% popshare20_30 = mean(datasmall.popshare(datasmall.timeLabel>=2020 & datasmall.timeLabel<2035 )); 

% MOM.US_Budget20_30           = globalcarbonBudget20_35*popshare20_30; % see whether code solves when can budget is for 2 periods

% new share of US gros emissions when target is reached
% shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector      
% starting from 2030s (3rd period in model)
% length_firstBu = (T-(length(USnetlimit_modelperiods)+(T*lengthh-(50-20))/lengthh));
E_vec          = [USnetlimit_modelperiods2035, USnetlimit_modelperiods, zeros(1,(T*lengthh-(50-20))/lengthh)]; 
                % from t=1 to t=29 E=30; from t=30 (2050) to t=T E=0

%- save gros emissions as average of baseline period expressed in model periods to
%   calibrate omega
MOM.grosemissionsUS2019 = sum(grosemissions1519_co2); 

%-- create statistic of how emission limit relates to 2019 net emissions
StatsEms.perRedu_relato2019_20to50 =1-[USnetlimit_annual2035_dynamic;USnetlimit_annual30s_dynamic]./netemissionsUS2019;
StatsEms.reduc5years =1-E_vec/(MOM.grosemissionsUS2019-deltaa);
StatsEms.USBudget_fromlimit = sum(E_vec);
StatsEms.BudgetUS_popshare = globalcarbonBudget20_50 * mean(datasmall(2020<=datasmall.timeLabel & datasmall.timeLabel<2050,:).popshare);
StatsEms.BAU_ems2050 = netemissionsUS2019*30;
StatsEms.BAU_shareglobBudget = StatsEms.BAU_ems2050./globalcarbonBudget20_50; % => 30\% but allowed: ~4\%
StatsEms.popshare=datasmall;

StatsEms.targetBidenGTCO2 =netCO2_2005*0.52; % for one year
StatsEms.targetBidenGTCO2_modelperiod=5*StatsEms.targetBidenGTCO2;
StatsEms.polRed_rel2019 = 1-StatsEms.targetBidenGTCO2/netemissionsUS2019;
StatsEms.totalnetEms_Biden = netemissionsUS2019*10+20*StatsEms.targetBidenGTCO2;
StatsEms.shareEquCapBudget = StatsEms.totalnetEms_Biden/StatsEms.BudgetUS_popshare;
StatsEms.Emslimit_constantEmsRat = [0.5*netemissionsUS2019*5*ones(1,4), zeros(1,6)];
StatsEms.totalUS_equalRed=20*0.5*netemissionsUS2019;
StatsEms.budgetUS_share_ifequalRed = StatsEms.totalUS_equalRed/globalcarbonBudget20_50;
StatsEms.Emslimit_constantEmsRat_Budget= [globalcarbonBudget20_50/30*5*StatsEms.contribUS2019*ones(1,6), zeros(1,6)];
StatsEms.RedPer_constantEmsRate_Budget = (globalcarbonBudget20_50/30*5*StatsEms.contribUS2019-netemissionsUS2019*5)/(netemissionsUS2019*5);
% globalRed=netemissionsGlobal2019*0.5;
% GDPshareUS1519=0.34; % from file GDP_USDOllar_World
% StatsEms.Emslimit_capabilityApp =[(netemissionsUS2019-globalRed*GDPshareUS1519)*5*ones(1,4), zeros(1,6)];

end