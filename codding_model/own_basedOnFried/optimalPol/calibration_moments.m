function MOM = calibration_moments()
% function reads in data and generates moments 
% to be calibrated: lambda, thetaf, thetag, Af0, Ag0, An0
% (thetan for now assumed shares as in population)


addpath('data')

%% - read in data

dataout =readtable('data/energy_consumption/combined.xlsx');
dataskill =readtable('data/skill/employment_sharesConsoli.xlsx');
datahours =readtable('data/hoursworked/ANHRS_03022021112723019.csv', 'ReadVariableNames', true);

%- average hours 2015-2019
MOM.targethour=1/5*(datahours.asShareOfTotalHours(datahours.TIME==2019)+datahours.asShareOfTotalHours(datahours.TIME==2018)+datahours.asShareOfTotalHours(datahours.TIME==2017)+datahours.asShareOfTotalHours(datahours.TIME==2016)+datahours.asShareOfTotalHours(datahours.TIME==2015));
%- extract data and variables names
    momsraw=dataout{:,end};
    vars=dataout{:,5};
    moms=momsraw(~isnan(momsraw));
    vars=string(vars(~isnan(momsraw)));
  
% moments for technology
MOM.EpeY = moms(vars=='Energy Expenditures as Share of GDP')/100;
MOM.Y    = 1; % normalised
MOM.FG   = moms(vars=='Total Fossil Fuels Consumption')./moms(vars=='Total Renewable Energy Consumption');
%MOM.F    = moms(vars=='Total Fossil Fuels Consumption'); % to match emissions
% Data on skill input and distribution
MOM.whwl = 1.9; % 
%MOM.hhhl = 1; % assume equal hours worked, as in Bick, Schünden, Lagakos
%MOM.hhnhln = zh./zl; % assume share of high skill in neutral sector matches skill distribtuion
%MOM.diff   = 1.3; % factor on skill share in green sector

% Gov Budget
MOM.Debt  = 0; % assume balanced budget in baseyear

%% skill distribution
dist_total=dataskill{1:22,6};
dist_greenlower=dataskill{1:22,15};
%dist_greenupper=dataskill{1:22,14};
dist_nongreen=dataskill{1:22,16};

indic_skill=dataskill{1:22,10};
green_total = dataskill{23,7};

%- share high skill
highskill_total = sum(dist_total(indic_skill==1)); 
highskill_green = sum(dist_greenlower(indic_skill==1));  
%highskill_greenUp = sum(dist_greenupper(indic_skill==1));  
highskill_nongreen = sum(dist_nongreen(indic_skill==1));  

%-save moments
MOM.hhg_hhghlg = highskill_green;
MOM.sharehighnongreen = highskill_nongreen;
MOM.hhehzh_total = highskill_total; 
MOM.hg_total   = green_total;

end