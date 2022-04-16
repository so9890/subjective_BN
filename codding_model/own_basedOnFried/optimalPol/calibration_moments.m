function MOM = calibration_moments(zh,zl)
% function reads in data and generates moments 
% to be calibrated: lambda, thetaf, thetag, Af0, Ag0, An0
% (thetan for now assumed shares as in population)


addpath('data')

%% - read in data

data =readtable('data/energy_consumption/combined.xlsx');

%- extract data and variables names
    momsraw=data{:,end};
    vars=data{:,5};
    moms=momsraw(~isnan(momsraw));
    vars=string(vars(~isnan(momsraw)));
  
% moments for technology
MOM.EpeY = moms(vars=='Energy Expenditures as Share of GDP');
MOM.Y    = 1; % normalised
MOM.FG   = moms(vars=='Total Fossil Fuels Consumption')./moms(vars=='Total Renewable Energy Consumption');
MOM.F    = moms(vars=='Total Fossil Fuels Consumption'); % to match emissions
% Data on skill input and distribution
MOM.whwl = 1.9; % 
MOM.hhhl = 1; % assume equal hours worked, as in Bick, Schünden, Lagakos
MOM.hhnhln = zh./zl; % assume share of high skill in neutral sector matches skill distribtuion
MOM.diff   = 1.3; % factor on skill share in green sector

% Gov Budget
MOM.Debt  = 0; % assume balanced budget in baseyear

end