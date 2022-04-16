function [An0, Af0, Ag0, thetaf, thetan, thetag, lambdaa]= calibration_matching(zh,zl)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear

% data
% Energy consumption share, GDP, Fossil to Green output share

% numeraire: Consumption Good

% balanced budget, 
% skills: match Consoli; or wage premia (then includes)

%% data
MOM = calibration_moments(zh,zl)



