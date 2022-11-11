function [deltaa, E_vec, MOM]= calibration_emissions(T, lengthh, MOM)

% 2019 is base year
% function calibrates emissions and targets

% data: us_ghg_inventory PDF
% uses data from 
%--- emissions and target emissions
emissionsUS2019             = 6.55835; 
% emissionsGlobal2019         = 34.2; % SOurce?
% share_US2019                = emissionsUS2019/emissionsGlobal2019;

% required reductions
% targetGlobal    = 30; % ASSUME for now: this is gros and not after deduction of sinks
% reductionGlobal = emissionsGlobal2019-targetGlobal;
% reductionUS     = reductionGlobal*share_US2019; % assume countries have to reduce the proportion which they contribute
% targetUS        = emissionsUS2019-reductionUS; 

% test if this amounts to
% sinks
netEmissionsUS2019  = 5.7691;

%new report 2022: reducing emissions by 50\% compared to 2019

deltaa         = emissionsUS2019-netEmissionsUS2019;  % sinks in the model given in Gt
nettargetUS    = emissionsUS2019*0.5-deltaa; % I assume the reduction to be related to gros emissions, so that the target can be higher by delta

% new share of US gros emissions when target is reached
% shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector      
% starting from 2030s (3rd period in model)
E_vec          = [repmat(nettargetUS,1, (50-30)/lengthh), zeros(1,(T*lengthh-(50-20))/lengthh)]; 
                % assume this is without sinks
                % emissions in Gt (so far, globally)
                % from t=1 to t=29 E=30; from t=30 (2050) to t=T E=0

MOM.emissionsUS2019 = emissionsUS2019; 

end