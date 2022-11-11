function [deltaa, E_vec, MOM]= calibration_emissions(T, lengthh, MOM)

% 2019 is base year
% function calibrates emissions and targets

% data: us_ghg_inventory PDF
% uses data from 
%--- emissions and target emissions
emissionsUS2019             = 6.55835; 
emissionsGlobal2019         = 34.2; % SOurce?
share_US2019                = emissionsUS2019/emissionsGlobal2019;

% required reductions
targetGlobal    = 30; % ASSUME for now: this is gros and not after deduction of sinks
reductionGlobal = emissionsGlobal2019-targetGlobal;
reductionUS     = reductionGlobal*share_US2019; % assume countries have to reduce the proportion which they contribute
targetUS        = emissionsUS2019*reductionUS; 

% sinks
netEmissionsUS2019  = 5.7691;
deltaa              = emissionsUS2019-netEmissionsUS2019;  % sinks in the model given in Gt

% new share of US gros emissions when target is reached
shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector                        
E_vec          = [repmat(targetUS-deltaa,1, (50-20)/lengthh-1), zeros(1,(T*lengthh-30)/lengthh+1)]; 
                % assume this is without sinks
                % emissions in Gt (so far, globally)
                % from t=1 to t=29 E=30; from t=30 (2050) to t=T E=0

MOM.emissionsUS2019 = emissionsUS2019; 

end