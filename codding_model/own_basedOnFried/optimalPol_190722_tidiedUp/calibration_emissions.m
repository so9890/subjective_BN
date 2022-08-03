function [deltaa, E_vec, MOM]= calibration_emissions(T, lengthh, MOM)

% 2019 is base year
% function calibrates emissions and targets


%-- read in data
% % data: us_ghg_inventory PDF
grosemissions1519= [6.6711, 6.5203, 6.4833, 6.6714, 6.5583]; % gros emissions by year in Gt
netemissions1519 = [5.9073, 5.6775, 5.7172, 5.8700, 5.7691];

%--- emissions and target emissions
% emissionsUS2019     = grosemissions1519(end); 
netemissionsUS2019  = netemissions1519(end);

%- calibrate sink capacity in US as average over initial period expressed
%  in model periods, i.e. times 5
deltaa         = mean(grosemissions1519-netemissions1519)*5;  % sinks in the model given in Gt

%-- emissin limits
%new report 2022: reducing emissions by 50\% compared to 2019
nettargetUS_annual30s    = netemissionsUS2019*0.5; 
nettargetUS_period30s    = nettargetUS_annual30s*5;

% new share of US gros emissions when target is reached
% shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector      
% starting from 2030s (3rd period in model)
E_vec          = [repmat(nettargetUS_period30s,1, (50-30)/lengthh), zeros(1,(T*lengthh-(50-20))/lengthh)]; 
                % from t=1 to t=29 E=30; from t=30 (2050) to t=T E=0

%- save gros emissions as average of baseline period expressed in model periods to
%   calibrate omega
MOM.grosemissionsUS2019 = mean(grosemissions1519)*5; 

end