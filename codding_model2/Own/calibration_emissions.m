function [targets_num, E_vec]= calibration_emissions(yd2019, symstargets,T)

% function calibrates emissions and targets
% uses total output in 2019 LF of dirty sector as an input.

%--- emissions and target emissions
emissionsUS2019 = 6.55835; 
emissionsGlobal2019 = 34.2; 
share_US2019= emissionsUS2019/emissionsGlobal2019;

% required reductions
targetGlobal    = 30; % ASSUME for now: this is gros and not after deduction of sinks
reductionGlobal = emissionsGlobal2019-targetGlobal;
reductionUS     = reductionGlobal*share_US2019; % assume countries have to reduce the proportion which they contribute
targetUS        = emissionsUS2019*reductionUS; 

% sinks
netEmissionsUS2019 = 5.7691;
deltaa   = emissionsUS2019-netEmissionsUS2019;  % sinks in the model given in Gt

% new share of US gros emissions when target is reached
shareUSnew  = targetUS/targetGlobal; % => the share reduces by 2 percentage points! 
                        % assume: share of total emissions= share of total
                        % necessary reduction
% emission target us as a vector                        
E_vec          = [repmat(targetUS-deltaa,1, 50-20-1), zeros(1,T-30+1)]; % assume this is without sinks
                        % emissions in Gt (so far, globally)
                        % from t=1 to t=29 E=30; from t=30 to t=T E=0
%E              = targetUS-deltaa;
%--- relation production and emissions  
                    
kappaa   = emissionsUS2019/yd2019; % share of dirty output which translates into emissions 
                                    % => to match total greenhouse gas emissions in Gt

% save                                    
targets_num =eval(symstargets);
end