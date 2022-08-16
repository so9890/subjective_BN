function LF_COUNT=compequ(T, list, params, init201519,   symms, LF_SIM, indic)
% Competitive equilibrium with policy optimal without spillovers
% DOES NOT SOLVE WITH ETAA ==1
% for version without emission target solve LF at (taul=0, taus=0, lambdaa=1, tauf=0)
% if Sparams.etaa~=1

% LF_SIM: initial guess

% output
% LF_COUNT: counterfactual competitive equilibrium allocation 

%- read in counterfactual policy

    if indic.tauf==1
        % version with taul set to zero but tauf as optimal
        LF_SIM(:,list.sepallvars=='taul')=zeros(size(LF_SIM(:,list.sepallvars=='taul')));
        LF_SIM(:,list.sepallvars=='taus')=zeros(size(LF_SIM(:,list.sepallvars=='taus')));
    elseif indic.tauf==0
        % version with tauf set to zero to study effect of taul allone
        LF_SIM(:,list.sepallvars=='tauf')=zeros(size(LF_SIM(:,list.sepallvars=='tauf')));
        LF_SIM(:,list.sepallvars=='taus')=zeros(size(LF_SIM(:,list.sepallvars=='taus')));    
    % for indic.tauf==2 there is no adjustment! 
    % useful when using optimal policy in model with xgr or nsk
    elseif indic.tauf==3 % uses benchmark model policy but sets tauf to zero
        LF_SIM(:,list.sepallvars=='tauf')=zeros(size(LF_SIM(:,list.sepallvars=='tauf')));
        LF_SIM(:,list.sepallvars=='taus')=zeros(size(LF_SIM(:,list.sepallvars=='taus')));    
    end

helper.LF_SIM=LF_SIM';

LF_COUNT=solve_LF_VECT(T, list, params,symms, init201519, helper, indic);
       
end

