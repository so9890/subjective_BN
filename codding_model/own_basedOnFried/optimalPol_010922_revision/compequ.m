function LF_COUNT=compequ(T, list, params, init201519,   symms, LF_SIM, indic, Ems, MOM)
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
        LF_SIM(:,list.allvars=='taul')=zeros(size(LF_SIM(:,list.allvars=='taul')));
        LF_SIM(:,list.allvars=='taus')=zeros(size(LF_SIM(:,list.allvars=='taus')));
    elseif indic.tauf==0
        % version with tauf set to zero to study effect of taul allone
        LF_SIM(:,list.allvars=='tauf')=zeros(size(LF_SIM(:,list.allvars=='tauf')));
        LF_SIM(:,list.allvars=='taus')=zeros(size(LF_SIM(:,list.allvars=='taus')));    
    % for indic.tauf==2 there is no adjustment! 
    % useful when using optimal policy in model with xgr or nsk
    elseif indic.tauf==3 % uses benchmark model policy but sets tauf to zero
        LF_SIM(:,list.allvars=='tauf')=zeros(size(LF_SIM(:,list.allvars=='tauf')));
        LF_SIM(:,list.allvars=='taus')=zeros(size(LF_SIM(:,list.allvars=='taus'))); 
    elseif indic.tauf==4 % uses benchmark model policy but sets tauf to zero and runs model with no know spils
        indic.noknow_spill=0; % effect of opt pol in model without kn spil in full model
    elseif indic.tauf ==5
        indic.limit_LF=1;
        taul= LF_SIM(:,list.allvars=='taul'); % save taul fixed
        % read in new variables as starting values
            helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
             indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
          LF_SIM=helper.opt_all;
        % replace taul
        LF_SIM(:,list.allvars=='taul')=taul; 
    elseif indic.tauf==6
        indic.limit_LF=0; % => in contrast to 5 tauf is given by joint optimal
        taul= LF_SIM(:,list.allvars=='taul'); % save taul fixed
        % read in joint optimal to fix tauf
        helper=load(sprintf('OPT_target_plus30_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
         indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
         LF_SIM=helper.opt_all;
        % replace taul
        LF_SIM(:,list.allvars=='taul')=taul; 
    end

helper.LF_SIM=LF_SIM';

LF_COUNT=solve_LF_VECT(T, list, params,symms, init201519, helper, indic, Ems, MOM);
       
end

