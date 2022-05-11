function [LF_SIM, pol, FVAL] = solve_LF_EMTARGET(T, list, pol, params, Sparams,  symms, x0LF, init, indexx, EMS)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise stuff
%-- to save results
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);
%-- initialise values
laggs   = init; % (init should refer to 2015-2020 period)
x0      = x0LF(list.choice~='F'); % initial guess from calibration to baseline period (2015-2020 period)
t       = 1; % number of periods: t=1: 2020-2024 => does not include base year period (in matrix on first row)

%%
while t<=T
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 

    guess_trans=trans_guess(indexx('EXPI'), x0, params, list.params);
    indic.tauffixed=0;
    
    % now add guess for policy variable, fixed
    if indic.tauffixed==1
        guess_trans=[guess_trans, log(1-0.8)];
        fixedvalue=0.9; 
    else
        guess_trans=[guess_trans, log((1-0.02)./0.02)];
        fixedvalue=0; 
    end
    % test
     f=laissez_faire_EMTARGET(guess_trans, params, list, pol, laggs, Ems(t), fixedvalue, indic);

    %% - solving model

    modFF = @(x)laissez_faire_EMTARGET(guess_trans, params, list, pol, laggs, Ems(t), fixedvalue, indic);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol, fval, exitf] = fsolve(modFF, guess_trans, options);

%     options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'trust-region');%, );%, );%, 'Display', 'Iter', );
%     [sol, fval, exitf] = fsolve(modFF, guess_trans, options);

    % pass to standard algorithm
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol, fval, exitf] = fsolve(modFF, sol, options);

    %- transform results to bounded variables
    LF=trans_allo_out(indexx('LF'), sol, params, list.params);
    
    %% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    LF_SIM(:,t)=aux_solutionLF(Sparams, SLF, pol, laggs, list, symms, indexx, params);
    FVAL(t)=max(abs(fval));
    %% - update for next round
    x0 = LF; % initial guess
        Af0= SLF.Af; % today's technology becomes tomorrow's lagged technology
        Ag0= SLF.Ag; 
        An0= SLF.An; 
    laggs=eval(symms.init);
    t=t+1;
end
end