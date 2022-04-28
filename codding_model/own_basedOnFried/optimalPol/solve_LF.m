function [LF_SIM, pol, FVAL] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init, indexx)
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
x0      = x0LF; % initial guess from calibration to baseline period (2015-2020 period)
t       = 1; % number of periods: t=1: 2020-2024 => does not include base year period (in matrix on first row)

while t<=T
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 

    guess_trans=trans_guess(indexx('LF'), x0, params, list.params);

    % test
    % f=laissez_faire(guess_trans, params, list, pol, laggs);

    %% - solving model

    modFF = @(x)laissez_faire(x, params, list, pol, laggs);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
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