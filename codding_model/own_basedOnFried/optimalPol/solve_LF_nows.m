function [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init, indexx)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise stuff
%-- to save results
% symms.allvars=[symms.allvars, gammasn, gammasg, gammasf];
% list.allvars=string(symms.allvars);
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);
x0      = x0LF;

%-- initialise values
laggs   = init; % (init should refer to 2010-2014 period)
t       = 1; % number of periods: t=1: 2015-2019 => does include base year period (in matrix on first row) but dont save!

while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 

    guess_trans=trans_guess(indexx('LF'), x0, params, list.params);

    % test
    f=laissez_faire_nows(guess_trans, params, list, pol, laggs);

    %% - solving model

    modFF = @(x)laissez_faire_nows(x, params, list, pol, laggs);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol, fval, exitf] = fsolve(modFF, guess_trans, options);

   options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'trust-region');%, );%, );%, 'Display', 'Iter', );
   [sol2, fval, exitf] = fsolve(modFF, sol, options);

    % pass to standard algorithm
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol3, fval, exitf] = fsolve(modFF, sol2, options);

    %- transform results to bounded variables
    LF=trans_allo_out(indexx('LF'), sol3, params, list.params);
    
    %% - save results
    % this part also checks correctness of results!
    cell_par=arrayfun(@char, symms.choice, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);
    if t>1
        LF_SIM(:,t)=aux_solutionLF(Sparams, SLF, pol, laggs, list, symms, indexx, params);
        FVAL(t)=max(abs(fval));
    end
    %% - update for next round
    x0 = LF; % initial guess
        Af0= SLF.Af; % today's technology becomes tomorrow's lagged technology
        Ag0= SLF.Ag; 
        An0= SLF.An; 
    laggs=eval(symms.init);
    t=t+1;
end
end