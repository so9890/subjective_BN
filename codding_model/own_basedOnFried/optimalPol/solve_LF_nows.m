function [LF_SIM, pol, FVAL] = solve_LF_nows(T, list, pol, params, Sparams,  symms, x0LF, init, indexx)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0      = [x0LF(list.choice~='ws'), 1,0,0,0]; % initial guess from calibration to baseline period (2015-2020 period)
syms lambdaa gammasf gammasn gammasg real
symms.choice= [symms.choice(list.choice~='ws'), lambdaa, gammasf, gammasn, gammasg ];
list.choice= string(symms.choice);

indexxLF.lab = boolean(zeros(size(list.choice)));
indexxLF.exp = boolean(zeros(size(list.choice)));
indexxLF.sqr = boolean(zeros(size(list.choice)));
indexxLF.oneab = boolean(zeros(size(list.choice)));
indexxLF.sci = boolean(zeros(size(list.choice)));

indexxLF.lab(list.choice=='hl'| list.choice=='hh')=1;
indexxLF.exp(list.choice~='sff'& list.choice~='sn' & list.choice~='sg'&list.choice~='hl'& list.choice~='hh' & list.choice~='gammall'& list.choice~='gammalh'& list.choice~='gammasg'& list.choice~='gammasn' & list.choice~='gammasf' )=1;
indexxLF.sqr(list.choice=='gammall'| list.choice=='gammalh'| list.choice=='gammasg'| list.choice=='gammasn' | list.choice=='gammasf' )=1;
indexxLF.sci(list.choice=='sff'& list.choice=='sn' & list.choice=='sg')=1;

% initialise stuff
%-- to save results
symms.allvars=[symms.allvars, gammasn, gammasg, gammasf];
list.allvars=string(symms.allvars);
LF_SIM=zeros(length(list.allvars),T); 
FVAL  = zeros(T,1);
%-- initialise values
laggs   = init; % (init should refer to 2015-2020 period)
t       = 1; % number of periods: t=1: 2020-2024 => does not include base year period (in matrix on first row)
phis    =1;

while t<=T
    fprintf('entering simulation of period %d', t);
    %% - transforming variables to unbounded variables
    %-- index for transformation 

    guess_trans=trans_guess(indexx('LF'), x0, params, list.params);
    guess_trans(indexxLF.sci)=log((params(list.params=='S')-x0(indexxLF.sci))./x0(indexxLF.sci));

    % test
    f=laissez_faire_nows(guess_trans, params, list, pol, laggs, phis);

    %% - solving model

    modFF = @(x)laissez_faire_nows(x, params, list, pol, laggs, phis);
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol, fval, exitf] = fsolve(modFF, guess_trans, options);

   options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5, 'Algorithm', 'trust-region');%, );%, );%, 'Display', 'Iter', );
   [sol2, fval, exitf] = fsolve(modFF, sol, options);

    % pass to standard algorithm
    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);% 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol3, fval, exitf] = fsolve(modFF, sol2, options);

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