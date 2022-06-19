function [LF_SIM, pol, FVAL, indexx] = solve_LF_nows_non(T, list, pol, params, Sparams,  symms, init, indexx, indic, Sall)
% simulate economy under laissez faire

% input: 
% pol: numeric policy vector

% output
% LF_SIM: matrix of simulated results in LF, rows= time period, column =
% variables as ordered in list.allvars

%% first run version with only one skill type
syms F pf Lg Af Ag sff sg wsf wsg gammalh gammasg gammasf real
symms.choice_small = [F pf Lg Af Ag sff sg wsf wsg gammalh gammasg gammasf];
list.choice_small=string(symms.choice_small);

%- indexx for transformation

indexxLF.lab = boolean(zeros(size(list.choice_small)));
indexxLF.exp = boolean(zeros(size(list.choice_small)));
indexxLF.sqr = boolean(zeros(size(list.choice_small)));
indexxLF.oneab = boolean(zeros(size(list.choice_small)));

indexxLF.lab(  list.choice_small=='sff'  | list.choice_small=='sg')=1;
indexxLF.exp(list.choice_small~='sff'&list.choice_small~='sg'...
    &list.choice_small~='gammalh'&list.choice_small~='gammasg'& list.choice_small~='gammasf')=1;
indexxLF.sqr(list.choice_small=='gammalh'| list.choice_small=='gammasg'|list.choice_small=='gammasf' )=1;

indexx('LF_noneutral_sep_noskill')=indexxLF;

%- read in solution as starting values
if indic.count_techgap==0
    helper=load('noskill_noneutral_growth_orgtechgap_1706', 'LF');
else
    helper=load('noskill_noneutral_growth_alttechgap_1706', 'LF');
end

F = helper.LF(list.choice_small=='F');
pf = helper.LF(list.choice_small=='pf');
Af = helper.LF(list.choice_small=='Af');
Ag = helper.LF(list.choice_small=='Ag');
sff = helper.LF(list.choice_small=='sff');
sg = helper.LF(list.choice_small=='sg');
wsf = helper.LF(list.choice_small=='wsf');
wsg = helper.LF(list.choice_small=='wsg');
Lg = helper.LF(list.choice_small=='Lg');
gammalh =helper.LF(list.choice_small=='gammalh');
gammasg = helper.LF(list.choice_small=='gammasg');
gammasf = helper.LF(list.choice_small=='gammasf');

%- starting from calibrated model
F = x0LF(list.choice=='F');
pf = x0LF(list.choice=='pf');
Af = x0LF(list.choice=='Af');
Ag = x0LF(list.choice=='Ag');
sff = x0LF(list.choice=='sff');
sg = x0LF(list.choice=='sg');
wsf = x0LF(list.choice=='ws');
wsg = x0LF(list.choice=='ws');
alphag=Sparams.alphag;
Lg = x0LF(list.choice=='G')./(Ag*(x0LF(list.choice=='pg').*alphag).^(alphag./(1-alphag)));
gammalh = x0LF(list.choice=='gammalh');
gammasg = x0LF(list.choice=='gammas');
gammasf = x0LF(list.choice=='gammas');

x0=eval(symms.choice_small);

%- initialise storing stuff
LF_SIM = zeros(length(list.allvars),T); 
FVAL   = zeros(T,1);

%-- initialise values
laggs   = init; % (init should refer to 2010-2014 period)
t       = 1; % number of periods: t=1: 2015-2019 => does include base year period (in matrix on first row) but dont save!

while t<=T+1 % because first iteration is base year
    fprintf('entering simulation of period %d', t);

    % transform data
    guess_trans=trans_guess(indexx('LF_noneutral_sep_noskill'), x0, params, list.params);
    % test
 f=laissez_faire_nows_sep_non_noskillSmall(guess_trans, params, list, pol, init, indic);

%     %- solve
%      lb=[];
%      ub=[];
%      objf=@(x)objectiveCALIBSCI(x);
%      
%      constrf = @(x)laissez_faire_nows_fmincon_sep(x, params, list, pol, init, indic);
%      
%     options = optimset('algorithm','active-set','TolCon', 1e-8,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     [sol3,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constrf,options);

    modFF = @(x)laissez_faire_nows_sep_non_noskillSmall(x, params, list, pol, laggs, indic);
    options = optimoptions('fsolve', 'TolFun', 10e-8, 'MaxFunEvals',8e4, 'MaxIter', 5e5,  'Algorithm', 'levenberg-marquardt', 'Display', 'Iter');%, );%, );%, );
    [sol2, fval, exitf] = fsolve(modFF, guess_trans, options);

    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e6, 'MaxIter', 5e5, 'Display', 'Iter');%, );%, );%, );
    [sol2, fval, exitf] = fsolve(modFF, sol2, options);

    % transform
    LF=trans_allo_out(indexx('LF_noneutral_sep_noskill'), sol2, params, list.params, indic);
    cell_par=arrayfun(@char, symms.choice_small, 'uniform', 0);
    SLF=cell2struct(num2cell(LF), cell_par, 2);

    if t>1
        LF_SIM(:,t-1)=aux_solutionLF_sep_non_small(Sparams, SLF,pol, laggs, list, symms, indexx, params, indic);
        FVAL(t-1)=max(abs(fval));
    end
    %% - update for next round
    x0 = LF; % initial guess
        Af0= SLF.Af; % today's technology becomes tomorrow's lagged technology
        Ag0= SLF.Ag; 
        An0 =zeros(size(Ag0));
    laggs=eval(symms.init);
    t=t+1;
end
end