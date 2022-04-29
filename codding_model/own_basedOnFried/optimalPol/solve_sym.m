function RAM = solve_sym(symms, list, Ftarget, indic)
% function to solve Ramsey problem resulting from symbolic derivation 
% using fsolve
% order of variables in Ram_Model as in list.optALL

%- guess: use sp allocation as starting point
if indic.target==1
    helper=load('SP_target');
    sp_all=helper.sp_all;
    x0=symms.optALL;
   
    x0(startsWith(list.optALL, 'hhf')) =sp_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(startsWith(list.optALL, 'hhg')) =(0.7)*sp_all(:,list.allvars=='hhg'); % hhg
    x0(startsWith(list.optALL, 'hlf')) =sp_all(:,list.allvars=='hlf'); % hlf
    x0(startsWith(list.optALL, 'hlg')) =(0.5)*sp_all(:,list.allvars=='hlg'); % hlg 
    x0(startsWith(list.optALL, 'C'))   =(0.8)*sp_all(:,list.allvars=='C');   % C
    x0(startsWith(list.optALL, 'F'))   =(0.7)*sp_all(:,list.allvars=='F');

    x0(startsWith(list.optALL, 'G'))    =(0.7)*sp_all(:,list.allvars=='G');   % G
    x0(startsWith(list.optALL, 'Af'))   =0.5*sp_all(:,list.allvars=='Af');  % Af
    x0(startsWith(list.optALL, 'Ag'))   =sp_all(:,list.allvars=='Ag');  % Ag
    x0(startsWith(list.optALL, 'An'))   =sp_all(:,list.allvars=='An');  % An
    x0(startsWith(list.optALL, 'HL'))   =sp_all(:,list.allvars=='hl');  % hl
    x0(startsWith(list.optALL, 'HH'))   =sp_all(:,list.allvars=='hh');   % hh

    % number of multipliers
    nm= sum(startsWith(list.optsym, 'mu_'));
    nkt= sum(startsWith(list.optsym, 'KT_'));
    x0(startsWith(list.optALL, 'mu_'))       = ones(nm*T,1);   % lagraneg multipliers (for emission target updated later)
    x0(startsWith(list.optALL, 'mu_target')) = 100;
    x0(startsWith(list.optALL, 'KT_'))       = zeros(nkt*T,1);  
else
    load('OPT_notarget.mat')
    % use result from numeric problem
    
%         taus=0;
%         tauf=0;
%         taul=0;
%         lambdaa=1; % balances budget with tauf= taul=0
%         pol=eval(symms.pol);
% 
%         if ~isfile('FB_LF_SIM_NOTARGET.mat')
%             [LF_SIM] = solve_LF(T, list, pol, params, Sparams,  symms, x0LF, init201519, indexx);   
%             save('FB_LF_SIM_NOTARGET','LF_SIM');
%         else
%             help=load('FB_LF_SIM_NOTARGET.mat');
%             LF_SIM=help.LF_SIM;
%         end
    
    x0 =symms.optALL;
   
    x0(startsWith(list.optALL, 'hhf')) =opt_all(1:T,list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(startsWith(list.optALL, 'hhg')) =opt_all(1:T,list.allvars=='hhg'); % hhg
    x0(startsWith(list.optALL, 'hlf')) =opt_all(1:T,list.allvars=='hlf'); % hlf
    x0(startsWith(list.optALL, 'hlg')) =opt_all(1:T,list.allvars=='hlg'); % hlg 
    x0(startsWith(list.optALL, 'C'))     =opt_all(1:T,list.allvars=='C');   % C
    x0(startsWith(list.optALL, 'F'))     =opt_all(1:T,list.allvars=='F');
    x0(startsWith(list.optALL, 'G'))     =opt_all(1:T,list.allvars=='G');   % G
    x0(startsWith(list.optALL, 'Af'))   =opt_all(1:T,list.allvars=='Af');  % Af
    x0(startsWith(list.optALL, 'Ag'))   =opt_all(1:T,list.allvars=='Ag');  % Ag
    x0(startsWith(list.optALL, 'An'))   =opt_all(1:T,list.allvars=='An');  % An
    x0(startsWith(list.optALL, 'HL'))   =opt_all(1:T,list.allvars=='hl');  % hl
    x0(startsWith(list.optALL, 'HH'))   =opt_all(1:T,list.allvars=='hh');  % hh

    %multipliers
    nm= sum(startsWith(list.optsym, 'mu_')); % number of multipliers
    nkt= sum(startsWith(list.optsym, 'KT_'));
    x0(startsWith(list.optALL, 'mu_'))       = ones(nm*T,1);   % lagraneg multipliers (for emission target updated later)
    x0(startsWith(list.optALL, 'KT_'))       = zeros(nkt*T,1);  
end
    
    x0=eval(x0);
    
%-- transform to unbounded variables 
indexxO=indexx('OPTSym');
guess_trans=trans_guess(indexxO,symms.optALL, params, list.params);
guess_trans=trans_guess(indexxO,x0, params, list.params);

% replace F manually
if indic.target==1
   guess_trans(startsWith(list.optALL, 'F'))= log((Ftarget-x0(startsWith(list.optALL, 'F')))./x0(startsWith(list.optALL, 'F')));
end

%-- test model and solve for optimal policy

if indic.target ==1
    model_trans = Ram_Model_target(guess_trans);
    modFF = @(x)Ram_Model_target(x);
else
    model_trans = Ram_Model_notarget(guess_trans);
    modFF = @(x)Ram_Model_notarget(x);
end

options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-12, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

[outt, fval, exitflag] = fsolve(modFF, guess_trans, options);

xx=real(outt);
[outt, fval, exitflag] = fsolve(modFF, xx, options);
fprintf('ramsey solved with exitflag %d ', exitflag );