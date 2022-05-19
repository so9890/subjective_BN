function RAM = solve_sym_sep(symms, list, Ftarget, indic, T, indexx, params)
% function to solve Ramsey problem resulting from symbolic derivation 
% using fsolve
% order of variables in Ram_Model as in list.optALL

%- guess: use sp allocation as starting point
if indic.target==0
    helper=load('SYMOPT_notarget_active_set_0505_spillover1_taus0_noskill0_notaul0_alt.mat');
else
    helper= load('OPT_notarget_active_set_1905_spillover0_taus0_noskill0_notaul0.mat');
end
    opt_all=helper.opt_all;
    x0=symms.optALL;

    x0(startsWith(list.optALL, 'hhf')) =opt_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
    x0(startsWith(list.optALL, 'hhg')) =opt_all(:,list.allvars=='hhg'); % hhg
    x0(startsWith(list.optALL, 'hlf')) =opt_all(:,list.allvars=='hlf'); % hlf
    x0(startsWith(list.optALL, 'hlg')) =opt_all(:,list.allvars=='hlg'); % hlg 
    x0(startsWith(list.optALL, 'C'))   =opt_all(:,list.allvars=='C');   % C
    x0(startsWith(list.optALL, 'F'))   =opt_all(:,list.allvars=='F');

    x0(startsWith(list.optALL, 'G'))   =opt_all(:,list.allvars=='G');   % G
    x0(startsWith(list.optALL, 'Af'))  =opt_all(:,list.allvars=='Af');  % Af
    x0(startsWith(list.optALL, 'Ag'))  =opt_all(:,list.allvars=='Ag');  % Ag
    x0(startsWith(list.optALL, 'An'))  =opt_all(:,list.allvars=='An');  % An
    x0(startsWith(list.optALL, 'HL'))  =opt_all(:,list.allvars=='hl');  % hl
    x0(startsWith(list.optALL, 'HH'))  =opt_all(:,list.allvars=='hh');   % hh
    x0(startsWith(list.optALL, 'sff'))   =opt_all(:,list.allvars=='sff');   % hh
    x0(startsWith(list.optALL, 'sn'))   =opt_all(:,list.allvars=='sn');
    x0(startsWith(list.optALL, 'sg'))   =opt_all(:,list.allvars=='sg');
    % number of multipliers
    nm= sum(startsWith(list.optsym, 'mu_'));
    %nkt= sum(startsWith(list.optsym, 'KT_'));
    x0(startsWith(list.optALL, 'mu_'))       = ones(nm*T,1);   % lagraneg multipliers (for emission target updated later)
    x0(startsWith(list.optALL, 'mu_target')) = 100;
    %x0(startsWith(list.optALL, 'KT_'))       = zeros(nkt*T,1);  
    
    x0=eval(x0);
    
%% -- transform to unbounded variables 
guess_trans=log(x0);

%- exceptions with upper bound; hl, hh, S, and F in case of target

if indic.noskill==0
    guess_trans(startsWith(list.optALL, 'HL'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'HL')))./...
        x0(startsWith(list.optALL, 'HL')));
    guess_trans(startsWith(list.optALL, 'HH'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'HH')))./...
        x0(startsWith(list.optALL, 'HH')));
else
        guess_trans(startsWith(list.optALL, 'H'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'H')))./...
        x0(startsWith(list.optALL, 'H')));
end
% initial value of scientists might be zero (in version with target; other transformation required)
if indic.target==1
    guess_trans(startsWith(list.optALL, 'S'))=sqrt(x0(startsWith(list.optALL, 'S')));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
else
       guess_trans(startsWith(list.optALL, 'sff'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'sff')))./...
        x0(startsWith(list.optALL, 'sff'))); 
      guess_trans(startsWith(list.optALL, 'sn'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'sn')))./...
        x0(startsWith(list.optALL, 'sn'))); 
      guess_trans(startsWith(list.optALL, 'sg'))=log((params(list.params=='upbarH')-x0(startsWith(list.optALL, 'sg')))./...
        x0(startsWith(list.optALL, 'sg'))); 
end

% guess_trans(startsWith(list.optALL, 'sg'))=sqrt(x0(startsWith(list.optALL, 'sg')));

if indic.target==1
    guess_trans(startsWith(list.optALL, 'F'))=log((Ftarget-x0(startsWith(list.optALL, 'F')))./...
     x0(startsWith(list.optALL, 'F')));
end
%- lagrange multipliers not transformed
guess_trans(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'))=x0(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'));

%% -- test model and solve for optimal policy

if indic.target ==1
    model_trans = Ram_Model_target(guess_trans);
    modFF = @(x)Ram_Model_target(x);
else
    model_trans = Ram_Model_notarget_sep_1905(guess_trans);
    modFF = @(x)Ram_Model_notarget_sep_1905(x);
end

options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-10, 'Display', 'Iter', 'Algorithm', 'levenberg-marquardt');%, );%, );%, );

[outt, fval, exitflag] = fsolve(modFF, guess_trans, options);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-10, 'Display', 'Iter'); %, 'Algorithm', 'levenberg-marquardt');%, );%, );%, );
[outt2, fval, exitflag] = fsolve(modFF, outt, options);

options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-10, 'Display', 'Iter'); %, 'Algorithm', 'levenberg-marquardt');%, );%, );%, );
[x, fval, exitflag] = fsolve(modFF, outt2, options);
options = optimoptions('fsolve', 'MaxFunEvals',20e6, 'MaxIter', 3e5, 'TolFun', 10e-10, 'Display', 'Iter');%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, );
[x14, fval, exitflag] = fsolve(modFF, x13, options);


% with fmincon
lb=[];
ub=[];
constf=@(x)sym_sep_fmincon(x);
objf=@(x)objectiveCALIBSCI(x);
options = optimset('algorithm','sqp','TolCon', 1e-11,'Tolfun',1e-26,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
[xsqp,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
ss=load(sprintf('opt_sym_1605_notarget_noskill%d_spillover%d_notaul%d_sep.mat', indic.noskill, indic.spillovers, indic.notaul));
%% transform
upbarH=params(list.params=='upbarH');
%  x=zeros(size(list.optALL));
%  for j =1:length(list.optALL)
%      jj=list.optALL(j);
%      x(list.optALL==string(jj))=ss.x14(ss.list.optALL==string(jj));
%  end
x=outt;
out_trans=exp(x);
if indic.noskill==0
    out_trans(contains(list.optALL,'HL'))=upbarH./(1+exp(x(contains(list.optALL,'HL'))));
 out_trans(contains(list.optALL,'HH'))=upbarH./(1+exp(x(contains(list.optALL,'HH'))));
else
 out_trans(contains(list.optALL,'H'))=upbarH./(1+exp(x(contains(list.optALL,'H'))));
end

if indic.target == 0
 out_trans(startsWith(list.optALL,'sn'))=upbarH./(1+exp(x(startsWith(list.optALL,'sn'))));
 out_trans(startsWith(list.optALL,'sg'))=upbarH./(1+exp(x(startsWith(list.optALL,'sg'))));
 out_trans(startsWith(list.optALL,'sff'))=upbarH./(1+exp(x(startsWith(list.optALL,'sff'))));

else
    out_trans(startsWith(list.optALL,'S')) = (x(startsWith(list.optALL,'S'))).^2;
end
out_trans(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'))=x(contains(list.optALL, 'mu')|contains(list.optALL, 'KT'));
% if indic.taus==1
%     out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
% end
if indic.target==1
    out_trans(startsWith(list.optALL,'F'))=Ftarget./(1+exp(x(startsWith(list.optALL,'F'))));
end

x=out_trans;
read_in_SYMMODEL_sep;
taus = zeros(size(pn));
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
gammasg = zeros(size(pn));
gammasf = zeros(size(pn));
gammasn = zeros(size(pn));

opt_all=eval(symms.sepallvars);

%% - test swf 
%  disc=repmat(betaa, 1,T);
%  expp=0:T-1;
%  vec_discount= disc.^expp;
% vec_discount*SWF
% test if opt solution is a laissez faire solution 
% function throws error if solution is not a solution to LF
helper.LF_SIM=opt_all';
f=test_LF_VECT(T, list,  params,symms, init201519, helper, indic);
 ind=1:length(f);
 ss=ind(abs(f)>1e-6);
 tt=floor(ss/T); 
% save results
if indic.target==1
    save(sprintf('SYMOPT_target_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d', indic.spillovers, indic.taus, indic.noskill, indic.notaul), 'opt_all')
else
    save(sprintf('SYMOPT_notarget_active_set_0505_spillover%d_taus%d_noskill%d_notaul%d_alt', indic.spillovers, indic.taus, indic.noskill, indic.notaul), 'opt_all')
end
end