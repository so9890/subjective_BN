%function [An0, Af0, Ag0, thetaf, thetan, thetag, lambdaa]= calibration_matching(MOM, symms, list)
% function to match moments to model equations
% i.e. solving model plus additional equations for paramters

% To be chosen: 
% An0, Af0, Ag0 (2019)
% thetaf, thetan, thetag : share of high skill in labour input 
% lambdaa to have gov budget = 0 in baseyear

% data
% Energy consumption share, GDP, Fossil to Green output share

% numeraire: Consumption Good

% balanced budget, 
% skills: match Consoli; or wage premia (then includes)

%% initial conditions

%- choiceCALIB variabes
syms hhf hhg hhn hln hlf hlg C F G Af Ag An ...
     hl hh sf sg sn wh wl ws pg pn pe pf gammalh gammall Y...
     Af0 Ag0 An0 thetan thetag thetaf lambdaa real

symms.choiceCALIB = [hhf, hhg, hhn, hln, hlf, hlg, C, F, G, ...
    Af, Ag, An, hl, hh,  sf, sg, sn, wh, wl, ws, pg, pn,...
    pe, pf, gammalh, gammall, Y, Af0, Ag0, An0, thetan, thetag, thetaf, lambdaa];

list.choiceCALIB  = string(symms.choiceCALIB);

%initial values
Af0     = 1.877; % Fried 
Ag0     = 0.9196;
An0     = 1; 
% for init code
syms Ag_lag Af_lag An_lag real
symms.laggs = [Ag_lag, Af_lag, An_lag];
list.laggs= string(symms.laggs);
laggs=[Ag0, Af0, An0];

thetan   = 0.3;
thetag   = 0.4;
thetaf   = thetag*0.5;
lambdaa  = 1; 
hhf =.02; % hhf
hhg =.04; % hhg
% hhn =.4; % hhn
% hln =.3; % hln
hlf =.02; % hlf
hlg =.01; % hlg 
F   = 2; % F
G   = F/MOM.FG; % G
Af  = Af0*1.02; % Af
Ag  = Ag0*1.002; % Ag
An  = An0*1.02; % An
hl  = 0.3; % hl
hh  = 0.2; % hh
% sg  = 0.2;
% sf  = 0.3;
% sn  = 0.4;
gammalh = 0;
gammall = 0;

paramss = eval(symms.params);
poll    = eval(symms.pol);
[x0,ni, checkk]=init(Af, An, Ag, hhg, hhf, hlg, hlf, hh, hl, F, G,  gammalh, gammall,paramss, list, poll, laggs, symms.choiceCALIB);


%% - transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choiceCALIB)));
indexx.exp = boolean(zeros(size(list.choiceCALIB)));
indexx.sqr = boolean(zeros(size(list.choiceCALIB)));
indexx.oneab = boolean(zeros(size(list.choiceCALIB)));

indexx.lab(list.choiceCALIB=='hl'| list.choiceCALIB=='hh')=1;
indexx.exp(list.choiceCALIB~='hl'& list.choiceCALIB~='hh' & list.choiceCALIB~='gammall'&...
    list.choiceCALIB~='gammalh'& list.choiceCALIB~='thetan' & list.choiceCALIB~='thetag'& ...
    list.choiceCALIB~='thetaf' )=1;
indexx.sqr(list.choiceCALIB=='gammall'| list.choiceCALIB=='gammalh')=1;
indexx.oneab(list.choiceCALIB=='thetan'| list.choiceCALIB=='thetag'| list.choiceCALIB=='thetaf') = 1;

%-- transform
guess_trans=trans_guess(indexx, x0, paramss, list);

%- test
f=target_equ(guess_trans, MOM, paramss, list, poll);

%% - solving model

modFF = @(x)target_equ(x, MOM, paramss, list, poll);
options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );

exitfl=0;
count=0;
countmax=100;
save=zeros(2,countmax);

while exitfl<=0 && count<countmax
    count=count+1;
    [sol, fval, exitf] = fsolve(modFF, guess_trans, options);
    guess_trans=sol;
    save(1,count)=max(fval);
    save(2,count)=exitf;
end


