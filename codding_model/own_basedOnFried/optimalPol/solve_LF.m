%function [LF, list, symms] = solve_LF(pol, list, params, targets, Ems )
%solves Laissez faire allocation taking policy as given 

%- lags relevant for equilibrium
syms Ag_lag Af_lag An_lag real
symms.laggs = [Ag_lag, Af_lag, An_lag];
list.laggs= string(symms.laggs);

%-- initialise laggs 
Ag0 = params(list.params=='Ag0');
Af0 = params(list.params=='Af0');
An0 = params(list.params=='An0');
laggs=[Ag0, Af0, An0];

%- choice variabes
syms hhf hhg hlf hlg C F G Af Ag An hl hh real

symms.choice = [hhf, hhg, hlf, hlg, C, F, G, Af, Ag, An, hl, hh];
list.choice  = string(symms.choice);

hhf =.02; % hhf
hhg =.04; % hhg
hlf =.03; % hlf
hlg =.01; % hlg 
C   = 0.2;  % C
F   = 0.3; % F
G   = 0.02; % G
Af  = Af0*1.02; % Af
Ag  = Ag0*1.02; % Ag
An  = An0*1.02; % An
hl  = 0.3; % hl
hh  = 0.2; % hh

read_in_params;
read_in_pol;
auxiliary_stuff;
x0=eval(symms.choice);

%- transforming variables to unbounded variables
%-- index for transformation 
indexx.lab = boolean(zeros(size(list.choice)));
indexx.exp = boolean(zeros(size(list.choice)));

indexx.lab(list.choice=='hl'| list.choice=='hh')=1;
indexx.exp(list.choice~='hl'& list.choice~='hh')=1;

guess_trans=trans_guess(indexx, x0, params, list);

%- solving model
f=laissez_faire(guess_trans, params, list, pol, laggs, targets, Ems);



modFF = @(x)laissez_faire(x, params, list, pol, laggs, targets, Ems);
options = optimoptions('fsolve', 'MaxFunEvals',8e5, 'MaxIter', 3e5, 'TolFun', 10e-12);%, 'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
[sol, fval, exitf] = fsolve(modFF, guess_trans, options);

% run ones more with solution as starting value
[sol2, fval, exitf] = fsolve(modFF, sol, options);

%- transform results to bounded variables
LF=trans_allo_out(indexx, sol2, params, list);
 hhf    = LF(list.choice=='hhf');
 hhg    = LF(list.choice=='hhg');
 hlf    = LF(list.choice=='hlf');
 hlg    = LF(list.choice=='hlg');
 C      = LF(list.choice=='C');
 F      = LF(list.choice=='F');
 G      = LF(list.choice=='G');
 Af     = LF(list.choice=='Af');
 Ag     = LF(list.choice=='Ag');
 An     = LF(list.choice=='An');
 hl     = LF(list.choice=='hl');
 hh     = LF(list.choice=='hh');
 
%end
% test market clearing
read_in_params;
read_in_pol;
auxiliary_stuff;
checkk=C+xn+xf+xg-Y