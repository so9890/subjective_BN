function [symms, list, opt_all]= OPT_solve_sep(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;

if indic.notaul==9 % taking labor tax from optimal policy without target
    helper  = load(sprintf('OPT_notarget_2112_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul4_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.sep, indic.extern , indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    taulFixed           =helper.opt_all(:,list.allvars=='taul');
else
    taulFixed=0;
end
% symbilic variables and lists
syms hhf hhg hlf hlg hhn hln C Ch Cl F G Af Ag An hl hh S sff sn sg Lf Lg h gammasg gammasf real

if indic.notaul >=7
        symms.opt = [hhf hhg hlf hlg hhn hln C F G hl hh sn sff sg S];    
else
    if indic.sep==1
        if indic.noskill==0
            symms.opt = [hhf hhg hlf hlg hhn hln C F G hl hh sn sff sg];
        else
                symms.opt = [Lf Lg C F G h sn sff sg];
        end
    elseif indic.sep==0
        if indic.noskill==0
        symms.opt = [hhf hhg hlf hlg hhn hln C F G hl hh sn sff sg S];
    else
            symms.opt = [Lf Lg C F G h sn sff sg S];
        end
    end
        
end

if indic.xgrowth==1
    if indic.noskill==0
        symms.opt= [hhf hhg hlf hlg hhn hln C F G hl hh];
    else
        symms.opt= [Lf Lg C F G h];
    end
end

list.opt  = string(symms.opt); 

nn= length(list.opt); % number of variables
%%
% if indic.sep==1   
%     list.allvars=list.sepallvars;
%     symms.allvars=symms.sepallvars;
% end

    % use correct list to read in variables!
%      saved=list.allvars;
%      hhelper= load(sprintf('params_0209_sep%d', 1));
% % 
%      list.allvars=hhelper.list.allvars;
%% Initial Guess %%% 
%%%%%%%%%%%%%%%%%%%%%
if indic.target==1
     
  %  helper= load(sprintf('OPT_notarget_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill,indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
  if indic.sigmaWorker~=0
    helper=load(sprintf('OPT_target_0509_sigmaW%d_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.sigmaWorker, indic.noknow_spill ,indic.noskill,5, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));
  elseif indic.Bop~=0
    helper=load(sprintf('OPT_target_0512_Bop%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Bop,indic.Sun, indic.noknow_spill ,indic.noskill,5, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));
  elseif indic.count_techgap==1
      helper=load(sprintf('OPT_target_countec_0509_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.noknow_spill ,indic.noskill,5, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));
   elseif   indic.notaul~=9
%        helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover0_knspil3_taus0_noskill%d_notaul4_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun,  plotts.nsk, indic.sep, plotts.xgr, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
       helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.noknow_spill,indic.noskill, 5, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));
  elseif indic.notaul==9
      helper=load(sprintf('COMPEquN_SIM_0501_taufopt%d_newCalib_target_emnet%d_Sun%d_knspil%d_spillover%d_notaul%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', 5,indic.targetWhat, indic.Sun, indic.noknow_spill,  indic.spillovers, 4, indic.noskill, indic.sep, indic.xgrowth, indic.PV,etaa));
      helper.opt_all=helper.LF_COUNT;
  else
      helper=load(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.noknow_spill ,indic.noskill, indic.notaul, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));

%    helper=load(sprintf('OPT_target_0509_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.noknow_spill ,indic.noskill, indic.notaul, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));
  end
        opt_all=helper.opt_all;
%     kappaa = [repmat(Ftarget(1), 1,percon) ,Ftarget']./opt_all(1:T,list.allvars=='F')'; % ratio of targeted F to non-emission
kappaa= Ftarget'./opt_all(1:T,list.allvars=='F')';    
kappaa = kappaa*(1-1e-10);
    
        x0 = zeros(nn*T,1);
    if indic.noskill==0
        x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =opt_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =opt_all(:,list.allvars=='hhg'); % hhg
        x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =opt_all(:,list.allvars=='hlf'); % hlf
        x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =opt_all(:,list.allvars=='hlg'); % hlg 
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =opt_all(:,list.allvars=='hl');  % hl
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =opt_all(:,list.allvars=='hh');  % hh
        x0(T*(find(list.opt=='hhn')-1)+1:T*(find(list.opt=='hhn'))) =opt_all(:,list.allvars=='hhn'); % hlf
        x0(T*(find(list.opt=='hln')-1)+1:T*(find(list.opt=='hln'))) =opt_all(:,list.allvars=='hln'); % hlg
    else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))) =opt_all(:,list.allvars=='hh'); % as starting value use hh
    end
    

        x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
        x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =kappaa'.*opt_all(:,list.allvars=='F'); %0.999*[Ftarget(1); Ftarget(1); Ftarget];%0.8999*0.1066/0.1159*opt_all(:,list.allvars=='F');
        x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
        
   if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))  =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
        x0(T*(find(list.opt=='wsg')-1)+1:T*(find(list.opt=='wsg')))   =opt_all(:,list.allvars=='wsg');  % An
         if indic.sep==0
                     x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))   =opt_all(:,list.allvars=='S');  % An
         end
   end  
elseif indic.target==0
    if indic.extern==0
       helper= load(sprintf('OPT_notarget_2112_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern , indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    elseif indic.extern==1
       helper= load(sprintf('OPT_notarget_0509_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul ,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    end
    opt_all=helper.opt_all;
    x0 = zeros(nn*T,1);
    if indic.noskill==0
        
        x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =opt_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =opt_all(:,list.allvars=='hhg'); % hhg
        x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =opt_all(:,list.allvars=='hlf'); % hlf
        x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =opt_all(:,list.allvars=='hlg'); % hlg 
        x0(T*(find(list.opt=='hhn')-1)+1:T*(find(list.opt=='hhn'))) =opt_all(:,list.allvars=='hhn'); % hlf
        x0(T*(find(list.opt=='hln')-1)+1:T*(find(list.opt=='hln'))) =opt_all(:,list.allvars=='hln'); % hlg 
        
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =opt_all(:,list.allvars=='hl');  % hl
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =opt_all(:,list.allvars=='hh');  % hh
    else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =opt_all(:,list.allvars=='Lf'); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =opt_all(:,list.allvars=='Lg'); % hhg
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))) =opt_all(:,list.allvars=='hh'); % as starting value use hh
    end
     x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
   
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =opt_all(:,list.allvars=='F');
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
    if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))   =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
        x0(T*(find(list.opt=='wsg')-1)+1:T*(find(list.opt=='wsg')))   =opt_all(:,list.allvars=='wsg');  % An

         if indic.sep==0
            x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))   =opt_all(:,list.allvars=='S');  % An
         end
    end 
end

  
%%
% Transform to unbounded variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- most of variables bounded by zero below
guess_trans=log(x0);

% guess_trans(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus')))=x0(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus'))) ;
%- exceptions with upper bound; hl, hh, S, and F in case of target

if indic.noskill==0
    guess_trans(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))))./...
        x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl'))));
    guess_trans(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))))./...
        x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh'))));
else
    guess_trans(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))))./...
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h'))));
end
% initial value of scientists might be zero (in version with target; other transformation required)

% if indic.sep==1
    guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))))./...
    x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))); 
    guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))))./...
    x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))); 
    guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))))./...
    x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))); 
 if indic.sep==0
%      guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=sqrt(x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
%      guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=sqrt(x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
%      guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=sqrt(x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
    guess_trans(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))=log((params(list.params=='upbarS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))))./...
    x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
 end

if indic.target==1
     guess_trans(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F'))));
%     guess_trans(T*(find(list.opt=='F')-1)+1+percon:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1+percon:T*(find(list.opt=='F'))))./...
%      x0(T*(find(list.opt=='F')-1)+1+percon:T*(find(list.opt=='F'))));
end

lb=[];
ub=[];

%%
% Test Constraints and Objective Function %%%
% percon=0;
f =  objective(guess_trans, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);
[c, ceq] = constraints_flexetaa(guess_trans, T, params, init201519, list, Ems, indic, MOM, percon, taulFixed);

% to examine stuff
% ind=1:length(ceq);
% ss=ind(abs(ceq)>1e-4);
% tt=floor(ss/T); 

%%% Optimize %%%
%%%%%%%%%%%%%%%%
%Note: Active-set algorithm is benchmark and assumed for calculations below 
%using Lagrange Multipliers. However, for convergence and accuracy, may
%first (need to) run scenario in interior-point (non-specified) and/or sqp algorithm,
%utilize results to generate initial point and then re-run active-set algorithm.

constf=@(x)constraints_flexetaa(x, T, params, init201519, list, Ems, indic, MOM, percon, taulFixed);
objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);

if indic.target==1
% if indic.notaul==1 && indic.noskill==1 && indic.xgrowth==0
%         options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% else
        options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% end
        [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
       save('incase_of_error_0201_tt', 'x')
%         [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
%        ll=load('incase_of_error_0201_3.mat', 'x')

        options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
%          save('incase_of_error_0201_4', 'x')
% ll=load('incase_of_error_0201_4', 'x');
elseif indic.target==0
%     if indic.notaul==0 %|| (indic.noskill==1 && indic.xgrowth==0)
%        options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     else
        options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     end
       [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
    %   if ~(indic.noskill==1 && indic.xgrowth==0)
       options = optimset('algorithm','active-set','TolCon',1e-11,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
        save('as_1001', "x")
        if ~ismember(exitflag, [1,4,5])
            error('active set did not solve sufficiently well')
        end
     %   end
end
% save('1101_x', "x");
% ll=load('1101_x.mat');
% x=ll.x;
%%
if indic.testT==1
    %  save(sprintf('2309_results_opt_main_notaul0_target%d', indic.target), 'x')
    % gg=  load('0308_results_opt_noskill_notaul0_notarget', 'x');
    % x=gg.x;

    %% add further direct optimization periods => goal: so that continuation value does not impact allocation anymore
    hhh= load(sprintf('2309_results_opt_main_notaul0_target%d', indic.target), 'x');
    x=hhh.x;
    T=12;
    count=0; % count number of iterations
    mm=10;
    Emsnew=Ems;
    Tinit=T;

    while mm>1e-5
        % save starting values
        Told=T;
        xold=x;
        %- transform and save for comparison
        out_trans=exp(x);
        if indic.noskill==0
            out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
            out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
        else
            out_trans((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)=upbarH./(1+exp(x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
        end

        if indic.sep==1
            out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
            out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
            out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
        elseif indic.sep==0
            out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
            out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
            out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
            out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
        end    
        if indic.target==1
            out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
        end
        out_transold=out_trans; 

        Emsnew=[Emsnew,0]; % add another net zero period to emissions limit
        Ftarget =  (Emsnew'+deltaa)/omegaa;

        % sequentially increase number of explicit periods
        count=count+1; 
        T=T+1;

        %- update initial values: add last period value as guess for new direct
        % optimization period
        x0 = zeros(nn*T,1);
            for ll=list.opt
                x0(T*(find(list.opt==ll)-1)+1:T*(find(list.opt==ll)))  = [xold(Told*(find(list.opt==ll)-1)+1:Told*(find(list.opt==ll)));xold(Told*(find(list.opt==ll)))]; 
            end
        %- optimization for new horizon
         constf=@(x)constraints_flexetaa(x, T, params, init201519, list, Emsnew, indic, MOM, percon, taulFixed);
         objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519, percon, MOM, taulFixed);

        if ~isfile(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d.mat',indic.target, count))

             options = optimset('algorithm','sqp','TolCon',1e-9,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
             [x,fval,exitflag,output,lambda] = fmincon(objf,x0,[],[],[],[],lb,ub,constf,options);
             save(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d',indic.target, count), 'x')
        else
           helper=load(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d', indic.target, count));
          x=helper.x;
        end
          abbs=zeros(length(list.opt),1);
            
          %- transform
            out_trans=exp(x);
            if indic.noskill==0
                out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
                out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
            else
                out_trans((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)=upbarH./(1+exp(x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
            end

            if indic.sep==1
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
            elseif indic.sep==0
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
                out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
            end    
        if indic.target==1
            out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
        end


          for ll=1:length(list.opt)
              %- slice by variable
              old = out_transold(Told*(ll-1)+1:Told*(ll));
              new = out_trans(T*(ll-1)+1:T*(ll));
              %- compares only first 12 periods
                abbs(ll)=max(abs((new(1:Tinit)-old(1:Tinit))./old(1:Tinit))); % x= new, xold=old
          end
          mm=max(abbs);
          fprintf('max deviation in percent %d, number of finished iterations %d', mm*100, count)

    end
end
%% transform
% helper=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus))
% x=xsqp;

out_trans=exp(x);
if indic.noskill==0
    out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
    out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
else
    out_trans((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)=upbarH./(1+exp(x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
end

%if indic.sep==1
    out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
    out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
    out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
if indic.sep==0
%     out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
%     out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
%     out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
    out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
end    


if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end

%% generate auxiliary variables
    if indic.noskill==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul,taus, lambdaa,...
            wln, wlg, wlf, SWF, S, GovCon, Tls, Tlsall, PV,PVSWF, objF]...
            = OPT_aux_vars_notaus_flex_newTauf(out_trans, list, params, T, init201519, indic, MOM, taulFixed);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
    else
       [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, GovCon, Tls, Tlsall, PV,PVSWF, objF]= OPT_aux_vars_notaus_skillHom(out_trans, list, params, T, init201519, indic, MOM);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
        taus = zeros(size(pn));

        wlf=w; wlg =w; wln=w; wh=w; wl=w; hh=h; hl=h; 
        hhf=zeros(size(pn));hhg=zeros(size(pn));hhn=zeros(size(pn));hlf=zeros(size(pn));hlg=zeros(size(pn));hln=zeros(size(pn));
    end
    
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
gammasn = zeros(size(pn));
obs =[PV, PVSWF, objF]; % save measures of objective function 
%% save stuff
% if indic.sep==1
%    opt_all=eval(symms.sepallvars);
% else
% list.allvars=saved;
opt_all=eval(symms.allvars);
% end

% additional government variabls
addGov=eval(symms.addgov);

% test if opt solution is a laissez faire solution 
% function throws error if solution is not a solution to LF
helper.LF_SIM=opt_all';
indic.limit_LF=0; % for testing no constraint on tauf
test_LF_VECT(T, list,  params,symms, init201519, helper, indic);

%%
if indic.elasE==0 && indic.sigmaWorker==0  && indic.Bop==0
if indic.count_techgap==0
    if indic.target==1
        if indic.notaul~=9
                    save(sprintf('OPT_target_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
 
        elseif indic.notaul==9
               if indic.targetWhat==0
                save(sprintf('OPT_target_taulFixed_2112_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
               else
                save(sprintf('OPT_target_taulFixed_2112_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
               end
        end

    else
        if indic.extern==1
            save(sprintf('OPT_notarget_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        else
           save(sprintf('OPT_notarget_2112_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        end
    end
else
    if indic.targetWhat~=0
        error('have not yet coded different target and counterfact technology gap')
    end
    if indic.target==1

        save(sprintf('OPT_target_countec_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
       
    else
        if indic.extern==1
            save(sprintf('OPT_notarget_0509_countec_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV, indic.sizeequ,indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        else
           save(sprintf('OPT_notarget_0509_countec_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        end
    end
end
elseif indic.elasE==1 && indic.sigmaWorker==0 && indic.Bop==0

 if indic.target==1
        save(sprintf('OPT_target_2112_elas10_emnet%d_Sun%_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
 else
        if indic.extern==1
            save(sprintf('OPT_notarget_0509_elas10_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        else
           save(sprintf('OPT_notarget_0509_elas10_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
        end
 end

elseif indic.elasE==0 && indic.sigmaWorker~=0 && indic.Bop==0
    if indic.targetWhat~=0
        error('have not yet coded different target and counterfact technology gap')
    end
    if indic.target==1
                save(sprintf('OPT_target_0509_sigmaW%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.sigmaWorker, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
    else
        error('not coded')
    end

elseif indic.elasE==0 && indic.sigmaWorker==0 && indic.Bop~=0
    if indic.targetWhat~=0
        error('have not yet coded different target and counterfact technology gap')
    end
    if indic.target==1
        save(sprintf('OPT_target_0512_Bop%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Bop,indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
    else
        save(sprintf('OPT_notarget_0512_Bop%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Bop,indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs')
    end
end

end