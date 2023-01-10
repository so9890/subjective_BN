function [symms, list, opt_all]= OPT_solve_sep_ExtT_direct(list, symms, params, x0LF, init201519, indexx, indic, T, Ems, MOM, percon, count)

% direct opotimization for count=30 from good starting values from without
% taul versions 
% pars
read_in_params;
Ems=[Ems,zeros(1,count)];
Ftarget =  (Ems'+deltaa)/omegaa;
T=T+count;

 taulFixed=0;
%     if isfile(sprintf('OPT_target_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')))
%     fprintf('file exists')
%      return
%     end
%     if isfile(sprintf('OPT_notarget_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')))
%         fprintf('file exists')
%          return
%     end
% symbilic variables and lists
syms hhf hhg hlf hlg hhn hln C Ch Cl F G Af Ag An hl hh S sff sn sg Lf Lg h gammasg gammasf real

if indic.taus ==1
    symms.opt = [hhf hhg hlf hlg hhn hln C F G Af Ag An hl hh s sg];    
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

    if indic.notaul==3
       helper=load(sprintf('OPT_target_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul-1, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    elseif indic.notaul == 0 || indic.notaul==4
     helper=load(sprintf('OPT_target_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            count,indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul+1, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
        else
       helper=load(sprintf('OPT_target_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    end
        opt_all=helper.opt_all_all;
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
         if indic.sep==0
                     x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S')))   =opt_all(:,list.allvars=='S');  % An
         end
   end  
elseif indic.target==0
    if indic.notaul==3
          helper=load(sprintf('OPT_notarget_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul-1,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    elseif indic.notaul == 0 || indic.notaul==4
         helper=load(sprintf('OPT_notarget_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
               count,indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul+1,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    else
        helper=load(sprintf('OPT_notarget_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
    end
    opt_all=helper.opt_all_all;
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

if indic.sep==1
    guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))))./...
    x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))); 
    guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))))./...
    x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))); 
    guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=log((params(list.params=='upbarH')-x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))))./...
    x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))); 
 elseif indic.sep==0
     guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=sqrt(x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
     guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=sqrt(x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
     guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=sqrt(x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
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
                options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         else
%                 options = optimset('algorithm','sqp','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         end
                [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
                 options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
elseif indic.target==0


        options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
       [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
       options = optimset('algorithm','active-set','TolCon',1e-9,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
        if ~ismember(exitflag, [1,4,5])
            error('active set did not solve sufficiently well')
        end
end
save('incase_of_error_Direct', 'x')
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

%% generate auxiliary variables
    if indic.noskill==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, taus, lambdaa,...
            wln, wlg, wlf, SWF, S, GovCon, Tls, Tlsall, PV,PVSWF, objF]...
            = OPT_aux_vars_notaus_flex_newTauf(out_trans, list, params, T, init201519, indic, MOM);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
    else
       [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, GovCon, Tls, PV,PVSWF, objF]= OPT_aux_vars_notaus_skillHom(out_trans, list, params, T, init201519, indic, MOM);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
      
        wlf=w; wlg =w; wln=w; wh=w; wl=w; hh=h; hl=h; 
        hhf=zeros(size(pn));hhg=zeros(size(pn));hhn=zeros(size(pn));hlf=zeros(size(pn));hlg=zeros(size(pn));hln=zeros(size(pn));
    end
    
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
gammasn = zeros(size(pn));
obs =[PV, PVSWF, objF]; % save measures of objective function %% save stuff
% if indic.sep==1
%    opt_all=eval(symms.sepallvars);
% else
% list.allvars=saved;
opt_all_all=eval(symms.allvars);
opt_all=opt_all_all(1:T-count,:);
% end

% additional government variabls
addGov=eval(symms.addgov);

% test if opt solution is a laissez faire solution 
% function throws error if solution is not a solution to LF
helper.LF_SIM=opt_all_all';
indic.limit_LF=0; % for testing no constraint on tauf
test_LF_VECT(T, list,  params,symms, init201519, helper, indic);

%%
if indic.count_techgap==0
    if indic.target==1
        save(sprintf('OPT_target_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
            count, indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
    else
        if indic.extern==1
            save(sprintf('OPT_notarget_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
        else
            save(sprintf('OPT_notarget_plus%d_0501_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
               count,indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
        end
    end
else
    if indic.target==1
    save(sprintf('OPT_target_countec_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.xgrowth,indic.PV,indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
    else
    if indic.extern==1
        save(sprintf('OPT_notarget_plus%d_0509_countec_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern,weightext, indic.xgrowth,indic.PV, indic.sizeequ,indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
    else
       save(sprintf('OPT_notarget_plus%d_0509_countec_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, indic.notaul,indic.sep, indic.extern, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')), 'opt_all', 'addGov', 'obs', 'opt_all_all')
    end
    end
end
end