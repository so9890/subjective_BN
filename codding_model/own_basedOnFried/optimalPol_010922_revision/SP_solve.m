function [symms, list, sp_all]=SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems, MOM, percon)

% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C Ch Cl hh hl F sff sg sn Lg Ln h real
if indic.noskill==0
    symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn];
else
    symms.sp = [Lg Ln xn xf xg An Af Ag C h F sff sg sn];
end

if indic.xgrowth==1
    if indic.noskill==0
        symms.sp= [hhf hhg hhn hlf hlg hln xn xf xg C hh hl F];
    else
        symms.sp= [Lg Ln xn xf xg C h F];
    end        
end
list.sp  = string(symms.sp); 
nn= length(list.sp); 

% to start from solution with sep==1
% saved=list.allvars;
%       hhelper= load(sprintf('params_0209_sep%d', 1));
% % 
%       list.allvars=hhelper.list.allvars;
%%
%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

%- use competitive equilibrium with policy (taus=0; tauf=0; taul=0)
    helper=load(sprintf('LF_SIM_spillover%d_knspil%d_noskill%d_xgr%d_sep%d_notaul0_GOV0_sizeequ%d_etaa%.2f.mat', indic.spillovers,0, indic.noskill, indic.xgrowth, 1, indic.sizeequ, params(list.params=='etaa')));
    LF_SIM=helper.LF_SIM';

if indic.count_techgap==1
    helper=load(sprintf('LF_countec_spillovers%d_knspil%d_noskill%d_sep%d_notaul0_etaa%.2f.mat', indic.spillovers,0, indic.noskill, 1, params(list.params=='etaa')));
    LF_SIM=helper.LF_SIM;
end


% %% - new list if uses separate market
% if indic.sep==1

%     list.allvars=list.sepallvars;
% end
%%
if indic.target==0
    x0 = zeros(nn*T,1);
    Ftarget = 0; % placeholder

    if ~isfile(sprintf('SP_notarget_1008_spillover%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers, indic.noskill, 1,indic.xgrowth, 1, params(list.params=='etaa')))
        if indic.noskill==0
            x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
            x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
            x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
            x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
            x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
            x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
            x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =LF_SIM(list.allvars=='hhn',1:T); % hhg
            x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =LF_SIM(list.allvars=='hln',1:T); % hlg 

        else
            x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg'))) =LF_SIM(list.allvars=='Lg',1:T); % hlg 
            x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln'))) =LF_SIM(list.allvars=='Ln',1:T); % hlg 
          %  x0(T*(find(list.sp=='Lf')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
            x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
        end
        
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =LF_SIM(list.allvars=='xf',1:T); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =LF_SIM(list.allvars=='xg',1:T); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =LF_SIM(list.allvars=='xn',1:T); % hlg 
        if indic.xgrowth==0
        x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
        x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
        x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An

        x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =LF_SIM(list.allvars=='sn',1:T);  % C
        x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % C
        x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =LF_SIM(list.allvars=='sff',1:T);  % C
        end
        
        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',1:T);  % C
        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =LF_SIM(list.allvars=='F',1:T);  % C
    
    else
       fprintf('using sp solution') 
       helper=load(sprintf('SP_notarget_2112_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers,0, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV,indic.sizeequ,  params(list.params=='etaa')));

%      helper=load(sprintf('SP_notarget_1008_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ0_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep,indic.xgrowth, indic.PV, params(list.params=='etaa')));
%          helper=load(sprintf('SP_notarget_1008_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, 1,indic.xgrowth, 1, params(list.params=='etaa')));
%        helper=load(sprintf('OPT_notarget_0308_spillover%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill,4, indic.sep, indic.extern, indic.xgrowth, indic.PV, etaa));

        sp_all=helper.sp_all;
            if indic.noskill==0
                x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
                x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
                x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =sp_all(1:T, list.allvars=='hhn'); % hhg
                x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
                x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
                x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =sp_all(1:T, list.allvars=='hln'); % hlg
                x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
                x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            else
                x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
                x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln')))   =sp_all(1:T, list.allvars=='Ln'); % hlg
                %x0(T*(find(list.sp=='Lf')-1)+1:T*(find(list.sp=='Lf')))   =sp_all(1:T, list.allvars=='Lf');  % hl
                x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))     =sp_all(1:T, list.allvars=='hh');  % hh
            end
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg 
            
            x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
                        
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            if indic.xgrowth==0
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C

                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
                x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            end
    end
    
elseif indic.target==1 
    
    Ftarget = (Ems+deltaa)/omegaa;
    x0 = zeros(nn*T,1);
        fprintf('using sp solution as initial value')
     %   helper=load(sprintf('OPT_target_0509_Bop%d_spillover0_knspil%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',indic.Bop, 0 ,indic.noskill, 4, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, indic.GOV, etaa));

%     helper= load(sprintf('OPT_target_plus30_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, 4, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));

%         helper= load(sprintf('SP_target_1008_bop1_sigmaa0_spillover%d_knspil0_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.PV, params(list.params=='etaa')));
       %helper=load(sprintf('OPT_target_0308_spillover0_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.noskill, 4 , indic.sep, indic.xgrowth, indic.PV, etaa));
    helper=load(sprintf('SP_target_2112_emnet%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.targetWhat, indic.spillovers, 0, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, params(list.params=='etaa')));

        sp_all=helper.sp_all;
        % with new emission target
        kappaa = [repmat(Ftarget(1),1, percon),Ftarget]./sp_all(1:T,list.allvars=='F')'; % ratio of targeted F to non-emission
        kappaa = kappaa*(1-1e-10);
        if indic.noskill==0
            x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
            x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
            x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =sp_all(1:T, list.allvars=='hhn'); % hhg
            x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
            x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
            x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =sp_all(1:T, list.allvars=='hln'); % hlg 
            x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
            x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            
        else
            x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
            x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln')))   =sp_all(1:T, list.allvars=='Ln'); % hlg
            x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))     =sp_all(1:T, list.allvars=='hh');  % hh
        end
        
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =sp_all(1:T, list.allvars=='xn'); % hlg
        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C

        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =kappaa.*sp_all(1:T, list.allvars=='F')';  % C
        if indic.xgrowth==0
            x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
            x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
            x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
            x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
            x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
        end

%     end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values for An0, Ag0, Af0 refer to 2015-2019=> first period in
% Laissez faire solution; do not take init which refers to 2010-2014!
% this version here under assumption of first best policy in initial policy
% that is: tauf=0, taul=0, taus=0, lambdaa to balage budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initOPT= init201519; % as calibrated under BAU policy
% Transform variables to unbounded vars => requires less constraints! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 guess_trans=log(x0);
 if indic.noskill==0
 guess_trans(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))))./...
     x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl'))));
 guess_trans(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))))./...
     x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh'))));
 else
      guess_trans(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))=log((params(list.params=='upbarH')-x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h'))))./...
     x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h'))));
 end
 if indic.xgrowth==0
     guess_trans(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff')))=sqrt(x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))));
     guess_trans(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))=sqrt(x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn'))));
     guess_trans(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))=sqrt(x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg'))));
     guess_trans(T*(find(list.sp=='S')-1)+1:T*(find(list.sp=='S')))=log((params(list.params=='upbarS')-x0(T*(find(list.sp=='S')-1)+1:T*(find(list.sp=='S'))))./...
     x0(T*(find(list.sp=='S')-1)+1:T*(find(list.sp=='S')))); 
 end
if indic.target==1
    % only from the third period onwards F is contrained
    
    guess_trans(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F')))=log((Ftarget'-x0(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F'))))./...
     x0(T*(find(list.sp=='F')-1)+1+percon:T*(find(list.sp=='F'))));
end
lb=[];
ub=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- reload parameters! to be correct!
% [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx]=get_params( T, indic, lengthh);

f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic, initOPT, percon);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic, percon, MOM);

fprintf('sigmaa = %f, and thetaa = %f, and v=%f', sigmaa, thetaa, ((thetaa-1)/sigmaa)/(1+(thetaa-1)/sigmaa+1/sigmaa));

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT, percon);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic, percon, MOM);

% if indic.target==1 && indic.noskill==0 && indic.xgrowth==0 && indic.sizeequ~=1
%     fprintf('skipping sqp')
%     options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%     [x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
% 
% else
    options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
    options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);
% end
  save('0512_SP_target')

%%
out_trans=exp(x);
if indic.noskill==0
    out_trans((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)=upbarH./(1+exp(x((find(list.sp=='hl')-1)*T+1:find(list.sp=='hl')*T)));
    out_trans((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)=upbarH./(1+exp(x((find(list.sp=='hh')-1)*T+1:find(list.sp=='hh')*T)));
else
    out_trans((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T)=upbarH./(1+exp(x((find(list.sp=='h')-1)*T+1:find(list.sp=='h')*T)));
end
if indic.xgrowth==0
    out_trans((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)=(x((find(list.sp=='sff')-1)*T+1:find(list.sp=='sff')*T)).^2;
    out_trans((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)=(x((find(list.sp=='sg')-1)*T+1:find(list.sp=='sg')*T)).^2;
    out_trans((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)=(x((find(list.sp=='sn')-1)*T+1:find(list.sp=='sn')*T)).^2;
end

if indic.target==1
    out_trans((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)=Ftarget'./(1+exp(x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
end

%%
if indic.noskill==0
    if indic.xgrowth==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV, PVSWF, objF ]= SP_aux_vars_2S(out_trans, list, params, T, initOPT, indic);
    else
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S_xgrowth(out_trans, list, params, T, initOPT, indic);
    end
else
            [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PV,PVSWF, objF]= SP_aux_vars_2S_noskill(out_trans, list, params, T,initOPT, indic);
        wh=w; wl=w;hhf=h; hhg=h; hhn=h; hlf=h; hlg=h; hln=h; hl=h; hh=h;
end
taus = zeros(size(pn));
ws=wsn;
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
gammasg = zeros(size(pn));
gammasf = zeros(size(pn));
gammasn = zeros(size(pn));
obs =[PV, PVSWF, objF]; % save measures of objective function 
%%
%list.allvars=saved;
% if indic.sep==0
sp_all=eval(symms.allvars);
% else
%     sp_all=eval(symms.sepallvars);
% end

%%
if indic.count_techgap==0
if indic.target==1
    save(sprintf('SP_target_2112_emnet%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.targetWhat, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs')
fprintf('saved')
else
    if indic.extern==1       
        save(sprintf('SP_notarget_1008_spillover%d_knspil%d_noskill%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, weightext, indic.xgrowth,  indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs')
    else   
        if indic.Bop==0
            save(sprintf('SP_notarget_2112_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV,indic.sizeequ,  params(list.params=='etaa')), 'sp_all', 'obs')
        elseif indic.Bop==1
            save(sprintf('SP_notarget_0512_bop%d_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.Bop, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV,indic.sizeequ,  params(list.params=='etaa')), 'sp_all', 'obs')
        end
   end
end
elseif indic.count_techgap==1
    if indic.target==1
        save(sprintf('SP_target_1008_bop%d_sigmaa%d_countec_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_zero%d_PV%d_sizeequ%d_etaa%.2f.mat',indic.Bop, indic.sigmaWorker, indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs')
    else
    if indic.extern==1       
        save(sprintf('SP_notarget_1008_countec_spillover%d_knspil%d_noskill%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, weightext, indic.xgrowth,  indic.PV, indic.sizeequ,params(list.params=='etaa')), 'sp_all',  'obs')
    else        
        save(sprintf('SP_notarget_1008_countec_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV, indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs')
    end
    end
elseif indic.sigmaWorker==1
        save(sprintf('SP_target_1008_sigmaW%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.sigmaWorker, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs')
end
end
