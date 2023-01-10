function [symms, list, sp_all]=SP_solve_extT(list, symms, params, count, init201519, indic, T, Ems, MOM, percon)

% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C Ch Cl hh hl F sff sg sn Lg Ln h real
if indic.noskill==0
    symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn];
else
    symms.sp = [Lg Ln xn xf xg An Af Ag C h F sff sg sn];
end

Ftarget =  (Ems+deltaa)/omegaa;

if indic.xgrowth==1
    if indic.noskill==0
        symms.sp= [hhf hhg hhn hlf hlg hln xn xf xg C hh hl F];
    else
        symms.sp= [Lg Ln xn xf xg C h F];
    end        
end
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%
%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

%%
if indic.target==0
    x0 = zeros(nn*T,1);
        helper=load(sprintf('SP_notarget_2112_spillover%d_knspil%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_sizeequ0_etaa%.2f.mat',...
            indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep,indic.xgrowth, indic.PV, params(list.params=='etaa')));

%         helper=load(sprintf('OPT_notarget_plus%d_0501_emnet%d_Sun%d_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',...
%              count,indic.targetWhat, indic.Sun, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, 4,indic.sep, indic.extern, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
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
    
elseif indic.target==1 
    
    x0 = zeros(nn*T,1);
       fprintf('using sp solution as initial value')
       helper=load(sprintf('SP_target_2112_emnet%d_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', indic.targetWhat, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'obs');      %helper=load(sprintf('OPT_target_plus%d_0509_spillover%d_knspil%d_taus%d_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_sizeequ%d_GOV%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.taus, indic.noskill, 4, indic.sep, indic.xgrowth,indic.PV, indic.sizeequ, indic.GOV, params(list.params=='etaa')));
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

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT, percon);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic, percon, MOM);

if indic.target==1 && indic.noskill==0 && indic.xgrowth==0 && indic.sizeequ~=1
    fprintf('skipping sqp')
    options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);

else
    options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
    options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
    [x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);
end
save('sp_target_nn0_nsk', 'x')
%- add further periods
    Emsnew=Ems;
    Tinit=T;

for cc=1:count        % number of additional periods
    fprintf('number of iteration %d', cc)
    indic
        Told=T;
        xold=x;

        Emsnew=[Emsnew,0]; % add another net zero period to emissions limit
        Ftarget =  (Emsnew'+deltaa)/omegaa;

        % new number of explicit optimization periods
        T=T+1;

        if ~isfile((sprintf('501_results_SP_main_notaul%d_target%d_Tplus%d_nsk%d_xgr%d_knspil%d.mat',indic.notaul, indic.target, cc, indic.noskill, indic.xgrowth, indic.noknow_spill)))

        %- update initial values: add last period value as guess for new direct
        % optimization period
        x0 = zeros(nn*T,1);
            for ll=list.sp
                x0(T*(find(list.sp==ll)-1)+1:T*(find(list.sp==ll)))  = [xold(Told*(find(list.sp==ll)-1)+1:Told*(find(list.sp==ll)));xold(Told*(find(list.sp==ll)))]; 
            end
        %- optimization for new horizon
       objfSP=@(x)objectiveSP(x,T,params, list, Ftarget', indic, initOPT, percon);
       constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Emsnew, indic, percon, MOM);

         options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         [x,fval,exitflag,output,lambda] = fmincon(objfSP,x0,[],[],[],[],lb,ub,constfSP,options);
         options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         [x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);
          save(sprintf('501_results_SP_main_notaul%d_target%d_Tplus%d_nsk%d_xgr%d_knspil%d',indic.notaul, indic.target, cc, indic.noskill, indic.xgrowth, indic.noknow_spill), 'x')
        else
         hhelper=load(sprintf('501_results_SP_main_notaul%d_target%d_Tplus%d_nsk%d_xgr%d_knspil%d',indic.notaul, indic.target, cc, indic.noskill, indic.xgrowth, indic.noknow_spill), 'x');
         x=hhelper.x;
        end
end
%% transform
% helper=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus))
% x=xsqp;
T=Tinit+count;
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
    out_trans((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)=Ftarget./(1+exp(x((find(list.sp=='F')-1)*T+1+percon:find(list.sp=='F')*T)));
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
sp_all_all=eval(symms.allvars);
sp_all=sp_all_all(1:Tinit,:);

%%
if indic.count_techgap==0
if indic.target==1
    save(sprintf('SP_target_plus%d_2112_emnet%d_pillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
        count, indic.targetWhat, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all', 'sp_all_all','obs')
fprintf('saved')
else
    if indic.extern==1       
        save(sprintf('SP_notarget_plus%d_2609_spillover%d_knspil%d_noskill%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, weightext, indic.xgrowth,  indic.PV,indic.sizeequ, params(list.params=='etaa')), 'sp_all','sp_all_all', 'obs')
    else        
        save(sprintf('SP_notarget_plus%d_2112_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',...
            count, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV,indic.sizeequ,  params(list.params=='etaa')), 'sp_all','sp_all_all', 'obs')
    end
end
elseif indic.count_techgap==1
    if indic.target==1
        save(sprintf('SP_target_plus%d_2609_countec_spillover%d_knspil%d_noskill%d_sep%d_xgrowth%d_zero%d_PV%d_sizeequ%d_etaa%.2f.mat',count, indic.spillovers, indic.noknow_spill, indic.noskill, indic.sep, indic.xgrowth, indic.PV, indic.sizeequ, params(list.params=='etaa')), 'sp_all','sp_all_all', 'obs')
    else
    if indic.extern==1       
        save(sprintf('SP_notarget_plus%d_2609_countec_spillover%d_knspil%d_noskill%d_sep%d_extern%d_weightext%.2f_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat', count, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, weightext, indic.xgrowth,  indic.PV, indic.sizeequ,params(list.params=='etaa')), 'sp_all', 'sp_all_all', 'obs')
    else        
        save(sprintf('SP_notarget_plus%d_2609_countec_spillover%d_knspil%d_noskill%d_sep%d_extern%d_xgrowth%d_PV%d_sizeequ%d_etaa%.2f.mat',count, indic.spillovers,indic.noknow_spill, indic.noskill, indic.sep, indic.extern, indic.xgrowth, indic.PV, indic.sizeequ, params(list.params=='etaa')), 'sp_all','sp_all_all', 'obs')
    end
    end
end
end
