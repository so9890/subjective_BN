function [symms, list, sp_all]=SP_solve(list, symms, params, Sparams, x0LF, init201014, init201519, indexx, indic, T, Ems)

% pars
read_in_params;

% function to find social planner allocation 
syms hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C Ch Cl hh hl F sff sg sn Lg Ln h real
if indic.noneutral==0
    if indic.ineq==0
    if indic.noskill==0
        symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag C hh hl F sff sg sn];
    else
        symms.sp = [Lg Ln xn xf xg An Af Ag C h F sff sg sn];
    end
    else
    if indic.noskill==0
        symms.sp = [hhf hhg hhn hlf hlg hln xn xf xg An Af Ag Ch Cl hh hl F sff sg sn];
    else
        symms.sp = [Lg Ln xn xf xg An Af Ag Ch Cl h F sff sg sn];
    end
    end

    if indic.xgrowth==1
        if indic.noskill==0
            symms.sp= [hhf hhg hhn hlf hlg hln xn xf xg C hh hl F];
        else
            symms.sp= [Lg Ln xn xf xg C h F];
        end        
    end
elseif indic.noneutral==1
    if indic.noskill==0
        if indic.xgrowth==0
            symms.sp=[hhf hhg hlf hlg xf xg Af Ag C hh hl F sff sg];

        elseif indic.xgrowth==1
           symms.sp=[hhf hhg hlf hlg xf xg C hh hl F];
        end
    elseif indic.noskill==1
        
        if indic.xgrowth==0
            symms.sp=[Lg xf xg Af Ag C h F sff sg];

        elseif indic.xgrowth==1
           symms.sp=[Lg xf xg C h F];
        end
    end
end
list.sp  = string(symms.sp); 
nn= length(list.sp); 

%%
%%%%%%%%%%%%%%%%%%%%%%
% Initial Guess %
%%%%%%%%%%%%%%%%%%%%%

%- use competitive equilibrium with policy (taus=0; tauf=0; taul=0)
if etaa ~=1 
    helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
    LF_SIM=helper.LF_SIM';
else
    helper= load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,indic.BN_red, params(list.params=='etaa')));
    LF_SIM=helper.LF_BAU';
end

if indic.xgrowth==1
    helper=load(sprintf('BAU_xgrowth_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
    LF_SIM=helper.LF_SIM;
end

if indic.count_techgap==1
    helper=load(sprintf('LF_countec_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
    LF_SIM=helper.LF_SIM;
end

if indic.noneutral==1 && indic.noskill==1
       helper=load(sprintf('LF_nonneutral_techgap%d_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_xgrowth%d_etaa%.2f.mat',indic.count_techgap, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, params(list.params=='etaa')));
       if indic.xgrowth==0
          LF_SIM=helper.LF_SIM';
       else
          LF_SIM=helper.LF_SIM;
       end
end
  %% - new list if uses separate market
  if indic.ineq==0
  if indic.sep==1
        list.allvars=list.sepallvars;
  end
  else
        if indic.sep==1
            list.allvars=list.sepallvars_ineq;
        else
            list.allvars=list.allvars_ineq;
        end
  end
%%
if indic.target==0
    x0 = zeros(nn*T,1);
    Ftarget = 0; % placeholder

    if ~isfile(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern0_xgrowth0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')))
        if indic.noneutral==0
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
        x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
        x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
        x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An

        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',1:T);  % C
        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =LF_SIM(list.allvars=='F',1:T);  % C
        x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =LF_SIM(list.allvars=='sn',1:T);  % C
        x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % C
        x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =LF_SIM(list.allvars=='sff',1:T);  % C
       
        elseif indic.noneutral==1
            if indic.noskill==0
                x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
                x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
                x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
                x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
                x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
                x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
            else
                x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg'))) =LF_SIM(list.allvars=='Lg',1:T); % hlg 
                x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
            end

            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =LF_SIM(list.allvars=='xf',1:T); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =LF_SIM(list.allvars=='xg',1:T); % hlg 
        
            x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =LF_SIM(list.allvars=='C',1:T);  % C
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =LF_SIM(list.allvars=='F',1:T);  % C
            
            if indic.xgrowth==0
                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag

                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =LF_SIM(list.allvars=='sff',1:T);  % C
            end
       end    
    else
       fprintf('using sp solution') 
        helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern0_xgrowth0_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
        sp_all=helper.sp_all;
        
        if indic.noneutral==0
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
            
            if indic.ineq==0
                x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            else
                x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch')))     =sp_all(1:T, list.allvars=='Ch');  % C
                x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl')))     =sp_all(1:T, list.allvars=='Cl');  % C
            end                
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            if indic.xgrowth==0
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C

                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
                x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
            end
       elseif indic.noneutral==1
            if indic.noskill==0
                x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
                x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
                x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
                x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
                x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
                x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            else
                x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
                x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))     =sp_all(1:T, list.allvars=='hh');  % hh
            end
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            
            if indic.ineq==0
                x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            else
                x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch')))     =sp_all(1:T, list.allvars=='Ch');  % C
                x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl')))     =sp_all(1:T, list.allvars=='Cl');  % C
            end                
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            if indic.xgrowth==0
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
            end
        end
    end
elseif indic.target==1 
    
    Ftarget = (Ems+deltaa)/omegaa;
    x0 = zeros(nn*T,1);
    kappaa = [Ftarget(1),Ftarget(1),Ftarget]./LF_SIM(list.allvars=='F',1:T); % ratio of targeted F to non-emission
    kappaa = kappaa*(1-1e-10);
    
    if ~isfile(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq,indic.BN_red,indic.xgrowth, params(list.params=='etaa')))
         fprintf('using LF solution as initial value')
        if indic.noskill==0
                       
        x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =kappaa.*LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
        x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =kappaa.*LF_SIM(list.allvars=='hhg',1:T); % hhg
        x0(T*(find(list.sp=='hhn')-1)+1:T*(find(list.sp=='hhn'))) =kappaa.*LF_SIM(list.allvars=='hhn',1:T); % hhg
        x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =kappaa.*LF_SIM(list.allvars=='hlf',1:T); % hlf
        x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =kappaa.*LF_SIM(list.allvars=='hlg',1:T); % hlg 
        x0(T*(find(list.sp=='hln')-1)+1:T*(find(list.sp=='hln'))) =kappaa.*LF_SIM(list.allvars=='hln',1:T); % hlg
        
        x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =kappaa.*LF_SIM(list.allvars=='hl',1:T);  % hl
        x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =kappaa.*LF_SIM(list.allvars=='hh',1:T);  % hh
        else
            x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg'))) =kappaa.*LF_SIM(list.allvars=='Lg',1:T); % hlg 
            x0(T*(find(list.sp=='Ln')-1)+1:T*(find(list.sp=='Ln'))) =kappaa.*LF_SIM(list.allvars=='Ln',1:T); % hlg
            x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))   =kappaa.*LF_SIM(list.allvars=='hh',1:T);  % hl
        end
        
        x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =kappaa.*LF_SIM(list.allvars=='xf',1:T); % hlf
        x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =kappaa.*LF_SIM(list.allvars=='xg',1:T); % hlg 
        x0(T*(find(list.sp=='xn')-1)+1:T*(find(list.sp=='xn')))   =kappaa.*LF_SIM(list.allvars=='xn',1:T); % hlg 
        x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =LF_SIM(list.allvars=='Af',1:T);  % Af
        x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =LF_SIM(list.allvars=='Ag',1:T);  % Ag
        x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =LF_SIM(list.allvars=='An',1:T);  % An
        x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =kappaa.*LF_SIM(list.allvars=='C',1:T);  % C
        x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =kappaa.*LF_SIM(list.allvars=='F',1:T);  % C
   else
        fprintf('using sp solution as initial value')
        helper= load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq,indic.BN_red,indic.xgrowth, params(list.params=='etaa')));
%        helper= load(sprintf('OPT_target_active_set_1905_spillover%d_noskill%d_notaul0_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')));
%        helper= load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_gammac.mat', indic.spillovers, indic.noskill, indic.sep));

        sp_all=helper.sp_all;
        if indic.noneutral==0
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
            
            if indic.ineq==0
                x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            else
                x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl')))     =sp_all(1:T, list.allvars=='Cl');  % C
                x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch')))     =sp_all(1:T, list.allvars=='Ch');  % C    
            end
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            if indic.xgrowth==0
                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
                x0(T*(find(list.sp=='An')-1)+1:T*(find(list.sp=='An')))   =sp_all(1:T, list.allvars=='An');  % An
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sn')-1)+1:T*(find(list.sp=='sn')))   =sp_all(1:T, list.allvars=='sn');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
            end
        elseif indic.noneutral==1
            if indic.noskill==0
                x0(T*(find(list.sp=='hhf')-1)+1:T*(find(list.sp=='hhf'))) =sp_all(1:T, list.allvars=='hhf'); % hhf; first period in LF is baseline
                x0(T*(find(list.sp=='hhg')-1)+1:T*(find(list.sp=='hhg'))) =sp_all(1:T, list.allvars=='hhg'); % hhg
                x0(T*(find(list.sp=='hlf')-1)+1:T*(find(list.sp=='hlf'))) =sp_all(1:T, list.allvars=='hlf'); % hlf
                x0(T*(find(list.sp=='hlg')-1)+1:T*(find(list.sp=='hlg'))) =sp_all(1:T, list.allvars=='hlg'); % hlg 
                x0(T*(find(list.sp=='hl')-1)+1:T*(find(list.sp=='hl')))   =sp_all(1:T, list.allvars=='hl');  % hl
                x0(T*(find(list.sp=='hh')-1)+1:T*(find(list.sp=='hh')))   =sp_all(1:T, list.allvars=='hh');  % hh
            
            else
                x0(T*(find(list.sp=='Lg')-1)+1:T*(find(list.sp=='Lg')))   =sp_all(1:T, list.allvars=='Lg'); % hlg
                x0(T*(find(list.sp=='h')-1)+1:T*(find(list.sp=='h')))     =sp_all(1:T, list.allvars=='hh');  % hh
            end
        
            x0(T*(find(list.sp=='xf')-1)+1:T*(find(list.sp=='xf')))   =sp_all(1:T, list.allvars=='xf'); % hlf
            x0(T*(find(list.sp=='xg')-1)+1:T*(find(list.sp=='xg')))   =sp_all(1:T, list.allvars=='xg'); % hlg 
            
            if indic.ineq==0
                x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))     =sp_all(1:T, list.allvars=='C');  % C
            else
                x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl')))     =sp_all(1:T, list.allvars=='Cl');  % C
                x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch')))     =sp_all(1:T, list.allvars=='Ch');  % C    
            end
            x0(T*(find(list.sp=='F')-1)+1:T*(find(list.sp=='F')))     =sp_all(1:T, list.allvars=='F');  % C
            if indic.xgrowth==0
                x0(T*(find(list.sp=='Af')-1)+1:T*(find(list.sp=='Af')))   =sp_all(1:T, list.allvars=='Af');  % Af
                x0(T*(find(list.sp=='Ag')-1)+1:T*(find(list.sp=='Ag')))   =sp_all(1:T, list.allvars=='Ag');  % Ag
                x0(T*(find(list.sp=='sg')-1)+1:T*(find(list.sp=='sg')))   =sp_all(1:T, list.allvars=='sg');  % C
                x0(T*(find(list.sp=='sff')-1)+1:T*(find(list.sp=='sff'))) =sp_all(1:T, list.allvars=='sff');  % C
            end
        end
        end
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
 
if indic.BN==1
    if indic.ineq==1
         guess_trans(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl')))=log((Bl-x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl'))))./...
         x0(T*(find(list.sp=='Cl')-1)+1:T*(find(list.sp=='Cl'))));
         guess_trans(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch')))=log((Bh-x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch'))))./...
         x0(T*(find(list.sp=='Ch')-1)+1:T*(find(list.sp=='Ch'))));
    else
        guess_trans(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C')))=log((B-x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C'))))./...
     x0(T*(find(list.sp=='C')-1)+1:T*(find(list.sp=='C'))));
    end
end
if indic.target==1
    % only from the third period onwards F is contrained
    guess_trans(T*(find(list.sp=='F')-1)+1+2:T*(find(list.sp=='F')))=log((Ftarget'-x0(T*(find(list.sp=='F')-1)+1+2:T*(find(list.sp=='F'))))./...
     x0(T*(find(list.sp=='F')-1)+1+2:T*(find(list.sp=='F'))));
end
lb=[];
ub=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Constraints and Objective Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- reload parameters! to be correct!
% [params, Sparams,  polCALIB,  init201014, init201519, list, symms, Ems,  Sall, x0LF, MOM, indexx]=get_params( T, indic, lengthh);

f =  objectiveSP(guess_trans,T,params, list, Ftarget, indic, initOPT);
[c, ceq] = constraintsSP(guess_trans, T, params, initOPT, list, Ems, indic);

% ind=1:length(ceq);
% ss=ind(abs(ceq)>1e-9);
% whe=floor(ss/T);

objfSP=@(x)objectiveSP(x,T,params, list, Ftarget, indic, initOPT);
constfSP=@(x)constraintsSP(x, T, params, initOPT, list, Ems, indic);


options = optimset('algorithm','sqp', 'TolCon',1e-8, 'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,guess_trans,[],[],[],[],lb,ub,constfSP,options);
     %save('sp_results_notarget_ext_sep1_spillover0_noskill1_etaa0.79_xgrowth_noneutral.mat')
%   ss=load('sp_results_target_sep1_spillover0_etaa079.mat')
options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
[x,fval,exitflag,output,lambda] = fmincon(objfSP,x,[],[],[],[],lb,ub,constfSP,options);

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
    out_trans((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T)=Ftarget'./(1+exp(x((find(list.sp=='F')-1)*T+1+2:find(list.sp=='F')*T)));
end

if indic.BN==1
    if indic.ineq==0
        out_trans((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T)=B./(1+exp(x((find(list.sp=='C')-1)*T+1:find(list.sp=='C')*T)));
    else
        out_trans((find(list.sp=='Ch')-1)*T+1:find(list.sp=='Ch')*T)=Bh./(1+exp(x((find(list.sp=='Ch')-1)*T+1:find(list.sp=='Ch')*T)));
        out_trans((find(list.sp=='Cl')-1)*T+1:find(list.sp=='Cl')*T)=Bl./(1+exp(x((find(list.sp=='Cl')-1)*T+1:find(list.sp=='Cl')*T)));

    end    
end
%%
if indic.noskill==0
    if indic.xgrowth==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, A_lag, S, SGov, Emnet, A,muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil, gammac]= SP_aux_vars_2S(out_trans, list, params, T, initOPT, indic);
    else
        [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, sff, sg, sn, ...
            F, N, G, E, Y, C, Ch, Cl, hl, hh, S, SGov, Emnet, A, muu, muuh, muul,...
            pn, pg, pf, pee, wh, wl, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_xgrowth(out_trans, list, params, T, initOPT, indic);
    end
else
    if indic.noneutral==0
            [ xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, w, wsn, wsf, wsg, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill(out_trans, list, params, T,initOPT, indic);
    elseif indic.noneutral==1
           [ xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag, sff, sg,  ...
            F, G, E, Y, C, h, A_lag, S, SGov, Emnet, A,muu,...
             pg, pf, pee, w, wsf, wsg, tauf, taul, lambdaa,...
             wlg, wlf, SWF, PVcontUtil]= SP_aux_vars_2S_noskill_noneutral(out_trans, list, params, T, initOPT, indic);
            wln=zeros(size(xf)); xn=zeros(size(xf));  An=zeros(size(xf)); Ln=zeros(size(xf)); sn=zeros(size(xf));pn=zeros(size(xf));
            N=zeros(size(xf)); wsn=zeros(size(xf));
    end
        
        wh=w; wl=w;hhf=h; hhg=h; hhn=h; hlf=h; hlg=h; hln=h; hl=h; hh=h;
end
taus = zeros(size(pn));
ws=wsn;
gammall = zeros(size(pn));
gammalh = zeros(size(pn));
gammasg = zeros(size(pn));
gammasf = zeros(size(pn));
gammasn = zeros(size(pn));

%%
if indic.sep==0
    sp_all=eval(symms.allvars);
else
    if indic.ineq==1
        sp_all=eval(symms.sepallvars_ineq);
    else
        sp_all=eval(symms.sepallvars);
    end
end

%%
if indic.count_techgap==0
if indic.target==1
    save(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth,indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
else
    if indic.extern==1       
        save(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.extern, weightext, indic.xgrowth, indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
    else        
        save(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.extern, indic.xgrowth,indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
    end
end
elseif indic.count_techgap==1
    if indic.target==1
    save(sprintf('SP_target_active_set_1705_countec_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, indic.zero,params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
    else
    if indic.extern==1       
        save(sprintf('SP_notarget_active_set_1705_countec_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.extern, weightext, indic.xgrowth, indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
    else        
        save(sprintf('SP_notarget_active_set_1705_countec_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.extern, indic.xgrowth, indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams')
    end
    end
end
end
