function [symms, list, opt_all]= OPT_solve_sep(list, symms, params, Sparams, x0LF, init201519, indexx, indic, T, Ems)

% pars
read_in_params;
Ftarget =  (Ems'+deltaa)/omegaa;
 
% symbilic variables and lists
syms hhf hhg hlf hlg hhn hln C Ch Cl F G Af Ag An hl hh S sff sn sg Lf Lg h gammasg gammasf real

if indic.taus ==1
    symms.opt = [hhf hhg hlf hlg hhn hln C F G Af Ag An hl hh s sg];    
else
    if indic.ineq==0
        if indic.noskill==0
            symms.opt = [hhf hhg hlf hlg hhn hln C F G hl hh sn sff sg];
        else
            if indic.noneutral==0
                symms.opt = [Lf Lg C F G h sn sff sg];
            else
                symms.opt = [Lf Lg C F G h sff sg];
            end
        end
    else
         if indic.noskill==0
             symms.opt = [hhf hhg hlf hlg hhn hln Ch Cl F G hl hh sn sff sg];
         else
            symms.opt = [Lf Lg Ch Cl F G h sg sff sn];
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
if indic.sep==1
    if indic.ineq==0
    list.allvars=list.sepallvars;
    symms.allvars=symms.sepallvars;
    else
     list.allvars=list.sepallvars_ineq;
    symms.allvars=symms.sepallvars_ineq;
    end
end

%% Initial Guess %%% 
%%%%%%%%%%%%%%%%%%%%%
if indic.target==1
    if isfile(sprintf('OPT_target_active_set_1905_spillover0_taus0_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_NEWems.mat',indic.noskill, 1, indic.sep, indic.BN,indic.ineq, indic.BN_red, indic.xgrowth,etaa))
     
        helper=load(sprintf('OPT_target_active_set_1905_spillover0_taus0_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_NEWems.mat',indic.noskill, 0, indic.sep, indic.BN,indic.ineq, indic.BN_red, indic.xgrowth,etaa));
        opt_all=helper.opt_all;

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
    
    if indic.ineq==0
        x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
    else
        x0(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch')))     =opt_all(:,list.allvars=='Ch');   % C
        x0(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl')))     =opt_all(:,list.allvars=='Cl');   % C
    end
        x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =0.9999*opt_all(:,list.allvars=='F'); %0.999*[Ftarget(1); Ftarget(1); Ftarget];%0.8999*0.1066/0.1159*opt_all(:,list.allvars=='F');
        x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
   if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))   =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
   end
    else

    helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth0_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red, Sparams.etaa));
    helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth,indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams');
    helper= load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.extern, weightext, indic.xgrowth, indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams');
helper =load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_etaa%.2f_nonneut%d_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth,indic.zero, params(list.params=='etaa'), indic.noneutral), 'sp_all', 'Sparams');

    sp_all=helper.sp_all;

     if ~isvarname('sp_all')
         error('did not load sp solution')
     end 

     x0 = zeros(nn*T,1);
     if indic.noskill==0
         x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =sp_all(:,list.allvars=='hhf'); % hhf; first period in LF is baseline
         x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =sp_all(:,list.allvars=='hhg'); % hhg
         x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =sp_all(:,list.allvars=='hlf'); % hlf
         x0(T*(find(list.opt=='hhn')-1)+1:T*(find(list.opt=='hhn'))) =sp_all(:,list.allvars=='hhn'); % hhg
         x0(T*(find(list.opt=='hln')-1)+1:T*(find(list.opt=='hln'))) =sp_all(:,list.allvars=='hln'); % hlf
         x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =sp_all(:,list.allvars=='hlg'); % hlg 
         x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =sp_all(:,list.allvars=='hl');  % hl
         x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =sp_all(:,list.allvars=='hh');  % hh
     else
         x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =sp_all(:,list.allvars=='Lf'); % hlf
         x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =sp_all(:,list.allvars=='Lg'); % hlg 
         x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h')))   =sp_all(:,list.allvars=='hh');  % hh    
     end
     if indic.ineq==0
        x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =sp_all(:,list.allvars=='C');   % C
     else
        x0(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch')))     =sp_all(:,list.allvars=='C');   % C
        x0(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl')))     =sp_all(:,list.allvars=='C');   % C
     end
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =0.9999*sp_all(:,list.allvars=='F');
     x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =sp_all(:,list.allvars=='G');   % G
     if indic.xgrowth==0
         x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))) =sp_all(:,list.allvars=='sff');  % Af
         x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =sp_all(:,list.allvars=='sg');  % Ag
         if indic.noneutral==0
            x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =0.999*sp_all(:,list.allvars=='sn');  % An
         else
            x0(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))  =ones(size(sp_all(:,list.allvars=='sn')));  % An
            x0(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))   =ones(size(sp_all(:,list.allvars=='sn')));  % An
         end
     end
  end  
elseif indic.target==0
       
    if isfile(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern0_xgrowth0_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq,indic.BN_red,  etaa))
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext0.01_xgrowth%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill,1,indic.sep,indic.BN, indic.ineq,indic.BN_red,indic.extern, indic.xgrowth, etaa));
   % helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_extern1_weightext0.1_xgrowth0_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill,indic.sep,indic.BN, indic.ineq,indic.BN_red, etaa));
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.extern,weightext, indic.xgrowth,indic.zero,  params(list.params=='etaa')), 'opt_all', 'Sparams')
     
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
    
    if indic.ineq==0
        x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =opt_all(:,list.allvars=='C');   % C
    else
        x0(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch')))     =opt_all(:,list.allvars=='Ch');   % C
        x0(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl')))     =opt_all(:,list.allvars=='Cl');   % C
    end
    x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =opt_all(:,list.allvars=='F');
    x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =opt_all(:,list.allvars=='G');   % G
    if indic.xgrowth==0
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))   =opt_all(:,list.allvars=='sff');  % Af
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =opt_all(:,list.allvars=='sg');  % Ag
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =opt_all(:,list.allvars=='sn');  % An
    end 
else    

    %- use competitive equilibrium with policy (taus=0; tauf=0; taul=0)
     if etaa ~=1 
        helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red, etaa));
        LF_SIM=helper.LF_SIM';
     else
        helper= load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, params(list.params=='etaa')));
        LF_SIM=helper.LF_BAU';
     end

    if indic.noneutral==1
            helper=load(sprintf('LF_nonneutral_techgap%d_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat',indic.count_techgap, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, params(list.params=='etaa')));
            LF_SIM=helper.LF_SIM';
    end
     x0 = zeros(nn*T,1);
    
     if indic.noskill==0
         x0(T*(find(list.opt=='hhf')-1)+1:T*(find(list.opt=='hhf'))) =LF_SIM(list.allvars=='hhf',1:T); % hhf; first period in LF is baseline
         x0(T*(find(list.opt=='hhg')-1)+1:T*(find(list.opt=='hhg'))) =LF_SIM(list.allvars=='hhg',1:T); % hhg
         x0(T*(find(list.opt=='hlf')-1)+1:T*(find(list.opt=='hlf'))) =LF_SIM(list.allvars=='hlf',1:T); % hlf
         x0(T*(find(list.opt=='hlg')-1)+1:T*(find(list.opt=='hlg'))) =LF_SIM(list.allvars=='hlg',1:T); % hlg 
         x0(T*(find(list.opt=='hl')-1)+1:T*(find(list.opt=='hl')))   =LF_SIM(list.allvars=='hl',1:T);  % hl
         x0(T*(find(list.opt=='hh')-1)+1:T*(find(list.opt=='hh')))   =LF_SIM(list.allvars=='hh',1:T);  % hh
         x0(T*(find(list.opt=='hhn')-1)+1:T*(find(list.opt=='hhn'))) =LF_SIM(list.allvars=='hhn',1:T); % hhg
         x0(T*(find(list.opt=='hln')-1)+1:T*(find(list.opt=='hln'))) =LF_SIM(list.allvars=='hln',1:T); % hlf

     else
        x0(T*(find(list.opt=='Lf')-1)+1:T*(find(list.opt=='Lf'))) =LF_SIM(list.allvars=='Lf',1:T); % hhf; first period in LF is baseline
        x0(T*(find(list.opt=='Lg')-1)+1:T*(find(list.opt=='Lg'))) =LF_SIM(list.allvars=='Lg',1:T); % hhg
        x0(T*(find(list.opt=='h')-1)+1:T*(find(list.opt=='h')))   =LF_SIM(list.allvars=='hh',1:T);
     end
     
     x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))     =LF_SIM(list.allvars=='C',1:T);   % C
     x0(T*(find(list.opt=='F')-1)+1:T*(find(list.opt=='F')))     =LF_SIM(list.allvars=='F',1:T);
     x0(T*(find(list.opt=='G')-1)+1:T*(find(list.opt=='G')))     =LF_SIM(list.allvars=='G',1:T);   % G
     x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))) =LF_SIM(list.allvars=='sff',1:T);  % Af
     x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))   =LF_SIM(list.allvars=='sg',1:T);  % Ag
     if indic.noneutral==0
         x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))   =LF_SIM(list.allvars=='sn',1:T);  % An
     else
         x0(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))   =LF_SIM(list.allvars=='gammasg',1:T);  % An
         x0(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))   =LF_SIM(list.allvars=='gammasf',1:T);  % An
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
if indic.target==1
     if etaa<1
        guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))))./...
        x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))); 
        guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))))./...
        x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))); 
        guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))))./...
        x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))));
     else
         guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=sqrt(x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
         guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=sqrt(x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));
         guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=sqrt(x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))));%;% sqrt(params(list.params=='upbS')-x0(T*(find(list.opt=='S')-1)+1:T*(find(list.opt=='S'))));

     end
else
    guess_trans(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff'))))./...
    x0(T*(find(list.opt=='sff')-1)+1:T*(find(list.opt=='sff')))); 
    guess_trans(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn'))))./...
    x0(T*(find(list.opt=='sn')-1)+1:T*(find(list.opt=='sn')))); 
    guess_trans(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))=log((params(list.params=='upbarH')*indic.minn-x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg'))))./...
    x0(T*(find(list.opt=='sg')-1)+1:T*(find(list.opt=='sg')))); 
end

if indic.target==1
    guess_trans(T*(find(list.opt=='F')-1)+1+2:T*(find(list.opt=='F')))=log((Ftarget-x0(T*(find(list.opt=='F')-1)+1+2:T*(find(list.opt=='F'))))./...
     x0(T*(find(list.opt=='F')-1)+1+2:T*(find(list.opt=='F'))));
end
if indic.noneutral==1
    guess_trans(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))= sqrt(x0(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf'))));    
    guess_trans(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))= sqrt(x0(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg'))));
end
if indic.BN==1
    if indic.ineq==0
         guess_trans(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C')))=log((B-x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C'))))./...
         x0(T*(find(list.opt=='C')-1)+1:T*(find(list.opt=='C'))));
    else
         guess_trans(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch')))=log((Bh-x0(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch'))))./...
         x0(T*(find(list.opt=='Ch')-1)+1:T*(find(list.opt=='Ch'))));
         guess_trans(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl')))=log((Bl-x0(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl'))))./...
         x0(T*(find(list.opt=='Cl')-1)+1:T*(find(list.opt=='Cl'))));

    end
end
lb=[];
ub=[];


%%
% Test Constraints and Objective Function %%%
f =  objective(guess_trans, T, params, list, Ftarget, indic, init201519);
[c, ceq] = constraints_flexetaa(guess_trans, T, params, init201519, list, Ems, indic);

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

    constf=@(x)constraints_flexetaa(x, T, params, init201519, list, Ems, indic);
    objf=@(x)objective(x, T, params, list, Ftarget, indic, init201519);
%if indic.target==0
%else
% options = optimset('algorithm','sqp','TolCon',1e-6,'Tolfun',1e-12,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% options = optimset('Tolfun',1e-6,'MaxFunEvals',1000000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
% THIS ONE DOES NOT WORK WELL WHEN OTHERS FIND SOLUTION:
% options = optimset('algorithm','active-set','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%end
if indic.target==1
%         options = optimset('algorithm','sqp','TolCon',1e-6,'Tolfun',1e-10,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
%         if indic.spillovers==1
%             options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         else
        options = optimset('algorithm','sqp','TolCon',1e-10,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
%         end
         [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
        options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
         [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
%           save('1905_opt_target_sep1_etaa079_emsnew_ineq0_notaul_BN_red075')
%            ll=load('1905_opt_notarget_sep1_etaa079_emsnew_notaul_BN_red1')
%         [xsqp,fval,exitflag,output,lambda] = fmincon(ss.objf,xsqp,[],[],[],[],ss.lb,ss.ub,ss.constf,options);
% 
% ss=load(sprintf('sqp_solu_notargetOPT_1005_spillover%d_taus%d_noskill%d.mat', indic.spillovers, indic.taus, indic.noskill))
%         [x,fval,exitflag,output,lambda] = fmincon(ss.objf,ss.x,[],[],[],[],ss.lb,ss.ub,ss.constf,options);
elseif indic.target==0
       options = optimset('algorithm','sqp','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);

    %    options = optimset('algorithm','active-set','TolCon',1e-6,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
       [x,fval,exitflag,output,lambda] = fmincon(objf,guess_trans,[],[],[],[],lb,ub,constf,options);
       save('opt_1905_etaa079_sep1_notarget_notaul0_BN0_ineq0_red0_noskill_nonneut')
%         options = optimset('algorithm','active-set','TolCon',1e-8,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        %[x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
        options = optimset('algorithm','active-set','TolCon',1e-11,'Tolfun',1e-6,'MaxFunEvals',500000,'MaxIter',6200,'Display','iter','MaxSQPIter',10000);
        [x,fval,exitflag,output,lambda] = fmincon(objf,x,[],[],[],[],lb,ub,constf,options);
%         ss=load(sprintf('active_set_solu_notargetOPT_505_spillover%d_taus%d_possible', indic.spillovers, indic.taus), 'x')
        if ~ismember(exitflag, [1,4,5])
            error('active set did not solve sufficiently well')
        end
end
%f =  objective(x, T, params, list, Ftarget, indic);

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

if indic.target == 0
    out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
    out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
    out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));

else
    if etaa<1
    out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
    out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
    out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH*indic.minn./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
    else
        out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
    out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
    out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
    end 
end
if indic.noneutral==1
    out_trans(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))= (x(T*(find(list.opt=='gammasf')-1)+1:T*(find(list.opt=='gammasf')))).^2;    
    out_trans(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))= (x(T*(find(list.opt=='gammasg')-1)+1:T*(find(list.opt=='gammasg')))).^2;
end
if indic.target==1
    out_trans((find(list.opt=='F')-1)*T+1+2:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+2:find(list.opt=='F')*T)));
end
if indic.BN==1
    if indic.ineq==0
        out_trans((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T)=B./(1+exp(x((find(list.opt=='C')-1)*T+1:find(list.opt=='C')*T)));
    else
       out_trans((find(list.opt=='Ch')-1)*T+1:find(list.opt=='Ch')*T)=Bh./(1+exp(x((find(list.opt=='Ch')-1)*T+1:find(list.opt=='Ch')*T)));
       out_trans((find(list.opt=='Cl')-1)*T+1:find(list.opt=='Cl')*T)=Bl./(1+exp(x((find(list.opt=='Cl')-1)*T+1:find(list.opt=='Cl')*T)));
    end
end
%% save results

    if indic.noskill==0
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S, gammac]= OPT_aux_vars_notaus_flex(out_trans, list, params, T, init201519, indic);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
    else
         if indic.noneutral==0
       [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom(out_trans, list, params, T, init201519, indic);
        gammasg = zeros(size(pn));
        gammasf = zeros(size(pn));
         elseif indic.noneutral==1
        [xf,xg,Ag, Af,...
            Lg, Lf, Af_lag, Ag_lag,sff, sg,  ...
            F, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pg, pf, pee,  ws, wsf, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, gammasg, gammasf]= OPT_aux_vars_notaus_skillHom_nn(out_trans, list, params, T, init201519, indic);
              wln=zeros(size(xf)); xn=zeros(size(xf));  An=zeros(size(xf)); Ln=zeros(size(xf)); sn=zeros(size(xf));pn=zeros(size(xf));
                N=zeros(size(xf)); wsn=zeros(size(xf));
         end
        wlf=w; wlg =w; wln=w; wh=w; wl=w; hh=h; hl=h; 
            hhf=zeros(size(pn));hhg=zeros(size(pn));hhn=zeros(size(pn));hlf=zeros(size(pn));hlg=zeros(size(pn));hln=zeros(size(pn));
    end
taus = zeros(size(pn));
gammall = zeros(size(pn));
gammalh = zeros(size(pn));

gammasn = zeros(size(pn));
%save('1006_target_taul_noskill_xgrowth_nonneutral')
%%
if indic.sep==1
    if indic.ineq==1
        opt_all=eval(symms.sepallvars_ineq);
    else
        opt_all=eval(symms.sepallvars);
    end
else
    opt_all=eval(symms.allvars);
end

% test if opt solution is a laissez faire solution 
% function throws error if solution is not a solution to LF
helper.LF_SIM=opt_all';
test_LF_VECT(T, list,  params,symms, init201519, helper, indic);
%%
if indic.count_techgap==0
    if indic.target==1
        save(sprintf('OPT_target_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_elasEN%d_etaa%.2f_NEWems.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth,indic.zero, indic.subs, params(list.params=='etaa')), 'opt_all', 'Sparams')
    else
        if indic.extern==1
            save(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_elasEN%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.extern,weightext, indic.xgrowth,indic.zero, indic.subs, params(list.params=='etaa')), 'opt_all', 'Sparams')
        else
           save(sprintf('OPT_notarget_active_set_1905_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_xgrowth%d_zero%d_elasEN%d_etaa%.2f_nonneut%d.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.extern, indic.xgrowth, indic.zero, indic.subs, params(list.params=='etaa'), indic.noneutral), 'opt_all', 'Sparams')
        end
    end
else
    if indic.target==1
    save(sprintf('OPT_target_active_set_1905_countec_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_elasEN%d_etaa%.2f_NEWems.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul, indic.sep, indic.BN, indic.ineq, indic.BN_red, indic.xgrowth, indic.zero, indic.subs, params(list.params=='etaa')), 'opt_all', 'Sparams')
    else
    if indic.extern==1
        save(sprintf('OPT_notarget_active_set_1905_countec_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_weightext%.2f_xgrowth%d_zero%d_elasEN%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.extern,weightext, indic.xgrowth, indic.zero, indic.subs, params(list.params=='etaa')), 'opt_all', 'Sparams')
    else
       save(sprintf('OPT_notarget_active_set_1905_countec_spillover%d_taus%d_noskill%d_notaul%d_sep%d_BN%d_ineq%d_red%d_extern%d_xgrowth%d_zero%d_elasEN%d_etaa%.2f.mat', indic.spillovers, indic.taus, indic.noskill, indic.notaul,indic.sep,indic.BN, indic.ineq, indic.BN_red, indic.extern, indic.xgrowth,indic.zero, indic.subs, params(list.params=='etaa')), 'opt_all', 'Sparams')
    end
    end
end
end