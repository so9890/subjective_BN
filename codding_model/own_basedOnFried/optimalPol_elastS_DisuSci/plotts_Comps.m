function []=plotts_Comps(list, T, etaa, indic, params)

if ~isfile('figures/testfig')
    mkdir(sprintf('figures/all_1705'));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.8 0.8 0.8];
if indic.ineq==0
    if indic.sep==1
        varlist=list.sepallvars;
    else
        varlist=list.allvars;
    end
else
    if indic.sep==1
        varlist=list.sepallvars_ineq;
    else
        varlist=list.allvars_ineq;
    end
end
% this script plots results

syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws wsg wsn wsf lambdaa Ch Cl C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An real
%- additional vars
syms AgAf sgsff GFF EY CY hhhl whwl Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
if indic.ineq==0
    symms.plotsvarsHH =[hh hl C SWF Emnet]; 
else
    symms.plotsvarsHH =[hh hl Ch Cl  SWF Emnet];  
end
symms.plotsvarsRes =[sn sff sg  Af Ag An];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [AgAf sgsff GFF EY CY hhhl whwl Utilcon Utillab Utilsci ];

if indic.sep==0
    symms.plotsvarsPri =[pg pf pee pn wh wl ws];  
else
    symms.plotsvarsPri =[pg pf pee pn wh wl wsg wsn wsf];  
end

% symms.plotspol = [tauf taul taus lambdaa]; 
% symms.plotsprices = [pg pn pf pee wh wl ws];

listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);
listt.plotsvarsAdd=string(symms.plotsvarsAdd);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%% read in results

%- read in version with compariosn across utility specification: no
%  reduction :
indic.BN_red=0; % cannot reduce with log utility!

for BN=0:1 % loop over utility specifications to read in both versions

    indic.BN=BN;

    helper=load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    bau=helper.LF_BAU';
    % helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
    % fb_lf=helper.LF_SIM';
    helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,   etaa));
    sp_t=helper.sp_all';
    helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,  indic.BN_red,  etaa));
    sp_not=helper.sp_all';
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_not_notaus=helper.opt_all';
    helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_t_notaus=helper.opt_all';

    RES = containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                            {bau,  sp_t, sp_not, opt_t_notaus, opt_not_notaus});
    %- add additional variables
    [RES]=add_vars(RES, list, params, indic, varlist, symms);

    %- results without taul
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_not_notaus_notaul=helper.opt_all';
    helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_t_notaus_notaul=helper.opt_all';


    RES_polcomp_full   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus, opt_not_notaus});
    RES_polcomp_notaul = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul});
    %- add additional variables
    RES_polcomp_full   =add_vars(RES_polcomp_full, list, params, indic, varlist, symms);
    RES_polcomp_notaul =add_vars(RES_polcomp_notaul, list, params, indic, varlist, symms);

    if indic.BN==1
        RES_BN1=RES;
        RES_polcomp_full_BN1=RES_polcomp_full;
        RES_polcomp_notaul_BN1=RES_polcomp_notaul;
    else
        RES_BN0=RES;
        RES_polcomp_full_BN0=RES_polcomp_full;
        RES_polcomp_notaul_BN0=RES_polcomp_notaul;
    end
end

%- read in versions with comparison across reduction
indic.BN=1; % cannot reduce with log utility!

for red=0:1 % loop over utility specifications to read in both versions

    indic.BN_red=red;

    helper=load(sprintf('LF_BAU_spillovers%d_noskill%d_sep%d_bn%d_ineq%d_red%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    bau=helper.LF_BAU';
    % helper=load(sprintf('FB_LF_SIM_NOTARGET_spillover%d_noskill0_sep0.mat', indic.spillovers));
    % fb_lf=helper.LF_SIM';
    helper=load(sprintf('SP_target_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_EMnew.mat', indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,   etaa));
    sp_t=helper.sp_all';
    helper=load(sprintf('SP_notarget_active_set_1705_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq,  indic.BN_red,  etaa));
    sp_not=helper.sp_all';
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_not_notaus=helper.opt_all';
    helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul0_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_t_notaus=helper.opt_all';

    RES = containers.Map({'BAU', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                            {bau,  sp_t, sp_not, opt_t_notaus, opt_not_notaus});
    %- add additional variables
    [RES]=add_vars(RES, list, params, indic, varlist, symms);

    %- results without taul
    helper=load(sprintf('OPT_notarget_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_not_notaus_notaul=helper.opt_all';
    helper=load(sprintf('OPT_target_active_set_1905_spillover%d_taus0_noskill%d_notaul1_sep%d_BN%d_ineq%d_red%d_etaa%.2f_NEWems.mat',indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa));
    opt_t_notaus_notaul=helper.opt_all';


    RES_polcomp_full   = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus, opt_not_notaus});
    RES_polcomp_notaul = containers.Map({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},{ opt_t_notaus_notaul, opt_not_notaus_notaul});
    %- add additional variables
    RES_polcomp_full   =add_vars(RES_polcomp_full, list, params, indic, varlist, symms);
    RES_polcomp_notaul =add_vars(RES_polcomp_notaul, list, params, indic, varlist, symms);

    if indic.BN_red==1
        RES_RED=RES;
        RES_polcomp_full_RED=RES_polcomp_full;
        RES_polcomp_notaul_RED=RES_polcomp_notaul;
    else
        RES_NORED=RES;
        RES_polcomp_full_NORED=RES_polcomp_full;
        RES_polcomp_notaul_NORED=RES_polcomp_notaul;
    end
end

%- update variable list
varlist=[varlist, string(symms.plotsvarsAdd)];

%- x axis stuff
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));


%% All figures single: Comparison with BN and red

 for compp=1 % loop over comparison versions (reduction or utility (BN))
     if compp==0 % Utility version
        RES=RES_BN1; % main results
        RES_gr = RES_BN0;
     elseif compp==1
         RES = RES_NORED;
         RES_gr = RES_RED;
     end
for lgdind=0:1

for i =keys(RES)
    %- loop
    ii=string(i);
    allvars= RES(ii);
    allvars_gr= RES_gr(ii);
%% 
fprintf('plotting %s',ii );
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);

    for v=1:length(plotvars)
        gcf=figure('Visible','off');


        varr=string(plotvars(v));
        main=plot(time,allvars(find(varlist==varr),:),time,allvars_gr(find(varlist==varr),:), 'LineWidth', 1.1);     
        if compp==0
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; grrey} )  
        else
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )  
        end
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
        if lgdind==1
            if compp==0
                lgd=legend('satiation point', 'log utility', 'Interpreter', 'latex');
            else
               lgd=legend('no reduction', 'voluntary reduction', 'Interpreter', 'latex');
            end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        end
    if compp==0
        path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_sep%d_BNCOMP_ineq%d_red0_etaa%.2f_lgd%d.png',  ii,varr, indic.spillovers, indic.sep, indic.ineq,  etaa, lgdind);
    else
       path=sprintf('figures/all_1705/Single_%s_%s_spillover%d_sep%d_BN1_ineq%d_redCOMP_etaa%.2f_lgd%d.png',  ii,varr, indic.spillovers, indic.sep, indic.ineq,  etaa, lgdind);
    end
    exportgraphics(gcf,path,'Resolution', 400)
    close gcf
    end
end
end
end
 end

%% comparison to efficient allocation % deviation with and without reduction 

eff= string({'SP_T', 'SP_NOT'});
opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

for i =[1,2]
    
    ie=eff(i);
    io=opt(i);
    allvars= RES_NORED(io);
    allvarseff=RES_NORED(ie); 
    
    allvars_red=RES_RED(io);
    allvarseff_red=RES_RED(ie);
    
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
        gcf=figure('Visible','off');
        varr=string(plotvars(v));
        if varr~=string('SWF')
            main=plot(time,(allvars(find(varlist==varr),:)-allvarseff(find(varlist==varr),:))./allvarseff(find(varlist==varr),:)*100, time, (allvars_red(find(varlist==varr),:)-allvarseff_red(find(varlist==varr),:))./allvarseff_red(find(varlist==varr),:)*100, 'LineWidth', 1.2);            
        else
            main=plot(time,(allvars(find(varlist==varr),:)-allvarseff(find(varlist==varr),:))./abs(allvarseff(find(varlist==varr),:))*100, time, (allvars_red(find(varlist==varr),:)-allvarseff_red(find(varlist==varr),:))./abs(allvarseff_red(find(varlist==varr),:))*100, 'LineWidth', 1.2);            
        end
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
       if lgdind==1
          lgd=legend('no reduction', 'voluntary reduction', 'Interpreter', 'latex');
          set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
       end
    path=sprintf('figures/all_1705/%s_CompEff%s_spillover%d_sep%d_BN1_ineq%d_redCOMP_etaa%.2f_lgd%d.png', varr, io, indic.spillovers, indic.sep, indic.ineq, etaa, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
    close gcf
    end
    end
end
end

%% comparison with and without taul with and without reduction and efficient
%- as households reduce, can the government reach a more efficient
%  allocation?; What if reduction is unequal?

RES=RES_polcomp_full_BN1; % main results
RES_gr = RES_polcomp_full_BN0;
for i =keys(RES_polcomp_full)
    
    ii=string(i);
         allvars= RES_polcomp_full(ii);
         allvarsnt=RES_polcomp_notaul(ii); 
         TableSWF_PV.NoTaul(TableSWF_PV.Allocation==ii)=vec_discount*allvarsnt(find(varlist=='SWF'),:)';

%% 
fprintf('plotting %s',ii );
for lgdind=1
for l =keys(lisst)
    ll=string(l);
    plotvars=lisst(ll);
    % number of figures in row in subplot
%     if ll~='VARS'
%         nn=2;
%     else 
%     nn=3;
%     end
    %%% with subplots
       for v=1:length(plotvars)
       gcf=figure('Visible','off');
        
     
            varr=string(plotvars(v));
%             subplot(floor(length(plotvars)/nn)+1,nn,v)
            main=plot(time,allvars(find(varlist==varr),:), time,allvarsnt(find(varlist==varr),:), 'LineWidth', 1.1);
            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
      
        if lgdind==1
          
        lgd=legend('full model', 'no income tax', 'Interpreter', 'latex');
           
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        end
%         sgtitle('Social Planner Allocation')
        path=sprintf('figures/all_1705/comp_notaul_%s_%s_spillover%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', ii, varr, indic.spillovers, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        % saveas(gcf,path)
        close gcf
       end % variables in group
end % variable group
end % legend
end      

%% comparison to BAU
% RES = containers.Map({'BAU', 'FB_LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
%                        {bau, fb_lf, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
for j=0
    indic.noskill=j;
    if indic.noskill==0
        bau=RES('BAU');
    else
        bau=RES_ns('BAU');
    end
for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
    ii=string(i);

    if indic.noskill==0
         allvars= RES(ii);
    else
         allvars= RES_ns(ii);
    end

%% 
fprintf('plotting %s',ii );
for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
    gcf=figure('Visible','off');
        

        varr=string(plotvars(v));
        %subplot(floor(length(plotvars)/nn)+1,nn,v)
        main=plot(time,allvars(find(varlist==varr),:), time,bau(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
        title(sprintf('%s', varr), 'Interpreter', 'latex')
       if lgdind==1
           if contains(ii, 'SP')
              lgd=legend('social planner', 'bau', 'Interpreter', 'latex');
           else
              lgd=legend('ramsey planner', 'bau', 'Interpreter', 'latex');
           end
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
       end

    path=sprintf('figures/all_1705/%s_BAUComp%s_spillover%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.sep, indic.BN, indic.ineq, indic.BN_red, etaa, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
    % saveas(gcf,path)
    close gcf
%     pause
    end
    end
  end
end      
end
%% comparison with and without target    
%- string to loop over 
ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
for j=0
    indic.noskill=j;
for i =[1,3]
    ii=ssr(i);
    %- read in data
    t=string(ssr(i));
    nt=string(ssr(i+1));
    
    if indic.noskill==0
         allvars= RES(t);
         allvarsnot=RES(nt); 
    else
         allvars= RES_ns(t);
         allvarsnot=RES_ns(nt);
    end

for l =keys(lisst) % loop over variable groups
    ll=string(l);
    plotvars=lisst(ll);
    for lgdind=0:1
    for v=1:length(plotvars)
        gcf=figure('Visible','off');
        varr=string(plotvars(v));
        main=plot(time,allvars(find(varlist==varr),:), time,allvarsnot(find(varlist==varr),:), 'LineWidth', 1.1);            
       set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
       xticks(txx)
       xlim([1, time(end)])

        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
       if lgdind==1
          lgd=legend('wih emission target', 'no emission target', 'Interpreter', 'latex');
          set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
       end
    path=sprintf('figures/all_1705/%s_TargetComp%s_spillover%d_sep%d_BN%d_ineq%d_red%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.sep,indic.BN, indic.ineq, indic.BN_red,etaa, lgdind);
    exportgraphics(gcf,path,'Resolution', 400)
   close gcf
    end
    end
  end
end      
end



end

       