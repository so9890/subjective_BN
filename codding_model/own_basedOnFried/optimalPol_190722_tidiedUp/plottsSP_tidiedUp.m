function []=plottsSP_tidiedUp(list, T, etaa, weightext,indic, params, Ems, plotts, percon)

% this script plots results

date="10Aout22";
if ~isfile(sprintf('figures/all_%s', date ))
    mkdir(sprintf('figures/all_%s', date));
end

%- color for graphs
orrange= [0.8500 0.3250 0.0980];
grrey = [0.6 0.6 0.6];

%- variables
syms hh hl Y F E N Emnet G pg pn pf pee tauf taul taus wh wl ws wsg wsn wsf lambdaa C Lg Lf Ln xn xg xf sn sff sg SWF Af Ag An A S real
syms analyTaul PV CEVv CEVvPV CEVvDy AgAf sgsff GFF EY CY hhhl whwl LgLf gAg gAf gAn gAagg Utilcon Utillab Utilsci real
symms.plotsvarsProd =[Y N E G F];
symms.plotsvarsHH =[hh hl C SWF Emnet]; 
symms.plotsvarsRes =[sn sff sg  S Af Ag An A];  
symms.plotsvarsProdIn =[xn xg xf Ln Lg Lf];  
symms.plotsvarsPol =[taus tauf taul lambdaa];  
symms.plotsvarsAdd = [analyTaul PV AgAf sgsff GFF EY CY hhhl whwl LgLf gAagg gAg gAf gAn Utilcon Utillab Utilsci];
% already exists: symms.addgov
symms.comp=[ CEVv CEVvDy CEVvPV ]; % for comparison of policy interventions, 

if indic.sep==0
    symms.plotsvarsPri =[pg pf pee pn wh wl ws];  
else
    symms.plotsvarsPri =[pg pf pee pn wh wl wsg wsn wsf];  
end

listt.plotsvarsProd=string(symms.plotsvarsProd);
listt.plotsvarsProdIn=string(symms.plotsvarsProdIn);
listt.plotsvarsHH=string(symms.plotsvarsHH);
listt.plotsvarsRes=string(symms.plotsvarsRes);
listt.plotsvarsPol=string(symms.plotsvarsPol);
listt.plotsvarsPri=string(symms.plotsvarsPri);
listt.plotsvarsAdd=string(symms.plotsvarsAdd);
list.comp=string(symms.comp);

lisst = containers.Map({'Prod', 'ProdIn','Res', 'HH', 'Pol', 'Pri', 'Add'}, {listt.plotsvarsProd, listt.plotsvarsProdIn, ...
    listt.plotsvarsRes,listt.plotsvarsHH,listt.plotsvarsPol, listt.plotsvarsPri, listt.plotsvarsAdd});
 
%- variables which to plot in one graph plus legend
lissComp = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {string([hh hl]), string([sff sg sn S]), string([gAf gAg  gAn gAagg]), string([Lf Lg])});
legg = containers.Map({'Labour', 'Science', 'Growth', 'LabourInp'}, {["high skill", "low skill"], ["fossil", "green", "non-energy", "total"], ["fossil", "green", "non-energy", "aggregate"], ["fossil", "green"]});

%- new list 
% varlist_polcomp=[list.allvars, list.addgov, string(symms.plotsvarsAdd)];

varlist=[list.allvars, string(symms.plotsvarsAdd)];
%% read in results

%-initialise cells to store matrices: one for each combination of xgr and nsk
OTHERPOL={}; % cell to save containers of different policy versions
OTHERPOL_xgr={}; 
OTHERPOL_nsk={};
OTHERPOL_xgr_nsk={}; 

% baseline results 
for xgr=0:1
for nsk =0:1
    indic.xgrowth=xgr;
    indic.noskill=nsk;
    %- sp solution independent of policy
    helper=load(sprintf('SP_target_1008_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,  indic.xgrowth, indic.PV, etaa));
    sp_t=helper.sp_all';
    helper=load(sprintf('SP_notarget_1008_spillover%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,  indic.PV,etaa));
    sp_not=helper.sp_all';

%- other results
    for i=0:5 % loop over policy versions
        if indic.xgrowth==0
            helper=load(sprintf('BAU_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,i, etaa));
            bau=helper.LF_BAU';
            helper=load(sprintf('LF_SIM_NOTARGET_spillover%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill,indic.sep, i, etaa));
            LF=helper.LF_SIM';
        else
            helper=load(sprintf('BAU_xgrowth_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep, i,etaa));
            bau=helper.LF_SIM;
            helper=load(sprintf('LF_xgrowth_spillovers%d_noskill%d_sep%d_notaul%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,i, etaa));
            LF=helper.LF_SIM;
        end

        helper=load(sprintf('OPT_notarget_1008_spillover%d_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill,i,  indic.sep, indic.xgrowth,indic.PV, etaa));
        opt_not_notaus=helper.opt_all';
        helper=load(sprintf('OPT_target_1008_spillover%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, i, indic.sep, indic.xgrowth,indic.PV, etaa));
        opt_t_notaus=helper.opt_all';

        RES = containers.Map({'BAU','LF', 'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'},...
                                {bau,  LF, sp_t, sp_not, opt_t_notaus, opt_not_notaus});
        %- add additional variables
        if xgr==0 && nsk==0
            OTHERPOL{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
        elseif xgr==0 && nsk==1
            OTHERPOL_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
        elseif xgr==1 && nsk==0
            OTHERPOL_xgr{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
        elseif xgr==1 && nsk==1
            OTHERPOL_xgr_nsk{i+1}=add_vars(RES, list, params, indic, list.allvars, symms);
        end
    end
end
end

%% Tables
if plotts.table==1
for xgr=0:1
for nsk= 0:1
    %- choose respective containers cell
    if xgr ==0 && nsk==0
    OTHERPOLL= OTHERPOL;
    elseif xgr ==1 && nsk==0
        OTHERPOLL= OTHERPOL_xgr;
    elseif xgr ==0 && nsk==1
        OTHERPOLL= OTHERPOL_nsk;
    elseif xgr ==1 && nsk==1
        OTHERPOLL= OTHERPOL_xgr_nsk;
    end
  %- discount vector
      betaa=params(list.params=='betaa');

      % version with T=11
      disc=repmat(betaa, 1,T+1);
      expp=0:T;
      vec_discount= disc.^expp;

     %- Table
     TableSWF_PV=table(keys(RES)',zeros(length(keys(RES)),1), zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1),zeros(length(keys(RES)),1));
     TableSWF_PV.Properties.VariableNames={'Allocation','Integrated', 'Scenario 1: as 0 but no taul', 'Scenario 2: Gov consu no taul', 'Scenario 3: Gov consu with taul', 'Scenario 4: Lump Sum with taul', 'Scenario 5: Lump Sum no taul'};

    %- all results
    for i =keys(RES)
         ii=string(i);
         count=0; % to keep track of which container is used
        for ccc=OTHERPOLL
            pp=ccc{1}; % to get from cell to its content
            count=count+1;
            allvars=pp(ii);
            TableSWF_PV{TableSWF_PV.Allocation==ii,count+1}=vec_discount*allvars(find(varlist=='SWF'),:)'+indic.PV*allvars(find(varlist=='PV'),1);
        end
    end
    %%
    save(sprintf('Table_SWF_%s_sep%d_noskill%d_etaa%.2f_xgrowth%d_PV%d_extern%d.mat',date, indic.sep, nsk, etaa, xgr, indic.PV, indic.extern), 'TableSWF_PV');
end
end
end
%% table CEV
if plotts.cev==1
    %- calculate CEV for a pair of policy regimes each
    if plotts.regime_gov==0
        h1= OTHERPOL{1}; % taul can be used
        h2= OTHERPOL{2}; % taul cannot be used
    elseif plotts.regime_gov==3
        h1= OTHERPOL{4}; % taul can be used
        h2= OTHERPOL{3}; % taul cannot be used
    end
        
    COMP = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);
    COMP
end


%% Plots
%- axes
time = 1:T;
txx=1:2:T; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

%- Pick main policy version for plots
if plotts.xgr ==0 && plotts.nsk==0
    OTHERPOLL= OTHERPOL;
elseif plotts.xgr ==1 && plotts.nsk==0
    OTHERPOLL= OTHERPOL_xgr;
elseif plotts.xgr ==0 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_nsk;
elseif plotts.xgr ==1 && plotts.nsk==1
    OTHERPOLL= OTHERPOL_xgr_nsk;
end
    
%% All figures single
if plotts.single_pol==1
    
    fprintf('plotting single graphs')

    %- read in variable container of chosen regime
    RES=OTHERPOLL{plotts.regime_gov+1};
%     RES=ccc;
  

    %- loop over economy versions
    for i = ["OPT_T_NoTaus" "OPT_NOT_NoTaus"]% only plotting polcies separately
        ii=string(i);
        allvars= RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);

        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ll=="HH" && varr=="Emnet"
                main=plot(time,allvars(find(varlist==varr),1:T),time(percon+1:end),Ems(1:T), 'LineWidth', 1.1);  
                set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
                lgd=legend('net emissions' , 'net emission limit',  'Interpreter', 'latex');
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

            else
                main=plot(time,allvars(find(varlist==varr),1:T), 'LineWidth', 1.1);   
                set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   

            end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if indic.count_techgap==0
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_etaa%.2f.png',date,  ii,varr , plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr,indic.extern, etaa);
            else
                path=sprintf('figures/all_%s/Single_%s_%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_countec_etaa%.2f.png',date,  ii, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, etaa);
            end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
end
%% figures single overlayed
if plotts.singov==1
    fprintf('plotting single overlayed graphs')
    RES=OTHERPOLL{plotts.regime_gov+1};
    varl=varlist;
    for lgdind=0:1
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]
        ii=string(i);
        allvars =RES(ii);
    fprintf('plotting %s',ii );
    for l =keys(lissComp) % loop over variable groups
        ll=string(l);
        plotvars=lissComp(ll); % here plotvars is a group of variable names which are to be plotted in the same graph

        gcf=figure('Visible','off');

        if length(plotvars)==2
            main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T), 'LineWidth', 1.1);    % plot vectors!        
            set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; 'k'} ) 
    %    elseif ll=="LabourInp"
    %           main=plot(time,allvars(find(varlist==plotvars(1)),:), time,allvars(find(varlist==plotvars(2)),:),time,allvars(find(varlist==plotvars(3)),:), 'LineWidth', 1.1);    % plot vectors!        
    %           set(main, {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; 'k'; 'k'} )   
       elseif ll=="Growth"
              main=plot(time(1:T),allvars(find(varl==plotvars(1)),1:T), time(1:T),allvars(find(varl==plotvars(2)),1:T),...
              time(1:T),allvars(find(varl==plotvars(3)),1:T),time(1:T),allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
              set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   

        elseif ll=="Science"
              main=plot(time,allvars(find(varl==plotvars(1)),1:T), time,allvars(find(varl==plotvars(2)),1:T),...
                  time,allvars(find(varl==plotvars(3)),1:T), time,allvars(find(varl==plotvars(4)),1:T),'LineWidth', 1.1);    % plot vectors!        
               set(main, {'LineStyle'},{'-'; '--'; ':'; '--'}, {'color'}, {'k'; 'k'; 'k'; grrey} )   
        end
           xticks(txx)
           xlim([1, time(end)])
            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
         if lgdind==1
            pp=legg(ll);
            if length(pp)==2
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),  'Interpreter', 'latex');
            elseif length(pp)==3
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)),  'Interpreter', 'latex');

            elseif length(pp)==4
                lgd=legend(sprintf('%s',pp(1)) ,sprintf('%s',pp(2)),sprintf('%s',pp(3)), sprintf('%s',pp(4)),  'Interpreter', 'latex');
            end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
         end
        if indic.count_techgap==0
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_extern%d_PV%d_etaa%.2f_lgd%d.png',date, plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,indic.extern, indic.PV,  etaa, lgdind);
        else
            path=sprintf('figures/all_%s/SingleJointTOT_regime%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_countec_extern%d_PV%d_etaa%.2f_lgd%d.png',date,  plotts.regime_gov, ii,ll, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.extern, indic.PV, etaa, lgdind);
        end
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
    end
    end
end

 %% comparison POLICY scenarios
if plotts.notaul==1
    fprintf('plotting comparison across policies') 
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    for nt =  bb
        pp = OTHERPOLL{nt};
        count = nt-1; % subtract 1 to get indicator of notaul
        
    for i =["OPT_T_NoTaus" "OPT_NOT_NoTaus"]

     ii=string(i);
     varl=varlist; 
     %- benchmark policy
     helpp=OTHERPOLL{plotts.regime_gov+1};
     allvars= helpp(ii);
     %- comparison
     allvarsnt=pp(ii); 
     
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =keys(lisst)
        ll=string(l);
      
        plotvars=lisst(ll);
      
        for v=1:length(plotvars)
           gcf=figure('Visible','off');

               varr=string(plotvars(v));
               main=plot(time,allvars(find(varl==varr),1:T), time,allvarsnt(find(varlist==varr),1:T), 'LineWidth', 1.1);

               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
               xticks(txx)
               xlim([1, time(end)])
               ax=gca;
               ax.FontSize=13;
               ytickformat('%.2f')
               xticklabels(Year10)

            if lgdind==1
                if count==0 % integrated policy
                    lgd=legend('benchmark', 'integrated policy, with income tax', 'Interpreter', 'latex');
                elseif count == 1 % integrated without income tax
                     lgd=legend('benchmark', 'integrated policy, no income tax', 'Interpreter', 'latex');
                elseif count == 2 % notaul =2 gov=tauf*F*pf, without income tax
                     lgd=legend('benchmark', 'no redistribution, no income tax', 'Interpreter', 'latex');
                elseif count == 3 % notaul=3
                     lgd=legend('benchmark', 'no redistribution, with income tax', 'Interpreter', 'latex');
                elseif count == 4 % Tls and taul
                     lgd=legend('benchmark', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                elseif count == 5 % Tls no taul
                     lgd=legend('benchmark', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
            path=sprintf('figures/all_%s/comp_benchregime%d_notaul%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', date, plotts.regime_gov, count, ii, varr, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end
    end
end

%% comparison with and without target    
%- string to loop over 
if plotts.comptarg==1
    fprintf('plotting comparison target graphs')
    RES=OTHERPOLL{plotts.regime_gov+1};
    ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,3] % to filter benchmark: policy with target
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));

         allvars= RES(t);
         allvarsnot=RES(nt); 


    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            main=plot(time,allvars(find(varlist==varr),1:T), time,allvarsnot(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('wih emission limit', 'no emission limit', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 18,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_TargetComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png', date, varr, ii, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end

%% comparison social planner and optimal policy benchmark and comparison policy
% graph incorporates with and without laissez faire allocation 
if plotts.compeff==1
    fprintf('plotting comparison efficient-optimal and non benchmark graphs')   

    %- read in container of results
    RES=OTHERPOLL{plotts.regime_gov+1};
    bb=1:length(OTHERPOLL);
    bb=bb(bb~=plotts.regime_gov+1); % drop benchmark policy
    
    for withlff=0
        lff=RES('LF');
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        varl=varlist;

         %- benchmark policy
         allvars= RES(io);
         allvarseff=RES(ie); 

     %- comparison         
     for nt =  bb % loop over policy scenarios but benchmark
        RES_help=OTHERPOLL{nt};
        count=nt-1;
        allvarsnotaul =RES_help(io);

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot(time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varl==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
           end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1

                   if count ==0
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy', 'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy', 'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy', 'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if count ==0
                        lgd=legend( 'efficient', 'benchmark policy', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'efficient', 'benchmark policy', 'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend( 'efficient', 'benchmark policy', 'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend( 'efficient', 'benchmark policy', 'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'efficient', 'benchmark policy', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'efficient', 'benchmark policy', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
                   
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_benchregime%d_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png',date, varr, io, plotts.regime_gov, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end


%% only social planner
if plotts.compeff1==1
    fprintf('plotting efficient')   

    %- read in container of results: any fine for social planner
    RES=OTHERPOLL{1};
    
    eff= ["SP_T" "SP_NOT"];   
    for i =[1,2]

        ie=eff(i);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           
           main=plot( time,allvarseff(find(varlist==varr),1:T));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
              lgd=legend('efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_onlyeff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d.png', date, varr, ie, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr, indic.count_techgap, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
end

%% comparison social planner and non-benchmark policy
if plotts.compeff2==1
    %- only efficient and no income tax
    fprintf('plotting comparison efficient-non benchmark optimal graphs')   
    
    for withlff=0
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
   
    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        io=opt(i);
    for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        lff=RES('LF');
        allvarsnotaul =RES(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),1:T), time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),1:T), time,allvarsnotaul(find(varlist==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
          
           if lgdind==1
               if withlff==1
                   if count ==0
                        lgd=legend('laissez-faire', 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend('laissez-faire', 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend('laissez-faire', 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend('laissez-faire', 'efficient', 'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               else
                  if count ==0
                        lgd=legend( 'efficient', 'integrated policy, with income tax', 'Interpreter', 'latex');
                   elseif count ==1
                        lgd=legend( 'efficient',  'integrated policy, no income tax',  'Interpreter', 'latex');
                   elseif count ==2
                        lgd=legend( 'efficient', 'no redistribution, no income tax', 'Interpreter', 'latex'); 
                   elseif count==3
                        lgd=legend( 'efficient',  'no redistribution, with income tax', 'Interpreter', 'latex');
                   elseif count==4
                        lgd=legend( 'efficient',  'lump-sum transfers, with income tax', 'Interpreter', 'latex');
                   elseif count==5
                        lgd=legend( 'efficient',  'lump-sum transfers, no income tax', 'Interpreter', 'latex');                        
                   end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_noopt_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, count, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
    RES=OTHERPOLL{plotts.regime_gov+1};

    for withlff=0
         lff=RES('LF');
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        
        allvars= RES(io);
        varl=varlist;
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varl==varr),1:T), time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varl==varr),1:T), time,allvars(find(varl==varr),1:T));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)

           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', 'optimal policy', 'Interpreter', 'latex');
               else
%                    if varr =="tauf"
%                       lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
%                    else
                     lgd=legend( 'efficient', ' optimal policy', 'Interpreter', 'latex');
%                    end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_%s/%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d_lff%d.png', date, varr, io, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep,plotts.xgr, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end

%% comparison to BAU
if plotts.bau==1
    fprintf('plotting bau graphs') 
for nt =  1:length(OTHERPOLL) % loop over policy regimes
        count=nt-1;
        RES=OTHERPOLL{nt};
        bau=RES('BAU');

    for i ={'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'}
        ii=string(i);
        allvars= RES(ii);

    %% 
    fprintf('plotting %s',ii );
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');


            varr=string(plotvars(v));
            %subplot(floor(length(plotvars)/nn)+1,nn,v)
            main=plot(time,allvars(find(varlist==varr),1:T), time,bau(find(varlist==varr),1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'status quo', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'status quo', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_BAUComp%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, ii, count, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end      
end
end

%%
if plotts.per_BAUt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end

%%
if plotts.per_LFt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to bau eff vs opt') 

        RES=OTHERPOLL{plotts.regime_gov+1};
        bau=RES('BAU');
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            rev = bau(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-rev)./rev,time,(allvars(find(varlist==varr),1:T)-rev)./rev, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFComp_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end
%%
 if plotts.per_effopt0==1
    % plot graphs in percent relative to LF, efficient/optimal world
    % without tagret dynamic and with t, percentage change over time relative to BAU scenario
    % in t=0;
     fprintf('plotting percentage relative to eff vs opt first period') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('OPT_NOT_NoTaus');
        revalleff =RES('SP_NOT');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            revopt = revall(find(varlist==varr), 1);
             reveff = revalleff(find(varlist==varr), 1);
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-reveff)./reveff,time,(allvars(find(varlist==varr),1:T)-revopt)./revopt, 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageEffOptFirstPeriod_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

%%
 if plotts.per_optd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to opt dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        revall =RES('OPT_NOT_NoTaus');
        
    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend( 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageOptDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_baud==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('BAU');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageBAUDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
end   
 
 %%
 % expressed in reduction relative to bau per period/ dynamic
if plotts.per_LFd==1
    % plot graphs in percent relative toefficient/optimal world
    % without tagret dynamic 
     fprintf('plotting percentage relative to bau dynamic') 

        RES=OTHERPOLL{plotts.regime_gov+1};
               
        allvars= RES('OPT_T_NoTaus');
        allvarseff= RES('SP_T');
        revall =RES('LF');

    %% 
    for l = keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
        gcf=figure('Visible','off');
            varr=string(plotvars(v));
            
            main=plot(time,(allvarseff(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T),time,(allvars(find(varlist==varr),1:T)-revall(find(varlist==varr), 1:T))./revall(find(varlist==varr), 1:T), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
           if lgdind==1
              lgd=legend('efficient', 'optimal policy', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_%s/%s_PercentageLFDyn_Target_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png',date, varr, plotts.regime_gov, indic.spillovers, plotts.nsk, indic.sep, plotts.xgr,  etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    
 end   

end     