function []=plottsSP_tidiedUp(list, T, etaa, weightext,indic, params, Ems, plotts, percon)

% this script plots results

date="Aout22";
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
    helper=load(sprintf('SP_target_0508_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat', indic.spillovers, indic.noskill, indic.sep,  indic.xgrowth, indic.PV, etaa));
    sp_t=helper.sp_all';
    helper=load(sprintf('SP_notarget_0508_spillover%d_noskill%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',  indic.spillovers, indic.noskill, indic.sep, indic.xgrowth,  indic.PV,etaa));
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

        helper=load(sprintf('OPT_notarget_0308_spillover%d_taus0_noskill%d_notaul%d_sep%d_extern0_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill,i,  indic.sep, indic.xgrowth,indic.PV, etaa));
        opt_not_notaus=helper.opt_all';
        helper=load(sprintf('OPT_target_0308_spillover%d_taus0_noskill%d_notaul%d_sep%d_xgrowth%d_PV%d_etaa%.2f.mat',indic.spillovers, indic.noskill, i, indic.sep, indic.xgrowth,indic.PV, etaa));
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
    h1= OTHERPOL{3+1}; % taul can be used
    h2= OTHERPOL{2+1}; % taul cannot be used
    COMP = comp_CEV(h1('OPT_T_NoTaus'),h2('OPT_T_NoTaus') , varlist, varlist, symms, list, params, T, indic);
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
    for nt =  1:length(OTHERPOL)
        pp = OTHERPOL{nt};
        count = nt;
        
    for i =keys(pp)

     ii=string(i);
     
     if plotts.regime_gov ==0
         bb=0;
        allvars= RES_polcomp_notaul0(ii);
        varl=varlist; 
     else
         bb=3;
         helpp=OTHERPOL{bb};
         allvars= helpp(ii);
         varl=varlist_polcomp;
     end
     allvarsnt=pp(ii); 
     
    %% 
    fprintf('plotting %s',ii );
    for lgdind=0:1
    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
           for v=1:length(plotvars)
           gcf=figure('Visible','off');

               varr=string(plotvars(v));
               main=plot(time,allvars(find(varl==varr),:), time,allvarsnt(find(varlist_polcomp==varr),:), 'LineWidth', 1.1);

               set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
               xticks(txx)
               xlim([1, time(end)])
               ax=gca;
               ax.FontSize=13;
               ytickformat('%.2f')
               xticklabels(Year10)

            if lgdind==1
                if count==1 % version with only no income tax, but integretaed transfers 
                    lgd=legend('benchmark policy', 'no income tax', 'Interpreter', 'latex');
                elseif count == 2 %
                     lgd=legend('benchmark policy', ['no redistribution,' newline 'no income tax'], 'Interpreter', 'latex');
                elseif count == 3 %
                     lgd=legend('benchmark policy', ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex');
                elseif count == 4 %
                     lgd=legend('benchmark policy', ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                elseif count == 5 %
                     lgd=legend('benchmark policy', ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');
                end
                set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
            end
    %         sgtitle('Social Planner Allocation')
            path=sprintf('figures/all_July22/comp_bb%d_notaul%d_%s_%s_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d.png', bb, count, ii, varr, indic.spillovers,indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa, lgdind);
            exportgraphics(gcf,path,'Resolution', 400)
            % saveas(gcf,path)
            close gcf
           end % variables in group
    end % variable group
    end % legend
    end
    end
end

%% comparison to BAU
if plotts.bau==1
    fprintf('plotting bau graphs')    
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
            main=plot(time,allvars(find(varlist==varr),:), time,bau(find(varlist==varr),:), 'LineWidth', 1.1);            
           set(main, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
%             title(sprintf('%s', varr), 'Interpreter', 'latex')
           if lgdind==1
               if contains(ii, 'SP')
                  lgd=legend('Social planner', 'bau', 'Interpreter', 'latex');
               else
                  lgd=legend('Ramsey planner', 'bau', 'Interpreter', 'latex');
               end
            set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
           end

        path=sprintf('figures/all_1705/%s_BAUComp%s_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.BN, indic.ineq, indic.BN_red,indic.xgrowth,  etaa, lgdind);
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
if plotts.comptarg==1
    fprintf('plotting comparison target graphs')    
    ssr= string({'SP_T', 'SP_NOT' ,'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,3]
        ii=ssr(i);
        %- read in data
        t=string(ssr(i));
        nt=string(ssr(i+1));

%         if indic.noskill==0
             allvars= RES(t);
             allvarsnot=RES(nt); 
%         else
%              allvars= RES_ns(t);
%              allvarsnot=RES_ns(nt);
%         end

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
        path=sprintf('figures/all_July22/%s_TargetComp%s_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d.png', varr, ii, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
       close gcf
        end
        end
    end
    end      
end
%% single alternative regimes with income tax
if plotts.regimes==1
    fprintf('plotting  other regimes with and without taul graphs')   
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt = [3] % 3 is version with gov consumes env revenues but income tax is available 
        
        % read in version to version without income tax
        RES_help=OTHERPOL{nt};
    for i =1:length(opt)
        io=opt(i);
        allvars =RES_help(io);

    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
                
           main=plot(time,allvars(find(varlist_polcomp==varr),:));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   

           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

        path=sprintf('figures/all_July22/%s_SingleAltPol%s_regime%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
    end
    end
    end
end


%% comparison social planner and optimal policy (with and without labour tax)
% graph incorporates with and without laissez faire allocation 
% ensure to use different variable list for policy scenario alternations
if plotts.compeff==1
    for withlff=0
    lff=RES('LF');
    fprintf('plotting comparison efficient-optimal graphs')   
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt =  1:length(OTHERPOL) % loop over policy scenarios
        RES_help=OTHERPOL{nt};
    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        if plotts.regime_gov ==0
            bb=0;
           allvars= RES(io);
           varl=varlist;
        else
            bb =3;
            helpp = OTHERPOL{bb};
            allvars = helpp(io);
            varl=varlist_polcomp;
        end
        allvarsnotaul =RES_help(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
           if withlff==1
               main=plot(time, lff(find(varlist==varr),1:T), time,allvarseff(find(varlist==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varlist_polcomp==varr),1:T));            
               set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
           else
               main=plot(time,allvarseff(find(varlist==varr),1:T), time,allvars(find(varl==varr),1:T), time,allvarsnotaul(find(varlist_polcomp==varr),1:T));            
               set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
           end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

%             if varr=="C"
%                 ylim([0.56, 0.72]);
%             elseif varr=="hh"
%                 ylim([0.46, 0.51]);
%             elseif varr=="hl"
%                 ylim([0.31, 0.335]);                
%             end
           if lgdind==1
               if withlff==1

                   if nt ==1
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  'no income tax', 'Interpreter', 'latex');
                   elseif nt ==2
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['no redistribution,' newline 'no income tax'],  'Interpreter', 'latex');
                   elseif nt ==3
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex'); 
                   elseif nt==4
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                   elseif nt==5
                        lgd=legend('laissez-faire', 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if nt ==1
                        lgd=legend( 'efficient', 'benchmark policy',  'no income tax', 'Interpreter', 'latex');
                   elseif nt ==2
                        lgd=legend( 'efficient', 'benchmark policy',  ['no redistribution,' newline 'no income tax'], 'Interpreter', 'latex');
                   elseif nt ==3
                        lgd=legend( 'efficient', 'benchmark policy',  ['no redistribution,' newline 'with income tax'], 'Interpreter', 'latex'); 
                   elseif nt==4
                        lgd=legend( 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'with income tax'], 'Interpreter', 'latex');
                   elseif nt==5
                        lgd=legend( 'efficient', 'benchmark policy',  ['lump-sum redistribution' newline 'no income tax'], 'Interpreter', 'latex');                        
                   end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
           end
        path=sprintf('figures/all_July22/%s_CompEff%s_bb%d_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_PV%d_etaa%.2f_lgd%d_lff%d.png', varr, io, bb, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, indic.PV, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
    end
end

%% comparison alternative regimes with and without income tax
if plotts.compeff_dd==1
    for withlff=0
    lff=RES('LF');
    fprintf('plotting comparison efficient- other regimes with and without taul graphs')   
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    % for withtaul=0:1
    for nt = [3,4] % 3 is version with gov consumes env revenues but income tax is available 
        
        % read in version to version without income tax
        if nt==3
            RES_dd=OTHERPOL{2}; % comparison is separate policy regime
        elseif nt ==4
            RES_dd=OTHERPOL{5}; % lump sum trans but no taul
        end
        RES_help=OTHERPOL{nt};
    for i =1:length(eff)

        ie=eff(i);
        io=opt(i);
        allvars= RES_dd(io);
        allvarswithtaul =RES_help(io);
        allvarseff=RES(ie); 

    for l =[keys(lisst), string('addgov')]
        ll=string(l);
        if ll == string('addgov')
            plotvars= list.addgov;
        else
            plotvars=lisst(ll);
        end
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
            if ismember(varr, ["tauf", "taul", "GovCon", "Tls"])
                 
               if withlff==1
                   main=plot(time, lff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, {1; 1.2; 1},  {'LineStyle'},{'--'; '--'; ':'}, {'color'}, {grrey; orrange; 'b'} )   
               else
                   main=plot( time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, { 1.2; 1},  {'LineStyle'},{ '--'; ':'}, {'color'}, { orrange; 'b'} )   
               end

           else
                
               if withlff==1
                   main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, {1; 1.2; 1.2; 1},  {'LineStyle'},{'--';'-'; '--'; ':'}, {'color'}, {grrey; 'k'; orrange; 'b'} )   
               else
                   main=plot(time,allvarseff(find(varlist==varr),:), time,allvars(find(varlist_polcomp==varr),:), time,allvarswithtaul(find(varlist_polcomp==varr),:));            
                   set(main,{'LineWidth'}, { 1.2; 1.2; 1},  {'LineStyle'},{'-'; '--'; ':'}, {'color'}, {'k'; orrange; 'b'} )   
               end
            end
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
            if ismember(varr, ["tauf","taul", "GovCon", "Tls"]) % then no efficient graph
            
               if withlff==1
                  if nt ==3
                    lgd=legend('laissez-faire', 'no income tax', 'with income tax', 'Interpreter', 'latex'); 
                  elseif nt ==4
                     lgd=legend('laissez-faire', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                   
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                   if nt ==3
                        lgd=legend( 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   elseif nt ==4
                        lgd=legend('no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
            else    % with efficient graph            
               if withlff==1
                   if nt==3
                      lgd=legend('laissez-faire', 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   elseif nt ==4
                      lgd=legend('laissez-faire', 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                   end
                  set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               else
                  if nt ==3
                        lgd=legend( 'efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  elseif nt ==4
                        lgd=legend('efficient', 'no income tax',  'with income tax', 'Interpreter', 'latex'); 
                  end
                 set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
               end
            end
           end
        path=sprintf('figures/all_July22/%s_DDCompEff%s_pol%d_spillover%d_noskill%d_sep%d_xgrowth%d_etaa%.2f_lgd%d_lff%d.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep, indic.xgrowth, etaa, lgdind, withlff);
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
%     lff=RES('LF');
    fprintf('plotting comparison efficient')   
    eff= string({'SP_T', 'SP_NOT'});   
    % for withtaul=0:1
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
           
           main=plot( time,allvarseff(find(varlist==varr),:));            
           set(main,{'LineWidth'}, {1.2},  {'LineStyle'},{'-'}, {'color'}, {'k'} )   
           xticks(txx)
           xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
              lgd=legend('efficient', 'Interpreter', 'latex');
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_July22/%s_CompEff%s_onlyeff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d.png', varr, ie, indic.spillovers, indic.noskill, indic.sep,indic.xgrowth, indic.count_techgap, etaa, lgdind);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
end
%% comparison social planner and optimal policy benchmark
if plotts.compeff3==1
    %- only efficient and benchmark
    fprintf('plotting comparison efficient-optimal graphs')   
    for withlff=0
         lff=RES('LF');
    
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});

    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        
        if plotts.regime_gov==0
            nt=0;
            allvars= RES(io);
            varl=varlist;
        else
            nt=3;
            helpp=OTHERPOL{nt};
            allvars=helpp(io);
            varl=varlist_polcomp;
        end
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvars(find(varl==varr),:));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),:), time,allvars(find(varl==varr),:));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
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
        path=sprintf('figures/all_July22/%s_CompEff%s_regime%d_opteff_spillover%d_noskill%d_sep%d_xgrowth%d_countec%d_etaa%.2f_lgd%d_lff%d.png', varr, io, nt, indic.spillovers, indic.noskill, indic.sep,indic.xgrowth, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end

%% comparison social planner and optimal policy (with and without labour tax)
if plotts.compeff2==1
    %- only efficient and no income tax
    fprintf('plotting comparison efficient-optimal graphs')   
    for withlff=0
         lff=RES('LF');
    
    if indic.zero==0
        eff= string({'SP_T', 'SP_NOT'});
        opt=string({'OPT_T_NoTaus', 'OPT_NOT_NoTaus'});
    else
        eff= string({'SP_NOT'});
         opt=string({'OPT_NOT_NoTaus'});
    end
    % for withtaul=0:1
    for i =[1,2]

        ie=eff(i);
        io=opt(i);
        allvars= RES(io);
        allvarsnotaul =RES_polcomp_notaul(io);
        allvarseff=RES(ie); 

    for l =keys(lisst) % loop over variable groups
        ll=string(l);
        plotvars=lisst(ll);
        for lgdind=0:1
        for v=1:length(plotvars)
            gcf=figure('Visible','off');
            varr=string(plotvars(v));
      if withlff==1
            main=plot(time, lff(find(varlist==varr),:), time,allvarseff(find(varlist==varr),:), time,allvarsnotaul(find(varlist==varr),:));            
           set(main, {'LineWidth'}, {1; 1.2; 1.2}, {'LineStyle'},{'--';'-'; '--'}, {'color'}, {grrey; 'k'; orrange} )   
      else
            main=plot(time,allvarseff(find(varlist==varr),:), time,allvarsnotaul(find(varlist==varr),:));            
           set(main, {'LineWidth'}, {1.2; 1.2}, {'LineStyle'},{'-'; '--'}, {'color'}, {'k'; orrange} )   
      end
      xticks(txx)
      xlim([1, time(end)])

            ax=gca;
            ax.FontSize=13;
            ytickformat('%.2f')
            xticklabels(Year10)
            if startsWith(varr,"gA")
               xlim([1, time(end-1)])
            else
               xlim([1, time(end)])
            end

            if varr=="C"
                ylim([0.56, 0.72]);
            elseif varr=="hh"
                ylim([0.46, 0.51]);
            elseif varr=="hl"
                ylim([0.31, 0.335]);                
            end
           if lgdind==1
               if withlff==1
                    lgd=legend('laissez-faire', 'efficient', ' no income tax', 'Interpreter', 'latex');
               else
                   if varr =="tauf"
                      lgd=legend( 'social cost of emissions', 'no income tax', 'Interpreter', 'latex');
                   else
                     lgd=legend( 'efficient', ' no income tax', 'Interpreter', 'latex');
                   end
               end
              set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 19,'Orientation', 'vertical');
           end
        path=sprintf('figures/all_1705/%s_CompEff%s_noopt_spillover%d_noskill%d_sep%d_BN%d_ineq%d_red%d_xgrowth%d_zero%d_countec%d_etaa%.2f_lgd%d_lff%d.png', varr, io, indic.spillovers, indic.noskill, indic.sep,indic.BN, indic.ineq, indic.BN_red,indic.xgrowth, indic.zero, indic.count_techgap, etaa, lgdind, withlff);
        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
        end
        end
      end
    end
    end
end

end

       