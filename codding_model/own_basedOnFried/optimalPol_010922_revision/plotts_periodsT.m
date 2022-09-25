% script to plot results with different optimization periods
date="13Sept22";
orrange= [0.8500 0.3250 0.0980];

indic.target=1;
Tinit=12;

time = 1:Tinit;
txx=1:2:Tinit; % reducing indices
%- using start year of beginning period 
Year =transpose(year(['2020'; '2025';'2030'; '2035';'2040'; '2045';'2050'; '2055'; '2060';'2065';'2070';'2075'],'yyyy'));
Year10 =transpose(year(['2020';'2030'; '2040'; '2050';'2060';'2070'],'yyyy'));

RR=containers.Map;
Ag=zeros(1,39);An=zeros(1,39);Af=zeros(1,39);sff=zeros(1,39);sn=zeros(1,39);sg=zeros(1,39);
Emnet=zeros(1,39);tauf=zeros(1,39);taul=zeros(1,39);S=zeros(1,39);
for count=1:39
        Emsnew=[Ems,zeros(1,count-1)]; % add another net zero period to emissions limit
        Ftarget =  (Emsnew'+deltaa)/omegaa;

    helper=load(sprintf('2309_results_opt_main_notaul0_target%d_Tplus%d', indic.target, count-1));
    x=helper.x;
    T=Tinit+count-1;
    hhelp =outtrans(x, list, T, Ftarget, indic, params , symms,init201519,MOM, percon);
    RR(sprintf('%d', count-1))=hhelp;
end
save('cont_mat_solutionPeriod2309', 'RR');
%%
for c=0:3:38
    SL1 = RR(sprintf('%d', c));
    SL2 = RR(sprintf('%d', c+1));
    SL3 = RR(sprintf('%d', c+2));
    for varr=["tauf", "taul", "Emnet", 'hh']
        gcf=figure('Visible','off');

        main=plot(time, SL1(1:Tinit, list.allvars==varr),time, SL2(1:Tinit,list.allvars==varr),time, SL3(1:Tinit,list.allvars==varr), 'LineWidth', 1.1);   
        set(main, {'LineStyle'},{'-';':'; '--'}, {'color'}, {'k'; orrange; 'b'} )   
        lgd=legend('0', '1', '2',  'Interpreter', 'latex');
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        xticks(txx)
        xlim([1, time(end)])
        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
        path=sprintf('figures/all_%s/compPeriodsSOl_%d_var%s.png',date,c, varr);

        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
end
%%
% longer distance
 for c=30
    SL1 = RR(sprintf('%d', c));
    SL2 = RR(sprintf('%d', c+5));
    SL3 = RR(sprintf('%d', c+8));
    for varr=["tauf", "taul", "sff" "hh"]
        gcf=figure('Visible','off');

        main=plot(time, SL1(1:Tinit, list.allvars==varr),time, SL2(1:Tinit,list.allvars==varr),time, SL3(1:Tinit,list.allvars==varr), 'LineWidth', 1.1);   
        set(main, {'LineStyle'},{'-';':'; '--'}, {'color'}, {'k'; orrange; 'b'} )   
        lgd=legend('0', 'plus 5', 'plus 8',  'Interpreter', 'latex');
        set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
        xticks(txx)
        xlim([1, time(end)])
        ax=gca;
        ax.FontSize=13;
        ytickformat('%.2f')
        xticklabels(Year10)
        path=sprintf('figures/all_%s/compPeriodsSOl_%d_var%s.png',date,c, varr);

        exportgraphics(gcf,path,'Resolution', 400)
        close gcf
    end
end
    

  