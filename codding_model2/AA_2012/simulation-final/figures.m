% ***************************************
% This script draws the different figures
%****************************************

% ********** Figure 1
subplot(3,2,1)
plot([0:5:395]',qe10rho15, '-g','LineWidth',4);
hold on
plot([0:5:395]',qe3rho01, '--b','LineWidth',4);
plot([0:5:395]',qe3rho15, ':r','LineWidth',4);
hold off
xlabel('year')
ylabel('clean innovation subsidy')
legend('e=10 & r=0.015','e=3 & r=0.001','e=3 & r=0.015')
subplot(3,2,2)
plot([0:5:395]',sce10rho15, '-g','LineWidth',4);
hold on
plot([0:5:395]',sce3rho01, '--b','LineWidth',4);
plot([0:5:395]',sce3rho15, ':r','LineWidth',4);
hold off
xlabel('year')
ylabel('scientists in clean innovation')
subplot(3,2,3)
plot([0:5:395]',taue10rho15, '-g','LineWidth',4);
hold on
plot([0:5:395]',taue3rho01, '--b','LineWidth',4);
plot([0:5:395]',taue3rho15, ':r','LineWidth',4);
hold off
xlabel('year')
ylabel('input tax')
subplot(3,2,4)
plot([0:5:395]',ratioe10rho15, '-g','LineWidth',4);
hold on
plot([0:5:395]',ratioe3rho01, '--b','LineWidth',4);
plot([0:5:395]',ratioe3rho15, ':r','LineWidth',4);
hold off
xlabel('year')
ylabel('share of clean inputs in total input production')
subplot(3,2,5)
plot([0:5:395]',te10rho15, '-g','LineWidth',4);
hold on
plot([0:5:395]',te3rho01, '--b','LineWidth',4);
plot([0:5:395]',te3rho15, ':r','LineWidth',4);
hold off
xlabel('year')
ylabel('increase in temperature')




%********** Figure B1 
subplot(4,2,1)
plot([0:5:395]',qe10rho15, '-.g','LineWidth',1.5);
hold on
plot([0:5:395]',qe3rho01, '-.b','LineWidth',1.5);
plot([0:5:395]',qe3rho15, '-.r','LineWidth',1.5);
plot([0:5:395]',exqe10rho15, '-g','LineWidth',3.5);
plot([0:5:395]',exqe3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exqe3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('clean innovation subsidy')
subplot(4,2,2)
plot([0:5:395]',sce10rho15, '-.g','LineWidth',1.5);
hold on
plot([0:5:395]',sce3rho01, '-.b','LineWidth',1.5);
plot([0:5:395]',sce3rho15, '-.r','LineWidth',1.5);
plot([0:5:395]',exsce10rho15, '-g','LineWidth',3.5);
plot([0:5:395]',exsce3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exsce3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('scientists in clean innovation')
subplot(4,2,3)
plot([0:5:395]',taue10rho15, '-.g','LineWidth',1.5);
hold on
plot([0:5:395]',taue3rho01, '-.b','LineWidth',1.5);
plot([0:5:395]',taue3rho15, '-.r','LineWidth',1.5);
plot([0:5:395]',extaue10rho15, '-g','LineWidth',3.5);
plot([0:5:395]',extaue3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',extaue3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('input tax')
subplot(4,2,4)
plot([0:5:395]',ratioe10rho15, '-.g','LineWidth',1.5);
hold on
plot([0:5:395]',ratioe3rho01, '-.b','LineWidth',1.5);
plot([0:5:395]',ratioe3rho15, '-.r','LineWidth',1.5);
plot([0:5:395]',exratioe10rho15, '-g','LineWidth',3.5);
plot([0:5:395]',exratioe3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exratioe3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('share of clean inputs in total input production')
subplot(4,2,5)
plot([0:5:395]',exthetae10rho15, '-g','LineWidth',3.5);
hold on
plot([0:5:395]',exthetae3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exthetae3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('resource tax')
legend('e=10 & r=0.015','e=3 & r=0.001','e=3 & r=0.015')
subplot(4,2,6)
plot([0:5:395]',exqstocke10rho15, '-g','LineWidth',3.5);
hold on
plot([0:5:395]',exqstocke3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exqstocke3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('resource stock')
subplot(4,2,7)
plot([0:5:395]',te10rho15, '-.g','LineWidth',1.5);
hold on
plot([0:5:395]',te3rho01, '-.b','LineWidth',1.5);
plot([0:5:395]',te3rho15, '-.r','LineWidth',1.5);
plot([0:5:395]',exte10rho15, '-g','LineWidth',3.5);
plot([0:5:395]',exte3rho01, '--b','LineWidth',3.5);
plot([0:5:395]',exte3rho15, ':r','LineWidth',3.5);
hold off
xlabel('year')
ylabel('increase in temperature')
legend('e=10 & r=0.015','e=3 & r=0.001','e=3 & r=0.015')




%******* Figure low epsilon
subplot(3,2,1)
plot([0:5:395]',qe3rho01, '-.or','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',qe3rho15, ':dc','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('clean innovation subsidy')
legend('epsilon=3 & rho=0.001','epsilon=3 & rho=0.015')
subplot(3,2,2)
plot([0:5:395]',sce3rho01, '-.or','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',sce3rho15, ':dc','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('scientists in clean innovation')
subplot(3,2,3)
plot([0:5:395]',taue3rho01, '-.or','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',taue3rho15, ':dc','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('input tax')
subplot(3,2,4)
plot([0:5:395]',ratioe3rho01, '-.or','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',ratioe3rho15, ':dc','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('share of clean inputs in total input production')
subplot(3,2,5)
plot([0:5:395]',te3rho01, '-.or','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',te3rho15, ':dc','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('increase in temperature')

%******** Figure high epsilon
subplot(3,2,1)
plot([0:5:395]',qe10rho01, '-sb','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',qe10rho15, '--xg','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('clean innovation subsidy')
legend('epsilon=10 & rho=0.001','epsilon=10 & rho=0.015')
subplot(3,2,2)
plot([0:5:395]',sce10rho01, '-sb','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',sce10rho15, '--xg','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('scientists in clean innovation')
subplot(3,2,3)
plot([0:5:395]',taue10rho01, '-sb','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',taue10rho15, '--xg','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('input tax')
subplot(3,2,4)
plot([0:5:395]',ratioe10rho01, '-sb','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',ratioe10rho15, '--xg','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('share of clean inputs in total input production')
subplot(3,2,5)
plot([0:5:395]',te10rho01, '-sb','LineWidth',1.5,'MarkerSize',3);
hold on
plot([0:5:395]',te10rho15, '--xg','LineWidth',1.5,'MarkerSize',3);
hold off
xlabel('year')
ylabel('increase in temperature')




