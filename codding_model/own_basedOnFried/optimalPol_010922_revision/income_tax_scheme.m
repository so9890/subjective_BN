%% Income tax function %% 

% function to investigate functioning of income tax scheme
lambdaa=[0.7, 1.2, 2]';
taul = 0.181;
pre_y=linspace(0,10);
post_y=lambdaa.*pre_y.^(1-taul);
T=pre_y-lambdaa.*pre_y.^(1-taul); % transfers = gov revenues
lin_post_y= lambdaa.*pre_y+T;

% linear versus non-linear tax schedule
plot(pre_y, post_y(2,:), pre_y, lin_post_y(2,:), pre_y, pre_y, 'LineWidth', 1.1)
lgd=legend('post-tax income non-linear', 'post-tax incomelinear','pre-tax income',  'Interpreter', 'latex');
set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
       
% effect of lambda
plot(pre_y, post_y, pre_y, pre_y, 'LineWidth', 1.1)
lgd=legend('post-tax income non-linear 0.7','post-tax income non-linear 1.2','post-tax income non-linear 2', 'pre-tax income',  'Interpreter', 'latex');
set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');

% effect taul 
taull =[-0.1, 0.1, 0.4]';
postt=2*pre_y.^(1-taull);
plot(pre_y, postt, 'LineWidth', 1.1)
lgd=legend('regressive', 'taul=0.1','taul=0.4',  'Interpreter', 'latex');
set(lgd, 'Interpreter', 'latex', 'Location', 'best', 'Box', 'off','FontSize', 20,'Orientation', 'vertical');
