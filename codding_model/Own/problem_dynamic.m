function [Obj_ram_dynamic,  optims, vecs, list ]=problem_dynamic(y, x, list, symms, E, Obj_ram , indic, pol, P )

% function reads in dynamic models: take explicitly whole time span until
% 2050 into account when maximizing. Not yet updated for a continuation
% value. 

% periods for which to solve directly
betaa=symms.params(list.params=='betaa');
%% generate vectors for each variable; save in matrix

if indic.approach==1
    list.joint=[list.y, list.x, list.ramsey_mu, string(E)];
elseif indic.approach==2
        list.joint=[list.y, list.x, list.optim, string(E)];
end

vecs=sym('a',[P, length([list.joint])]);
    for s = [list.joint]  % loop over list entries
        vecs(:,[list.joint]==s)=sym(sprintf('%s%d',s),  [P,1]); 
    end

%% replace indexed variables in static problem
% note that script "primal_problem" reads in the sum of swf and constraints
% for one period => replace variables with dynamic ones to get sum of welfare functions
% and add up to get dynamic problem
% remains to connect variables=> use growth rates; should be linked automatically through  

% can now take derivatives if static, or do value function iteration if constant tax rate!
% if dynamic policy => look at sum over next 30 years

% vector of objective function by period, already substituted technology,
% important to variables linked across periods: which they are via A

vec_Obj_ram=sym('a',[P,1]);

if indic.approach==1
    for i=1:P
        vec_Obj_ram(i)=subs(Obj_ram, [y,x,symms.ramsey_mu, E], vecs(i,:) );
    end
elseif indic.approach ==2
    for i=1:P
        vec_Obj_ram(i)=subs(Obj_ram, [y,x, list.optim, E], vecs(i,:) );
    end
end


%% LOM for technology: replace future technology as function of today's

% vectore to be substituted for Ac
Ac_help=repmat(1+pol(list.pol=='vc'),1,P).^(0:P-1)*x(list.x=='Ac');
Ad_help=repmat(1+pol(list.pol=='vdd'),1,P).^(0:P-1)*x(list.x=='Ad');

vec_Obj_ram_LOM= subs(vec_Obj_ram, [vecs(:, list.joint=='Ac'), vecs(:, list.joint=='Ad')], [transpose(Ac_help), transpose(Ad_help)]);

%% new problem: discounted sum of period objective function
%- create discount vector
disc=repmat(betaa, 1,P);
expp=1:P;
vec_discount= disc.^expp;
Obj_ram_dynamic= vec_discount*vec_Obj_ram_LOM;
%% variables over which to optimise

optims = vecs(:,ismember(list.joint,list.optim));
