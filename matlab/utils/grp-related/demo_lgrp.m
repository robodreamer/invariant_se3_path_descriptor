ccc
%% Shuffle a certain region 
ccc

% First sample a smooth trajectory 
eps = 1e-2;
hyp = [1.0, 1/2, 0.1]; % [gain/length/alpha]
sig2w = 1e-4;
n_test = 200;
grp = get_grp([1,eps,1-eps,1]',[0,0,0,0]',linspace(0,1,n_test)',ones(4,1),ones(n_test,1),...
        'kernel_levrq',hyp,sig2w,hyp,sig2w);
trajs = sample_grp(grp,1); traj = trajs{1};

% Retrain GRP
l_anchor = ones(n_test,1);
l_anchor(30:80) = 0.8;
l_anchor(150:180) = 0.8;
hyp = [10.0, 1/5, 0.1]; % [gain/length/alpha]
sig2w = 1e-8;
lgrp = get_grp(linspace(0,1,n_test)',traj,linspace(0,1,n_test)',l_anchor,ones(n_test,1),...
        'kernel_levrq',hyp,sig2w,hyp,sig2w);
n_sample = 50;    
sampled_trajs = sample_grp(lgrp,n_sample);

% Plot
fig = figure(1); set_fig_size(fig,[0.1,0.3,0.5,0.4]); 
xm = 0.1; ym = 0.15; subaxes(fig,1,1,1,xm,ym); hold on;
colors = jet(n_sample);
for i = 1:length(sampled_trajs)
    sampled_traj = sampled_trajs{i};
    plot(linspace(0,1,n_test)',sampled_traj,'-','Color',colors(i,:));
end
plot(linspace(0,1,n_test)',traj,'k-','LineWidth',2);
ylim([-1,+1]);
set(gcf,'Color','w');
