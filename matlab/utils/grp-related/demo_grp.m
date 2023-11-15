ccc
%%
ccc
init_pose = 0; final_pose = 0;
eps = 1e-2;
t_anchor = [0,eps,1-eps,1]';
x_anchor = [init_pose'; init_pose'; final_pose'; final_pose'];
n_test = 100;
t_test = linspace(0,1,n_test)';
l_anchor = [1,1,1,1]';
l_test = ones(n_test,1);
kfun_str = 'kernel_levrq';
hyp_mu = [10.0, 1/2, 1.0]; % [gain/length/alpha]
sig2w_mu = 1e-4;
% Increasing gain => more deviations
% Decreasing length => more curvature
% Decreasing alpha => more curvature
hyp_vars  = {[10.0, 1/1, 5.0],[10.0, 1/2, 5.0],...
    [10.0, 1/5, 5.0],[10.0, 1/5, 0.1]};
fig = figure(1); 
set_fig_size(fig,[0.05,0.1,0.7,0.8]);
for h_idx = 1:length(hyp_vars)
    hyp_var = hyp_vars{h_idx}; % [gain/length/alpha]
    sig2w_var = 1e-4;
    
    % Get GRP
    grp = get_grp(t_anchor,x_anchor,t_test,l_anchor,l_test,...
        kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var);
    title_str = sprintf('Gain:[%.2f] Length:[%.2f] Alpha:[%.2f]',hyp_var(1),hyp_var(2),hyp_var(3));
    
    % Plot
    n_sample = 100;
    xm = 0.05; ym = 0.1;
    subaxes(fig,2,2,h_idx,xm,ym); hold on;
    colors = hsv(grp.xdim);
    
    % GRP varaince with shaded araes
    for d_idx = 1:grp.xdim
        h_fill = fill([grp.t_test;flipud(grp.t_test)],...
            [grp.mu(:,d_idx)-2*grp.std;flipud(grp.mu(:,d_idx))+2*grp.std],...
            colors(d_idx,:),'LineStyle','none'); % grp 2std
        set(h_fill,'FaceAlpha',0.2);
    end
    % GRP sampled paths
    sampled_path_list = sample_grp(grp,n_sample);
    for i = 1:n_sample
        for d_idx = 1:grp.xdim
            plot(grp.t_test,sampled_path_list{i}(:,d_idx),...
                '-','Color',0.8*colors(d_idx,:),'linewidth',0.5);
        end
    end
    % Data points (circles when #data<100, otherwise line)
    for d_idx = 1:grp.xdim
        if grp.n_anchor > 100
            % If #data is bigger than 100, do not plot
            plot(grp.t_anchor,grp.x_anchor(:,d_idx),'-',...
                'MarkerSize',12,'linewidth',2,'Color','k'); % plot data
        else
            plot(grp.t_anchor,grp.x_anchor(:,d_idx),'ko',...
                'MarkerSize',12,'linewidth',2,'MarkerFaceColor',colors(d_idx,:)); % plot data
        end
    end
    
    % GRP mu
    for d_idx = 1:grp.xdim
        hs(d_idx) = plot(grp.t_test,grp.mu(:,d_idx),'-',...
            'linewidth',3,'color',colors(d_idx,:)); % plot GRP mu
        strs{d_idx} = sprintf('%d',d_idx);
    end
    set(gcf,'color','w');
    % Title
    if title_str
        title(title_str,'fontsize',14,'interpreter','none');
    end
    grid on;
    ylim([-hyp_var(1),hyp_var(1)]);
    drawnow;
end
