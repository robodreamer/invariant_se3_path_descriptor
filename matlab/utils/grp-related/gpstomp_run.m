function [ opt_solution ] = gpstomp_run( algorithm_param, cost, waypts, traj_len, init_traj)
% run GP-STOMP algorithm 
% this algorithm uses GRP as basis function to sample random trajectories. 
% updated in July 2019

%% initialization
%= initialize the variables
global quiet;
max_iter = algorithm_param.max_iter;
% Ndim = numel(init_pose);
sample_size = algorithm_param.sample_size;

% set the stepsize
alpha = algorithm_param.alpha;

% set the sensitivity
h = algorithm_param.h;

% update grp mean with optimal delta traj at each iteration
update_grp_mean_traj = true;

%% Get GRP for init traj
eps = 1e-2;
% t_anchor = [0,eps,1-eps,1]';
% x_anchor = [init_pose'; init_pose'; final_pose'; final_pose'];

%- automatically determine anchor points and time from waypoints
% waypoints should be Npts x Ndim and Npts > 2 
Npts = size(waypts,1);
Ndim = size(waypts,2);
Nanchors = Npts + 2;

%- get time anchors
t_anchor = zeros(Nanchors, 1);
t_anchor(1) = 0;
if (Npts > 2)
    t_anchor(2:end-1) = linspace(eps, 1-eps, Npts);
elseif (Npts == 2)
    t_anchor(2) = eps;
    t_anchor(end-1) = 1-eps;
end
t_anchor(end) = 1;

%- get position anchors
x_anchor = zeros(Nanchors, Ndim);
x_anchor(1,:) = waypts(1,:);
if (Npts > 2)
    x_anchor(2:end-1,:) = waypts;
elseif (Npts == 2)
    x_anchor(2,:) = waypts(1,:);
    x_anchor(end-1,:) = waypts(end,:);
end
x_anchor(end,:) = waypts(end,:);

n_test = traj_len;
t_test = linspace(0,1,n_test)';
% l_anchor = [1,1,1,1]';
l_anchor = ones(Nanchors,1);
l_test = ones(n_test,1);
kfun_str = 'kernel_levrq';
hyp_mu = [5.0, 1/4, 1.0]; % [gain/length/alpha]
sig2w_mu = 1e-4;
% hyp_var = [10.0, 1/4, 10.0]; % [gain/length/alpha] % default
hyp_var = [10.0, 1/10, 10.0]; % [gain/length/alpha] % experiment
% hyp_var = [10.0, 0.8, 10.0]; % [gain/length/alpha] % simplified shape

% varying length depending on Npts
varyLengthWithNpts = true;
if(varyLengthWithNpts)
    hyp_var(2) = 1/Npts;
end

% here increasing length simplifies the shapes
sig2w_var = 1e-4;

if ~exist('init_traj','var')
    % if init trajectory is not provided, use GRP
%     x_anchor = [init_pose'; init_pose'; final_pose'; final_pose'];
    grp1 = get_grp(t_anchor,x_anchor,t_test,l_anchor,l_test,...
    kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var);
%     fprintf('GRP for init traj ready.\n');
    
    %= Sample an initial trajectory
    out = sample_init_traj_grp(grp1, algorithm_param.init_sample_size, cost);
    init_traj = out.init_traj;
    init_cost = out.init_cost;
else % when init traj is provided
    init_cost = cost(init_traj);
end

% initialize the state variables
curr_traj = init_traj;
curr_cost = init_cost;

% print out the initial cost
fprintf('GP-STOMP -- initial cost: %.2f \n', init_cost);

tic;

if (quiet == 0)
    fprintf('[1/%d]init_cost is :%.2f\n', max_iter, curr_cost);
end

%% Get GRP for random noise trajectory sampling
% x_anchor_zero = [zeros(1,joint_dim); zeros(1,joint_dim); zeros(1,joint_dim); zeros(1,joint_dim)];
x_anchor_zero = zeros(Nanchors, Ndim);
grp2 = get_grp(t_anchor,x_anchor_zero,t_test,l_anchor,l_test,...
    kfun_str,hyp_mu,sig2w_mu,hyp_var,sig2w_var);
% fprintf('GRP for random sample traj ready.\n');

%% Loop
%= initialize the variables updated in the loop

% learning curve stores the cost in each iteration
learning_curve = zeros(1,max_iter);
running_time = zeros(1,max_iter);

learning_curve(1) = curr_cost;
running_time(1) = toc;

opt_traj = curr_traj;
opt_cost = curr_cost;

opt_traj_data = cell(max_iter,1);
opt_traj_data{1} = init_traj;

if (init_cost < 1e-5) % if the init cost is zero, then skip optimization
    
    for iter = 2:max_iter
        opt_traj_data{iter} = opt_traj;
        learning_curve(iter) = opt_cost;
        running_time(iter) = toc;
    end
    
else
    
    for iter = 2:max_iter
        
        %= Sample paths
        delta_traj_data = cell(sample_size,1);
        sample_cost = zeros(1, sample_size);
        chol_K = chol(grp2.K, 'lower');
        for i = 1:sample_size % for each samples
            sampled_path = zeros(grp2.n_test,grp2.xdim);
            for d_idx = 1:grp2.xdim % for each dimension
                sampled_path(:,d_idx) = chol_K*randn(grp2.n_test,1);
%                 sampled_path(:,d_idx) = mvnrnd(zeros(grp2.n_test,1),grp2.K,1)';
            end
            delta_traj_tmp = sampled_path + grp2.mu;
            delta_traj_data{i} = delta_traj_tmp';
            delta_traj_data{i} = alpha*delta_traj_data{i};
            sample_traj = curr_traj + delta_traj_data{i};
            % compute cost
            sample_cost(i) = cost(sample_traj); % this takes second most time
        end
        
        %= compute probabilities of samples and normalize them
        probs_tmp = zeros(1,sample_size);
        probs = probs_tmp;
        Smax = max(sample_cost);
        Smin = min(sample_cost);
        Sbar = (sample_cost - Smin)/(Smax-Smin);
        probs_tmp = exp(-h*Sbar);
        probs = probs_tmp;        
        probs = probs/sum(probs);
        
        %= plotting delta traj
        if(0)
            figure(1); clf;
            for i = 1:sample_size
                hold on; plot(delta_traj_data{i});
            end
            disp('Press Enter');
            pause();
        end
        
        %= merge delta trajectories based on the probabilities
        delta_traj = zeros(Ndim,traj_len);
        for i = 1:sample_size
            delta_traj = delta_traj + probs(i)*(delta_traj_data{i});
        end
        
        %= plotting the merged delta traj
        if(0)
            figure(1); clf;
            plot(delta_traj)
        end
        
        %= update the current trajectories with delta trajectories
        curr_traj = curr_traj + delta_traj;
        curr_cost = cost(curr_traj);
        
        %= update optimal trajectory if the cost gets improved
        if opt_cost > curr_cost
            opt_cost = curr_cost;
            opt_traj = curr_traj;
            
            % update grp mu with the new delta_traj
            if (update_grp_mean_traj)
                grp2.mu = delta_traj';
            end
            
            if quiet == 0
                fprintf('[%d/%d]opt_cost is updated!:%.2f\n',iter,max_iter,opt_cost);
            end
        else
            %         fprintf('[%d/%d]curr_cost is updated!:%.2f\n',iter,max_iter,curr_cost);
        end
        
        % store the information from the current iteration
        opt_traj_data{iter} = opt_traj;
        learning_curve(iter) = opt_cost;
        running_time(iter) = toc;
        
        % if the current cost is zero, terminate the optimization
        if (curr_cost < 1e-5)
            for iterTmp = iter:max_iter
                opt_traj_data{iterTmp} = opt_traj;
                learning_curve(iterTmp) = opt_cost;
                running_time(iterTmp) = toc;
            end
            iter = iterTmp;
        end
    end
end

%% store result
opt_solution.init_traj = init_traj;
opt_solution.init_cost = init_cost;
opt_solution.opt_traj = opt_traj;
opt_solution.opt_cost = opt_cost;
opt_solution.traj_data = opt_traj_data;
opt_solution.learning_curve = learning_curve;
opt_solution.running_time = running_time;

% print out the final cost
fprintf('GP-STOMP -- final cost: %.2f \n', opt_cost);

end

function out = sample_init_traj_grp(grp, sample_size, cost)
%
% Sample an initial path from the GRP distribution.
%

sample_path_list = cell(sample_size,1);
sample_cost_list = cell(sample_size,1);
min_cost = 1e5;
chol_K = chol(grp.K, 'lower');
for i = 1:sample_size % for each samples
    
    % sample a path
    sample_path_tmp = zeros(grp.n_test,grp.xdim);
    for d_idx = 1:grp.xdim % for each dimension
        sample_path_tmp(:,d_idx) = chol_K*randn(grp.n_test,1);
%         sample_path_tmp(:,d_idx) = mvnrnd(zeros(grp.n_test,1),grp.K,1)';
    end
    sample_path_tmp2 = sample_path_tmp + grp.mu;    
    sample_path_list{i} = sample_path_tmp2;
    
    % get current cost
    curr_cost_tmp = cost(sample_path_tmp2');
    sample_cost_list{i} = curr_cost_tmp;
    
    % store the sample traj with the lowest cost
    if curr_cost_tmp < min_cost
        min_cost = curr_cost_tmp;
        min_cost_traj = sample_path_tmp2;
    end
end

%= plotting init traj
if(0)
    figure(1); clf;
    plot(min_cost_traj);
end

% store the results
out = struct();
out.init_cost = min_cost;
out.init_traj = min_cost_traj';
out.sample_path_list = sample_path_list;
out.sample_cost_list = sample_cost_list;
end
