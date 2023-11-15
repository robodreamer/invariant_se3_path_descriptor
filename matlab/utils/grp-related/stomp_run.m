function [ opt_solution ] = stomp_run( algorithm_param, cost, init_pose, final_pose, traj_len, init_traj)
% run STOMP algorithm for given intial trajectories
% @originally written by Sungjoon Choi 
% @modified by Andy Park (andypark.purdue@gmail.com)

%% initialization
%= initialize the variables
global quiet;
max_iter = algorithm_param.max_iter;
joint_dim = size(init_pose,1);
sample_size = algorithm_param.sample_size;

% this length is used to ensure the initial and final positions remain
% fixed
traj_len_mid = traj_len-2;

% set the randomness (1e-4)
alpha = algorithm_param.alpha;

% set the sensitivity
h = algorithm_param.h;

tic;

%% Get GRP for init traj
eps = 1e-2;
t_anchor = [0,eps,1-eps,1]';
% x_anchor = [zeros(1,joint_dim); zeros(1,joint_dim); zeros(1,joint_dim); zeros(1,joint_dim)];
x_anchor = [init_pose'; init_pose'; final_pose'; final_pose'];
n_test = traj_len;
t_test = linspace(0,1,n_test)';
l_anchor = [1,1,1,1]';
l_test = ones(n_test,1);
kfun_str = 'kernel_levrq';
hyp_mu = [5.0, 1/4, 1.0]; % [gain/length/alpha]
sig2w_mu = 1e-4;
hyp_var = [10.0, 1/4, 10.0]; % [gain/length/alpha]
sig2w_var = 1e-4;

if ~exist('init_traj','var')
    % if init trajectory is not provided, use GRP
    x_anchor = [init_pose'; init_pose'; final_pose'; final_pose'];
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

% print out the initial cost
fprintf('STOMP -- initial cost: %.2f \n', init_cost);

if quiet == 0
    fprintf('[1/%d]init_cost is :%.2f\n',max_iter,init_cost);
end

%% Precomputation
%- generate a finite difference matrix A
% A*theta produces acceleration (theta is a trajectory vector for 1 dim)
% With R = A'*A, theta'*R*theta is the sum of squared accelerations along
% the trajectory
% A is finite difference matrix

% first fill A with central difference method
A = (diag(ones(1,traj_len_mid-1),1) + diag(-2*ones(1,traj_len_mid)) + diag(ones(1,traj_len_mid-1),-1));

Rinv = alpha*inv(A'*A); %1e-4

% Rinv2 = padarray(Rinv,[1 1]); % padded with zeros for plotting purposes
stdR = Rinv^(0.5);

%- generate M matrix used to smooth the cost gradient estimate
M = Rinv; 

%= scaling -- with each column scaled such that the maximum element is 1/traj_len
for c_idx = 1:traj_len_mid
    M(:,c_idx) = M(:,c_idx)/max(abs(M(:,c_idx)))/traj_len;
end

%% Optimization Loop

%= initialize the variables updated in the loop

% curr_traj = init_traj;
curr_traj = init_traj(:,2:end-1); % exclude the first and last positions
curr_cost = init_cost;

% learning curve stores the cost in each iteration
learning_curve = zeros(1,max_iter);
running_time = zeros(1,max_iter);

learning_curve(1) = curr_cost;
running_time(1) = toc;

opt_traj = curr_traj;
opt_cost = curr_cost;

opt_traj_data = cell(max_iter,1);
opt_traj_data{1} = init_traj;

%= loop
% sample size indicate K in the algorithm
for iter = 2:max_iter
    sample_traj = cell(1,sample_size);
    sample_cost = zeros(sample_size,traj_len_mid);
    
    %= Generate noise to the current trajectory with N(0, Rinv)
    delta_traj_data = cell(1,sample_size);
    for i = 1:sample_size        
        sample_traj{i} = zeros(joint_dim,traj_len_mid);
        for k = 1:joint_dim
            delta_traj_tmp = randn(1,traj_len_mid)*stdR';
            sample_traj{i}(k,:) = curr_traj(k,:) + delta_traj_tmp;
            delta_traj_data{i}(k,:) = delta_traj_tmp;
        end
        % compute cost
        [~,sample_cost(i,:)] = cost(sample_traj{i}); % this takes second most time
    end    
    
    %= compute probabilities of samples and normalize them
    probs_tmp = zeros(size(sample_cost));     
    probs = zeros(size(sample_cost));     
    for t = 1:traj_len_mid
        Smax = max(sample_cost(:,t));
        Smin = min(sample_cost(:,t));
        if (Smax>Smin)
            Sbar = (sample_cost(:,t) - Smin)/(Smax-Smin);
        else
            Sbar = 0;
        end
        probs_tmp(:,t) = exp(-h*Sbar);
        probs(:,t) = probs_tmp(:,t);
        
        % It works better without the following line at the expense of losing scale invariance
        probs(:,t) = probs(:,t)/sum(probs(:,t)); 
    end
    
    %= plotting delta traj
    if(0)
        figure(1); clf;
        for i = 1:sample_size
            delta_traj_data{i} = [0; delta_traj_data{i}'; 0];
            hold on; plot(delta_traj_data{i});
        end
    end

    %= merge delta trajectories based on the probabilities
    delta_traj = zeros(joint_dim,traj_len_mid);    
    for i = 1:sample_size
        for t = 1:traj_len_mid
            delta_traj(:,t) = delta_traj(:,t) + probs(i,t)*(sample_traj{i}(:,t)-curr_traj(:,t));
        end
    end
        
    %= project delta trajectory onto the basis of M matrix (M*delta_traj)
    delta_traj_proj = delta_traj;    
    for k = 1:joint_dim
        delta_traj_proj(k,:) = delta_traj(k,:)*M';
    end
    
    %= update the current trajectories with delta trajectories
    curr_traj = curr_traj + delta_traj_proj;
    curr_cost = cost(curr_traj);
    
    %= update optimal trajectory if the cost gets improved
    if opt_cost > curr_cost
        opt_cost = curr_cost;
        opt_traj = curr_traj;
        
        if quiet == 0
            fprintf('[%d/%d]opt_cost is updated!:%.2f\n',iter,max_iter,opt_cost);
        end
    else
        %         fprintf('[%d/%d]curr_cost is updated!:%.2f\n',iter,max_iter,curr_cost);
    end
    
    % update current trajectory
%     curr_traj_data{iter} = curr_traj;
    opt_traj_full = [init_pose, opt_traj, final_pose];
%     curr_traj_full = [init_traj(:,1), opt_traj, init_traj(:,end)];
    opt_traj_data{iter} = opt_traj_full;    
    
    learning_curve(iter) = opt_cost;
    running_time(iter) = toc;
end

% store the results
opt_solution.init_traj = init_traj;
opt_solution.init_cost = init_cost;
opt_solution.opt_traj = opt_traj_full;
opt_solution.opt_cost = opt_cost;
opt_solution.traj_data = opt_traj_data;
opt_solution.learning_curve = learning_curve;
opt_solution.running_time = running_time;

% print out the final cost
fprintf('STOMP -- final cost: %.2f \n', opt_cost);

end

function out = sample_init_traj_grp(grp, sample_size, cost)
%
% Sample an initial path from the GRP distribution.
%

sample_path_list = cell(sample_size,1);
sample_cost_list = cell(sample_size,1);
min_cost = 1e5;
for i = 1:sample_size % for each samples
    
    % sample a path
    sample_path_tmp = zeros(grp.n_test,grp.xdim);
    for d_idx = 1:grp.xdim % for each dimension
        sample_path_tmp(:,d_idx) = mvnrnd(zeros(grp.n_test,1),grp.K,1)';
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