% Adapt the input trajectory for a given input and target pose
% this method solves an NLP to find another path in SE3 that best preserves
% the shape description represented by DHB invariant se3 path shape descriptor.
function result = generate_trajectory(inputPoseData, params, T_init, T_final)

    result = struct();

    % get the position and rotation data from the input
    pos_data_orig = inputPoseData.pos_data;
    rotm_data_orig = inputPoseData.rotm_data;
    N = size(pos_data_orig,1);

    rvec_data_orig = zeros(N,3);
    for i=1:N
        rvec_data_orig(i, :) = rotationMatrixToVector(rotm_data_orig(:,:,i));
    end

    % resample the path
    Nframes = params.Nframes;
    timeNew = linspace(1, size(pos_data_orig,1), Nframes-3);
    timeNew = [timeNew timeNew(end) timeNew(end) timeNew(end)];
    pos_data = SpatialRobotModel.cubicInterp(pos_data_orig, 1:N, timeNew);
    rvec_data = SpatialRobotModel.cubicInterp(rvec_data_orig, 1:N, timeNew);
    N_orig = Nframes;
    N = N_orig;

    % resize the data
    rotm_data = zeros(3,3,N);
    for i=1:N
        rvec = rvec_data(i,:);
        rotm = rotationVectorToMatrix(rvec);
        rotm_data(:,:,i) = rotm;
    end

    % smooth the data a bit
    if (params.smoothing)
        pos_data_orig = pos_data;
        pos_data = smoothdata(pos_data, "sgolay", 10);
        rvec_data = smoothdata(rvec_data, "sgolay", 10);
    end

    % store the resampled input data
    result.pos_data_demo = pos_data;
    result.rotm_data_demo = rotm_data;
    result.rvec_data_demo = rvec_data;

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    % Compute position based DHB invariants
    [pos_invariant, rot_invariant, p_frame_init, R_frame_init] = ...
        computeDHBMex_mex(pos_diff_data, rvec_data(1:end-1,:), 'pos', T_init);
    invariants_demo = [pos_invariant, rot_invariant];

    % update the number of output
    N = N_orig - 3;

    %% Try trajectory adaptation with OCP

    %%% construct an ocp with casadi
    opti = casadi.Opti();

    % Initial values
    init_traj = [pos_data(1:N, :), rvec_data(1:N, :)];

    % constraints
    constraints = struct();
    constraints.start_pose = [T_init(1:3,4)', ...
        rotationMatrixToVector(T_init(1:3,1:3))];
    constraints.target_pose = [T_final(1:3,4)', ...
        rotationMatrixToVector(T_final(1:3,1:3))];

    % Define system states
    p_DHB_var = cell(1,N);
    rvec_DHB_var = cell(1,N);
    U_var = opti.variable(N, 6);

    for k=1:N
        p_DHB_var{k} = opti.variable(1,3); % position
        rvec_DHB_var{k} = opti.variable(1,3); % rotation vector
    end

    % Initialize linear and angular frames
    linear_frame = p_frame_init;
    angular_frame = R_frame_init;

    % Constraints on the start pose
    opti.subject_to(p_DHB_var{1} == constraints.start_pose(1:3));
    opti.subject_to(rvec_DHB_var{1} == constraints.start_pose(4:6));

    % Apply dynamic constraints using the single step method
    for k = 1:N-1
        [new_linear_frame, new_angular_frame, new_position, new_rotation] = ...
            reconstructTrajectorySingleStep(U_var(k,:), linear_frame, angular_frame, 'pos');

        % Update frames for next iteration
        linear_frame = new_linear_frame;
        angular_frame = new_angular_frame;

        % Set the dynamic constraints
        opti.subject_to(p_DHB_var{k+1} == new_position);
        opti.subject_to(rvec_DHB_var{k+1} == new_rotation);
    end

    % Constraints on the end pose
    opti.subject_to(p_DHB_var{end} == constraints.target_pose(1:3));
    opti.subject_to(rvec_DHB_var{end} == constraints.target_pose(4:6));

    % Construct objective
    objective = 0;
    for k=1:N
        e = U_var(k,:) - invariants_demo(k,:); % invariants error
        e_weighted = sqrt(params.weights).*e';
        objective = objective + e_weighted'*e_weighted;
    end

    % Initialize states and controls
    for k=1:N
        opti.set_initial(p_DHB_var{k}, init_traj(k,1:3));
        opti.set_initial(rvec_DHB_var{k}, init_traj(k,4:6));
        opti.set_initial(U_var(k,:), invariants_demo(k,:));
    end

    opti.minimize(objective);
    opti.solver('ipopt', struct(), struct('tol', 1e-5));

    % Solve the NLP
    sol = opti.solve();

    % get the updated invariants
    result.invariants_demo = invariants_demo;
    result.invariants = zeros(size(invariants_demo));
    result.invariants = sol.value(U_var);

    % compare the DHB invariants before and after optimization

    if (params.plot_comparison_invariants)
        % Plot the invariants
        figure('Name', 'Comparison on the invariants');

        dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2', ...
            'm_\omega' '\theta_\omega^1' '\theta_\omega^2'};
        for i=1:6
            subplot(6,1,i)
            plot(invariants_demo(:,i))
            hold on;
            plot(result.invariants(:,i))
            ylabel(dhbInvNames{i});
            legend('original', 'optimized(casadi)');
            grid on
        end
    end

    %-- get the state data
    pos_r_opt_data = zeros(N,3);
    rotm_r_opt_data = zeros(3,3,N);
    rvec_r_opt_data = zeros(N,3);
    for i=1:N
        pos_r_opt_data(i,:) = sol.value(p_DHB_var{i});
        rvec = sol.value(rvec_DHB_var{i});
        rotm = rotationVectorToMatrix(rvec);
        rotm_r_opt_data(:,:,i) = rotm;
        rvec_r_opt_data(i,:) = rvec;
    end

    % construct the result
    result.pos_data = pos_r_opt_data;
    result.rotm_data = rotm_r_opt_data;
    result.rvec_data = rvec_r_opt_data;
end
