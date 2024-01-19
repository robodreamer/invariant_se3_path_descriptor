%%% Testing the concept of DHB invariants
% Andy Park, Nov 2023

select_program = 7;

%%% 11-10-23
% 1) show the invariance over affine transforms
% 2) visualize the reconstructed path with different initial poses
% 3) fit a spline
if (select_program == 1)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw
    % note: you can also use marker data together with the function markers2pose.m

    load data/vive_data.mat
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    T_data = zeros(4,4,N);
    rvec_data = zeros(N, 3);
    qt_data = zeros(N, 4);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        angvec = rotm2axang(rotm); % angle and axis
        rvec = angvec(4) * angvec(1:3); % rotation vector
        T = [rotm, pos; 0 0 0 1];

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
        T_data(:,:,i) = T;
        qt_data(i,:) = qt;
    end

    % Get the initial transform
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    %% Compute DHB invariants

    % Compute position based DHB invariants
    [linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = ...
        computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    invariants_pos = [linear_motion_invariant, angular_motion_invariant];

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};
    for i=1:6
        subplot(2,3,i)
        plot(invariants_pos(:,i),'LineWidth',2)
        ylabel(dhbInvNames{i});
        grid on
    end

    %% Reconstruct the original trajectory

    % position
    [pos_r_data, rvec_r_data] = reconstructTrajectory(invariants_pos, linear_frame_initial, angular_frame_initial, 'pos');
    % [pos_r_data, rvec_r_data] = reconstructTrajectoryCasadi(invariants_pos, linear_frame_initial, angular_frame_initial, 'pos');

    % compute reconstruction errors
    errSP_pos = [(pos_r_data - pos_data(1:end-3,:)).^2 (rvec_r_data - rvec_data(1:end-3,:)).^2];

    % Compute rmse error
    err_pos = zeros(1,6);
    for i=1:6
        err_pos(i) = sum(errSP_pos(:,i));
    end

    RMSE_pos = sqrt([sum(err_pos(1:3)) sum(err_pos(4:6))]./(N-2));
    disp(['Reconstruction errors in pose (RMSE): ' num2str(RMSE_pos)])

    %% Plot the original and reconstructed paths

    % Plot original and reconstructed paths -- position
    pose_data = [rvec_data, pos_data]; % original
    pose_r_data = [rvec_r_data, pos_r_data]; % reconstructed

    figure('NumberTitle', 'off', 'Name', 'DHB to Cartesian pose');
    for i=1:6
        subplot(2,3,i)
        plot(pose_data(3:end,i),'g','LineWidth',4)
        hold on;
        plot(pose_r_data(:,i),'b','LineWidth',2)
        if(i<4)
            ylabel(['r_' num2str(i)])
        else
            ylabel(['p_' num2str(i-3)])
        end
        grid on
    end

    %% Plot the original and reconstructed paths in 3D

    % Convert the rotation vector data into rotation matrix and
    % reconstruct the transform data.
    rotm_r_data = zeros(3, 3, N-3);
    T_r_data = zeros(4, 4, N-3);
    qt_r_data = zeros(N-3, 4);
    for i=1:size(rvec_r_data,1)
        p_r = pos_r_data(i,:);
        rvec_r = rvec_r_data(i,:);
        rotm_r = rotationVectorToMatrix(rvec_r);
        qt_r = rotm2quat(rotm_r);
        T_r = [rotm_r, p_r'; 0 0 0 1];

        rotm_r_data(:,:,i) = rotm_r;
        T_r_data(:,:,i) = T_r;
        qt_r_data(i, :) = qt_r;
    end

    % plot 3D trajectory
    figure('NumberTitle', 'off', 'Name', 'Reconstructed Path with Affine Transform');
    plot3(pos_data(1:N-3,1), pos_data(1:N-3,2), pos_data(1:N-3,3), 'b');
    hold on;
    plot3(pos_r_data(:,1), pos_r_data(:,2), pos_r_data(:,3), 'r');
    plotTransforms(pos_r_data(1,:), qt_r_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_r_data(end,:), qt_r_data(end,:), 'FrameSize', 0.05);
    grid on;

    %% Apply affine transform to the pos data and show invariance

    % define affine transforms
    % pos_offset = 0.0*[1 2 3];
    pos_offset = -pos_data(1,:);
    eul = [pi/3 pi/5 pi/4];
    scale = 1;

    % apply the transform

    rotmZYX_offset = eul2rotm(eul);
    pos_tr_data = pos_data' + pos_offset';
    pos_tr_data = rotmZYX_offset * pos_tr_data;
    pos_tr_data = scale * pos_tr_data';

    % rotational transform on the rotational vector
    rvec_tr_data = (rotmZYX_offset * rvec_data')';

    % get position and rotation diff data
    pos_tr_diff_data = diff(pos_tr_data);

    % get the initial pose with transform
    T0_tr = T0;
    T0_tr(1:3,4) = pos_tr_data(1,:)';
    T0_tr(1:3,1:3) = rotmZYX_offset * rotm_data(1:3,1:3,1);

    % Compute position based DHB invariants
    [linear_motion_invariant_tr, angular_motion_invariant_tr, linear_frame_initial_tr, angular_frame_initial_tr] = computeDHB(pos_tr_diff_data, rvec_tr_data(1:end-1,:), 'pos', T0_tr);
    invariants_pos_tr = [linear_motion_invariant_tr, angular_motion_invariant_tr];

    % to check the scale invariance, normalize m in translation and rotation
    lin_inv_normalized = linear_motion_invariant;
    ang_inv_normalized = angular_motion_invariant;
    lin_inv_normalized(:, 1) = lin_inv_normalized(:, 1) / sum(lin_inv_normalized(:, 1));
    ang_inv_normalized(:, 1) = ang_inv_normalized(:, 1) / sum(ang_inv_normalized(:, 1));
    invariants_pos_normalized = [lin_inv_normalized, ang_inv_normalized];

    lin_inv_tr_normalized = linear_motion_invariant_tr;
    ang_inv_tr_normalized = angular_motion_invariant_tr;
    lin_inv_tr_normalized(:, 1) = lin_inv_tr_normalized(:, 1) / sum(lin_inv_tr_normalized(:, 1));
    ang_inv_tr_normalized(:, 1) = ang_inv_tr_normalized(:, 1) / sum(ang_inv_tr_normalized(:, 1));
    invariants_pos_tr_normalized = [lin_inv_tr_normalized, ang_inv_tr_normalized];

    % compute reconstruction errors
    err_invariant = [(lin_inv_tr_normalized - lin_inv_normalized).^2 (ang_inv_tr_normalized - ang_inv_normalized).^2];

    % Compute rmse error
    err_invariant_sum = zeros(1,6);
    for i=1:6
        err_invariant_sum(i) = sum(err_invariant(:,i));
    end

    RMSE_invariant = sqrt([sum(err_invariant_sum(1:3)) sum(err_invariant_sum(4:6))]./(N-3));
    disp(['Errors (RMSE) between invariants with affine transforms: ' num2str(RMSE_invariant)])

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB (w/ transform)');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};
    for i=1:6
        subplot(2,3,i)
        plot(invariants_pos_normalized(:,i),'r','LineWidth',2)
        hold on;
        plot(invariants_pos_tr_normalized(:,i),'b','LineWidth',2)
        ylabel(dhbInvNames{i});
        grid on
    end


    %% Reconstruction after affine transforms

    % position
    [pos_r_tr_data, rvec_r_tr_data] = reconstructTrajectory(invariants_pos_tr, linear_frame_initial_tr, angular_frame_initial_tr, 'pos');

    % Convert the rotation vector data into rotation matrix and
    % reconstruct the transform data.
    rotm_r_data = zeros(3, 3, N-3);
    T_r_data = zeros(4, 4, N-3);
    qt_r_data = zeros(N-3, 4);
    for i=1:size(rvec_r_data,1)
        p_r = pos_r_tr_data(i,:);
        rvec_r = rvec_r_tr_data(i,:);
        rotm_r = rotationVectorToMatrix(rvec_r)';
        qt_r = rotm2quat(rotm_r);
        T_r = [rotm_r, p_r'; 0 0 0 1];

        rotm_r_data(:,:,i) = rotm_r;
        T_r_data(:,:,i) = T_r;
        qt_r_data(i, :) = qt_r;
    end

    % plot 3D trajectory

    figure('NumberTitle', 'off', 'Name', 'Reconstructed Path with Affine Transform');

    % original path with starting at the origin
    pos_data_temp = pos_data - pos_data(1,:);
    plot3(pos_data_temp(1,1), pos_data_temp(1,2), pos_data_temp(1,3), 'o', 'MarkerSize', 10);
    hold on;
    plot3(pos_data_temp(1:N-3,1), pos_data_temp(1:N-3,2), pos_data_temp(1:N-3,3), 'b');
    plotTransforms(pos_data_temp(1,:), qt_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_data_temp(end,:), qt_data(end,:), 'FrameSize', 0.05);

    % reconstructed path
    plot3(pos_r_tr_data(:,1), pos_r_tr_data(:,2), pos_r_tr_data(:,3), 'r');
    plot3(pos_r_tr_data(1,1), pos_r_tr_data(1,2), pos_r_tr_data(1,3), 'o', 'MarkerSize', 10);
    plotTransforms(pos_r_tr_data(1,:), qt_r_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_r_tr_data(end,:), qt_r_data(end,:), 'FrameSize', 0.05);
    grid on;

    %% Fit a Spline to generate time-invariant data

    numData = size(invariants_pos, 1);
    tInit = 0;
    tFinal = 1;
    tDataOrig = linspace(tInit, tFinal, numData);

    %= fit a spline
    refTraj = spline(tDataOrig, invariants_pos');

    % get the velocity and acceleration splines
    pp = struct();
    pp.q = refTraj;
    pp.qd = ppDer(pp.q);
    pp.qdd = ppDer(pp.qd);

end

%%% 11-14-23
% 1) simple example of trajectory adaptation -- position only
% 2) normalize the cost function for each parameter
if (select_program == 2)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw
    % note: you can also use marker data together with the function markers2pose.m

    load data/vive_data.mat
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    T_data = zeros(4,4,N);
    rvec_data = zeros(N, 3);
    qt_data = zeros(N, 4);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        angvec = rotm2axang(rotm); % angle and axis
        rvec = angvec(4) * angvec(1:3); % rotation vector
        T = [rotm, pos; 0 0 0 1];

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
        T_data(:,:,i) = T;
        qt_data(i,:) = qt;
    end

    % Get the initial transform
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    %% Compute DHB invariants

    % Compute position based DHB invariants
    [pos_invariant_orig, rot_invariant_orig, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    invariants_pos_orig = [pos_invariant_orig, rot_invariant_orig];

    % Reconstruct the original trajectory and plot it in 3D

    % Modify the invariants
    modify_invariants = true;
    invariants_pos = invariants_pos_orig;
    temp_path = linspace(0.5, 1, size(invariants_pos, 1))';
    if (modify_invariants)
        % magnitute
        invariants_pos(:,1) = invariants_pos(:, 1) .* temp_path;
        % shape paramter 1
        % invariants_pos(:,2) = invariants_pos(:, 2) + temp_path;
        % shape parameter 2
        % invariants_pos(:,3) = invariants_pos(:, 3) + temp_path;
    end

    % position
    [pos_r_data, rvec_r_data] = reconstructTrajectory(invariants_pos, linear_frame_initial, angular_frame_initial, 'pos');

    % compute reconstruction errors
    errSP_pos = [(pos_r_data - pos_data(1:end-3,:)).^2 (rvec_r_data - rvec_data(1:end-3,:)).^2];

    % Compute rmse error
    err_pos = zeros(1,6);
    for i=1:6
        err_pos(i) = sum(errSP_pos(:,i));
    end

    RMSE_pos = sqrt([sum(err_pos(1:3)) sum(err_pos(4:6))]./(N-2));
    disp(['Reconstruction errors in pose (RMSE): ' num2str(RMSE_pos)])

    % Plot the original and reconstructed paths in 3D

    % Convert the rotation vector data into rotation matrix and
    % reconstruct the transform data.
    qt_r_data = rvecdata2qt(rvec_r_data);

    % plot 3D trajectory
    figure('NumberTitle', 'off', 'Name', 'Reconstructed Path with Modified Invariants');
    hold on;
    plotTransforms(pos_data(1,:), qt_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_data(N-3,:), qt_data(N-3,:), 'FrameSize', 0.05);
    plot3(pos_data(1:N-3,1), pos_data(1:N-3,2), pos_data(1:N-3,3), 'b');

    plot3(pos_r_data(:,1), pos_r_data(:,2), pos_r_data(:,3), 'r');
    plotTransforms(pos_r_data(1,:), qt_r_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_r_data(end,:), qt_r_data(end,:), 'FrameSize', 0.05);
    grid on;

    %% A simple 3D position target adaptation

    % set a desired target location
    pathsize = N-3;
    goal_orig = pos_data(pathsize, :);
    goal_offset = [-0.5, -1.5, 0.5];
    % goal_offset = zeros(1,3);
    goal_new = goal_orig + goal_offset;

    % create an init trajectory with a simple transformation
    pos_delta_data = SpatialRobotModel.jtraj(zeros(1,3), goal_offset, N);
    init_traj = pos_data + pos_delta_data;

    % plot 3D trajectory with original and new target positions
    plot_transformed_path = false;
    if (plot_transformed_path)
        figure('NumberTitle', 'off', 'Name', 'Reconstructed Path with Target Adaptation');
        hold on;
        plot3(pos_data(:,1), pos_data(:,2), pos_data(:,3), 'b');
        plot3(goal_orig(1), goal_orig(2), goal_orig(3), 'oy');
        plot3(goal_new(1), goal_new(2), goal_new(3), 'or');
        plot3(init_traj(:,1), init_traj(:,2), init_traj(:,3));
        view([-134.131 22.691]);

        grid on;
    end

    %% generate a GRP path that connects the initial path to the new target and optimize it

    Nframes = size(pos_data,1);

    % specify algorithms parameters
    algorithm_params = struct('init_sample_size',100,'sample_size',100,'max_iter',100,'h',10,'alpha',5e-1);

    % run trajectory optimization
    init_pose = pos_data(1, :)';
    final_pose = goal_new';

    tic;
    global quiet;
    quiet = 0;
    cost_fun = @(traj)cost_shape_descriptor(traj, rvec_data, T0, pos_invariant_orig);
    result = tromp_run(algorithm_params, cost_fun, init_pose, final_pose, Nframes, init_traj');
    elapsed_time = toc;  % Capture the elapsed time

    % Print out the elapsed time
    fprintf('Elapsed time is %.6f seconds.\n', elapsed_time);

    %% plot the results
    % plot the current trajectory
    plot_animation = 1;
    Niter = numel(result.learning_curve);
    if (plot_animation)
        figure('NumberTitle', 'off', 'Name', 'Progress trajectory in adaptation');
        for iter=1:Niter
            cla;
            curr_traj = result.traj_data{iter}';
            opt_cost = result.learning_curve(iter);
            plot3(pos_data(1:N-3,1), pos_data(1:N-3,2), pos_data(1:N-3,3),'g');
            hold on;
            plot3(init_traj(:,1), init_traj(:,2), init_traj(:,3),'m');
            plot3(curr_traj(:,1), curr_traj(:,2), curr_traj(:,3),'y');
            title(strcat('Optimized path ',sprintf(' [%d/%d], cost: (%.3f)',iter, Niter, opt_cost)));

            plot3(goal_orig(1), goal_orig(2), goal_orig(3), 'og');
            plot3(goal_new(1), goal_new(2), goal_new(3), 'om');
            legend('Original', 'Transformed', 'Optimized', 'Original Target', 'New Target');
            grid on;
            view([-134.131 22.691]);

            if(iter == 1)
                disp('press Enter to start');
                pause();
            else
                pause(1e-1);
                drawnow;
            end
        end
    end

    %% compare the DHB invariants before and after optimization
    % Compute DHB invariants for trajectory before optimization
    [pos_invariant_before, rot_invariant_before, linear_frame_initial, angular_frame_initial] = computeDHB(diff(init_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Compute DHB invariants for trajectory after optimization
    final_traj = result.traj_data{end}';
    [pos_invariant_after, rot_invariant_after, linear_frame_initial, angular_frame_initial] = computeDHB(diff(final_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2'};
    for i=1:3
        subplot(3,1,i)
        plot(invariants_pos_orig(:,i))
        hold on;
        plot(pos_invariant_before(:,i))
        plot(pos_invariant_after(:,i))
        ylabel(dhbInvNames{i});
        legend('original', 'transformed', 'optimized');
        grid on
    end


end

%%% 11-16-23
% 1) mexify the computeDHB method
% 2) compare the data with the one from eFSI representation
% 3) smooth the data a bit
if (select_program == 3)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw
    % note: you can also use marker data together with the function markers2pose.m

    load data/vive_data.mat
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    T_data = zeros(4,4,N);
    rvec_data = zeros(N, 3);
    qt_data = zeros(N, 4);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        angvec = rotm2axang(rotm); % angle and axis
        rvec = angvec(4) * angvec(1:3); % rotation vector
        T = [rotm, pos; 0 0 0 1];

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
        T_data(:,:,i) = T;
        qt_data(i,:) = qt;
    end

    % Get the initial transform
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    % Compute DHB invariants

    % Compute position based DHB invariants
    [pos_invariant_orig, rot_invariant_orig, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    invariants_pos_orig = [pos_invariant_orig, rot_invariant_orig];

    %% plot the original data and smoothed one

    % smooth the data a bit
    pos_data_orig = pos_data;
    pos_data_smooth = smoothdata(pos_data_orig, "sgolay", 5);

    % Compute position based DHB invariants
    [pos_invariant_smooth, rot_invariant_smooth, linear_frame_initial, angular_frame_initial] = computeDHB(diff(pos_data_smooth), rvec_data(1:end-1,:), 'pos', T0);
    invariants_pos_smooth = [pos_invariant_smooth, rot_invariant_smooth];

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2'};
    for i=1:3
        subplot(3,1,i)
        plot(pos_invariant_orig(:,i),'LineWidth',2)
        hold on;
        plot(pos_invariant_smooth(:,i),'LineWidth',2)
        ylabel(dhbInvNames{i});
        legend('Original', 'Smoothed')
        grid on
    end

    figure('NumberTitle', 'off', 'Name', 'Original Path with Smoothed Data');
    plot3(pos_data_orig(1:N-3,1), pos_data_orig(1:N-3,2), pos_data_orig(1:N-3,3));
    hold on;
    plot3(pos_data_smooth(:,1), pos_data_smooth(:,2), pos_data_smooth(:,3));
    legend('Original', 'Smoothed');
    grid on;
    view([-134.131 22.691]);

    %% A simple 3D position target adaptation

    use_smooth_data = true;
    if (use_smooth_data)
        pos_data = pos_data_smooth;
    else
        pos_data = pos_data_orig;
    end


    % set a desired target location
    pathsize = N-3;
    goal_orig = pos_data(pathsize, :);
    goal_offset = [-0.5, -1.5, 0.5];
    % goal_offset = zeros(1,3);
    goal_new = goal_orig + goal_offset;

    % create an init trajectory with a simple transformation
    pos_delta_data = SpatialRobotModel.jtraj(zeros(1,3), goal_offset, N);
    init_traj = pos_data + pos_delta_data;

    %% generate a GRP path that connects the initial path to the new target and optimize it

    Nframes = size(pos_data,1);

    % set the seed
    % rng(2423423);

    % specify algorithms parameters
    % h - sensitivity
    % alpha - stepsize
    algorithm_params = struct('sample_size',50,'max_iter',100,'h',1e1,'alpha',5e-1);
    weights = [2 1 1];
    weights = weights/norm(weights);

    % run trajectory optimization
    init_pose = pos_data(1, :)';
    final_pose = goal_new';

    % check the elapsed time

    tic;
    global quiet;
    quiet = 0;
    cost_fun = @(traj)cost_shape_descriptor_mex(traj, rvec_data, T0, pos_invariant_orig, weights);
    result = tromp_run(algorithm_params, cost_fun, init_pose, final_pose, Nframes, init_traj');
    elapsed_time = toc;  % Capture the elapsed time

    % Print out the elapsed time
    fprintf('Elapsed time is %.6f seconds.\n', elapsed_time);

    %% plot the data comparing to the result with eFSI representation

    % get the optimal trajectory
    grp_traj = result.opt_traj';

    % get the data from eFSI representation
    load eFSI_result_231117.mat
    eFSI_traj = optim_gen_result.Obj_location;

    figure('NumberTitle', 'off', 'Name', 'Generalized Trajectory (DHB vs eFSI)');
    plot3(pos_data(1:N-3,1), pos_data(1:N-3,2), pos_data(1:N-3,3));
    hold on;
    plot3(init_traj(:,1), init_traj(:,2), init_traj(:,3));
    plot3(grp_traj(:,1), grp_traj(:,2), grp_traj(:,3));
    plot3(eFSI_traj(:,1), eFSI_traj(:,2), eFSI_traj(:,3));

    plot3(goal_orig(1), goal_orig(2), goal_orig(3), 'og');
    plot3(goal_new(1), goal_new(2), goal_new(3), 'om');
    legend('Original', 'Transformed', 'Optimized (GRP)', 'Optimized (eFSI)', 'Original Target', 'New Target');
    grid on;
    view([-134.131 22.691]);

    %% compare the DHB invariants before and after optimization
    % Compute DHB invariants for trajectory before optimization
    [pos_invariant_before, rot_invariant_before, linear_frame_initial, angular_frame_initial] = computeDHB(diff(init_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Compute DHB invariants for trajectory after optimization
    final_traj = result.opt_traj';
    [pos_invariant_after, rot_invariant_after, linear_frame_initial, angular_frame_initial] = computeDHB(diff(final_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Compute DHB invariants for trajectory after optimization with eFSI
    [pos_invariant_after2, rot_invariant_after2, linear_frame_initial, angular_frame_initial] = computeDHB(diff(eFSI_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2'};
    for i=1:3
        subplot(3,1,i)
        plot(invariants_pos_orig(:,i))
        hold on;
        plot(pos_invariant_before(:,i))
        plot(pos_invariant_after(:,i))
        plot(pos_invariant_after2(:,i))
        ylabel(dhbInvNames{i});
        legend('original', 'transformed', 'optimized(GRP)', 'optimized(eFSI)');
        grid on
    end

    %% reconstruct the data from eFSI
    % position
    [pos_r_tr_data, rvec_r_tr_data] = ...
        reconstructTrajectory([pos_invariant_after2, rot_invariant_after2], ...
        linear_frame_initial, angular_frame_initial, 'pos');

    % Convert the rotation vector data into rotation matrix and
    % reconstruct the transform data.
    rotm_r_data = zeros(3, 3, N-3);
    T_r_data = zeros(4, 4, N-3);
    qt_r_data = zeros(N-3, 4);
    for i=1:size(pos_invariant_after2,1)
        p_r = pos_r_tr_data(i,:);
        rvec_r = rvec_r_tr_data(i,:);
        rotm_r = rotationVectorToMatrix(rvec_r)';
        qt_r = rotm2quat(rotm_r);
        T_r = [rotm_r, p_r'; 0 0 0 1];

        rotm_r_data(:,:,i) = rotm_r;
        T_r_data(:,:,i) = T_r;
        qt_r_data(i, :) = qt_r;
    end

    % plot 3D trajectory

    figure('NumberTitle', 'off', 'Name', 'Reconstructed Path with eFSI data');

    % original path with starting at the origin
    pos_data_temp = pos_data - pos_data(1,:);
    plot3(pos_data_temp(1,1), pos_data_temp(1,2), pos_data_temp(1,3), 'o', 'MarkerSize', 10);
    hold on;
    plot3(pos_data_temp(1:N-3,1), pos_data_temp(1:N-3,2), pos_data_temp(1:N-3,3));
    plotTransforms(pos_data_temp(1,:), qt_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_data_temp(end,:), qt_data(end,:), 'FrameSize', 0.05);

    % reconstructed path starting at the origin
    pos_r_tr_data = pos_r_tr_data - pos_r_tr_data(1,:);
    plot3(pos_r_tr_data(:,1), pos_r_tr_data(:,2), pos_r_tr_data(:,3));
    plot3(pos_r_tr_data(1,1), pos_r_tr_data(1,2), pos_r_tr_data(1,3), 'o', 'MarkerSize', 10);
    plotTransforms(pos_r_tr_data(1,:), qt_r_data(1,:), 'FrameSize', 0.05);
    plotTransforms(pos_r_tr_data(end,:), qt_r_data(end,:), 'FrameSize', 0.05);
    grid on;


end

%%% Tue Nov 21 05:37:56 AM EST 2023
% 1) update the cost function
% 2) use a resampled path
% 3) use another plotting function for trajectory
if (select_program == 4)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw
    % note: you can also use marker data together with the function markers2pose.m

    load data/vive_data.mat
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % get the data from eFSI representation
    load eFSI_result_231117.mat
    eFSI_pos_traj_orig = optim_gen_result.Obj_location;
    eFSI_rot_traj_orig = optim_gen_result.Obj_frames;

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    rvec_data = zeros(N, 3);
    rvec_eFSI_data = zeros(N,3);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        rvec = rotationMatrixToVector(rotm);

        rotm_eFSI = eFSI_rot_traj_orig(:,:,i);
        rvec_eFSI_data(i,:) = rotationMatrixToVector(rotm_eFSI);

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
    end
    rotm_data_orig = rotm_data;
    rvec_eFSI_data_orig = rvec_eFSI_data;

    % Get the initial transform
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];

    % resample the path -- 100 points
    Nframes = 100;
    timeNew = linspace(1, size(pos_data,1), Nframes);
    pos_data = SpatialRobotModel.cubicInterp(pos_data, 1:N, timeNew);
    rvec_data = SpatialRobotModel.cubicInterp(rvec_data, 1:N, timeNew);
    eFSI_pos_traj = SpatialRobotModel.cubicInterp(eFSI_pos_traj_orig, 1:N, timeNew);
    rvec_eFSI_data = SpatialRobotModel.cubicInterp(rvec_eFSI_data, 1:N, timeNew);
    N = Nframes;

    % resize the data
    rotm_data = zeros(3,3,N);
    rotm_eFSI_data = zeros(3,3,N);
    for i=1:N
        rvec = rvec_data(i,:);
        rotm = rotationVectorToMatrix(rvec)';
        rotm_data(:,:,i) = rotm;

        rvec_eFSI = rvec_eFSI_data(i,:);
        rotm_eFSI = rotationVectorToMatrix(rvec_eFSI)';
        rotm_eFSI_data(:,:,i) = rotm_eFSI;
    end
    eFSI_rot_traj = rotm_eFSI_data;

    % smooth the data a bit
    pos_data_orig = pos_data;
    pos_data = smoothdata(pos_data, "sgolay", 5);

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    % Compute position based DHB invariants
    [pos_invariant, rot_invariant, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    invariants_pos = [pos_invariant, rot_invariant];

    %% A simple 3D position target adaptation

    % set a desired target location
    pathsize = N-3;
    goal_orig = pos_data(pathsize, :);
    goal_offset = [-0.5, -1.5, 0.5];
    % goal_offset = zeros(1,3);
    goal_new = goal_orig + goal_offset;

    % create an init trajectory with a simple transformation
    pos_delta_data = SpatialRobotModel.jtraj(zeros(1,3), goal_offset, N);
    init_traj = pos_data + pos_delta_data;
    % init_traj = eFSI_pos_traj;

    %%% generate a GRP path that connects the initial path to the new target and optimize it

    Nframes = size(pos_data,1);

    % set the seed
    % rng(2423423);

    % specify algorithms parameters
    % h - sensitivity
    % alpha - stepsize
    algorithm_params = struct('sample_size',50,'max_iter',100,'h',1e1,'alpha',5e-1);
    weights = [2 1 1 2];  % the last weight for the path length
    weights = weights/norm(weights);

    % run trajectory optimization
    init_pose = pos_data(1, :)';
    final_pose = goal_new';

    % check the elapsed time

    tic;
    global quiet;
    quiet = 1;
    cost_fun = @(traj)cost_shape_descriptor_mex2(traj, rvec_data, T0, pos_invariant, weights);
    result = tromp_run(algorithm_params, cost_fun, init_pose, final_pose, Nframes, init_traj');
    elapsed_time = toc;  % Capture the elapsed time

    % Print out the elapsed time
    fprintf('Elapsed time is %.6f seconds.\n', elapsed_time);

    %%% plot the data comparing to the result with eFSI representation

    % get the optimal trajectory
    grp_traj = result.opt_traj';

    select_plot_method = 2;

    if (select_plot_method == 1)
        figure('NumberTitle', 'off', 'Name', 'Generalized Trajectory (DHB vs eFSI)');
        plot3(pos_data(1:N-3,1), pos_data(1:N-3,2), pos_data(1:N-3,3));
        hold on;
        plot3(init_traj(:,1), init_traj(:,2), init_traj(:,3));
        plot3(grp_traj(:,1), grp_traj(:,2), grp_traj(:,3));
        plot3(eFSI_pos_traj(:,1), eFSI_pos_traj(:,2), eFSI_pos_traj(:,3));

        plot3(goal_orig(1), goal_orig(2), goal_orig(3), 'og');
        plot3(goal_new(1), goal_new(2), goal_new(3), 'om');
        legend('Original', 'Transformed', 'Optimized (GRP)', 'Optimized (eFSI)', 'Original Target', 'New Target');
        grid on;
        view([-134.131 22.691]);

    elseif (select_plot_method == 2)
        % trajectory 1
        grp_path = struct();
        grp_path.pos_data = grp_traj;
        grp_path.rot_data = rotm_data;

        % trajectory 2
        eFSI_path = struct();
        eFSI_path.pos_data = eFSI_pos_traj;
        eFSI_path.rot_data = eFSI_rot_traj;

        plot_trajectory(grp_path, eFSI_path, 'Result with DHB-GRP (blue), Result with eFSI(red)', true);

    end

    %% compare the DHB invariants before and after optimization
    % Compute DHB invariants for trajectory before optimization
    [pos_invariant_before, rot_invariant_before, linear_frame_initial, angular_frame_initial] = computeDHB(diff(init_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Compute DHB invariants for trajectory after optimization
    final_traj = result.opt_traj';
    [pos_invariant_after, rot_invariant_after, linear_frame_initial, angular_frame_initial] = computeDHB(diff(final_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Compute DHB invariants for trajectory after optimization with eFSI
    [pos_invariant_after2, rot_invariant_after2, linear_frame_initial, angular_frame_initial] = computeDHB(diff(eFSI_pos_traj), rvec_data(1:end-1,:), 'pos', T0);

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2'};
    for i=1:3
        subplot(3,1,i)
        plot(invariants_pos(:,i))
        hold on;
        plot(pos_invariant_before(:,i))
        plot(pos_invariant_after(:,i))
        plot(pos_invariant_after2(:,i))
        ylabel(dhbInvNames{i});
        legend('original', 'transformed', 'optimized(GRP)', 'optimized(eFSI)');
        grid on
    end

end

%%% Mon Nov 27 05:18:09 AM EST 2023
% 1) try using an OCP approach to solve trajectory adaptation problem
% 2) compare the results with GRP, eFSI
if (select_program == 5)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw
    % note: you can also use marker data together with the function markers2pose.m

    load data/vive_data.mat
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % get the data from eFSI representation
    load eFSI_result_231117.mat
    eFSI_pos_traj_orig = optim_gen_result.Obj_location;
    eFSI_rot_traj_orig = optim_gen_result.Obj_frames;

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    rvec_data = zeros(N, 3);
    rvec_eFSI_data = zeros(N,3);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        rvec = rotationMatrixToVector(rotm);

        rotm_eFSI = eFSI_rot_traj_orig(:,:,i);
        rvec_eFSI_data(i,:) = rotationMatrixToVector(rotm_eFSI);

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
    end
    rotm_data_orig = rotm_data;
    rvec_eFSI_data_orig = rvec_eFSI_data;

    % Get the initial transform
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];

    % resample the path -- 50 points
    Nframes = 50;
    timeNew = linspace(1, size(pos_data,1), Nframes);
    pos_data = SpatialRobotModel.cubicInterp(pos_data, 1:N, timeNew);
    rvec_data = SpatialRobotModel.cubicInterp(rvec_data, 1:N, timeNew);
    eFSI_pos_traj = SpatialRobotModel.cubicInterp(eFSI_pos_traj_orig, 1:N, timeNew);
    rvec_eFSI_data = SpatialRobotModel.cubicInterp(rvec_eFSI_data, 1:N, timeNew);
    N_orig = Nframes;
    N = N_orig;

    % resize the data
    rotm_data = zeros(3,3,N);
    rotm_eFSI_data = zeros(3,3,N);
    for i=1:N
        rvec = rvec_data(i,:);
        rotm = rotationVectorToMatrix(rvec)';
        rotm_data(:,:,i) = rotm;

        rvec_eFSI = rvec_eFSI_data(i,:);
        rotm_eFSI = rotationVectorToMatrix(rvec_eFSI)';
        rotm_eFSI_data(:,:,i) = rotm_eFSI;
    end
    eFSI_rot_traj = rotm_eFSI_data;

    % smooth the data a bit
    pos_data_orig = pos_data;
    pos_data = smoothdata(pos_data, "sgolay", 5);

    % get the pos diff data
    pos_diff_data = diff(pos_data);

    % Compute position based DHB invariants
    [pos_invariant, rot_invariant, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    path_invariants = [pos_invariant, rot_invariant];

    % transform the goal position
    N = N_orig - 3;

    % set a desired target location
    goal_orig = pos_data(N, :);
    goal_offset = [-0.5, -1.5, 0.5];
    goal_new = goal_orig + goal_offset;

    %% Trajectory adaptation with GRP

    % create an init trajectory with a simple transformation
    pos_delta_data = SpatialRobotModel.jtraj(zeros(1,3), goal_offset, N_orig);
    init_traj = pos_data + pos_delta_data;
    % init_traj = eFSI_pos_traj;

    %%% generate a GRP path that connects the initial path to the new target and optimize it

    Nframes = size(pos_data,1);

    % set the seed
    % rng(2423423);

    % specify algorithms parameters
    % h - sensitivity
    % alpha - stepsize
    algorithm_params = struct('sample_size',50,'max_iter',100,'h',1e1,'alpha',5e-1);
    weights = [2 1 1 2];  % the last weight for the path length
    weights = weights/norm(weights);

    % run trajectory optimization
    init_pose = pos_data(1, :)';
    final_pose = goal_new';

    % check the elapsed time

    tic;
    global quiet;
    quiet = 1;
    cost_fun = @(traj)cost_shape_descriptor_mex2(traj, rvec_data, T0, pos_invariant, weights);
    result = tromp_run(algorithm_params, cost_fun, init_pose, final_pose, Nframes, init_traj');
    elapsed_time = toc;  % Capture the elapsed time

    % Print out the elapsed time
    fprintf('Elapsed time is %.6f seconds.\n', elapsed_time);

    %%% plot the data comparing to the result with eFSI representation

    % get the optimal trajectory
    grp_traj = result.opt_traj';

    %% Try trajectory adaptation with OCP

    %%% construct an ocp with casadi

    % Import statement that loads the Casadi module with its functions
    import casadi.*
    clc;

    opti = casadi.Opti();

    % Initial values
    init_traj = [pos_data(1:N, :), rvec_data(1:N, :)];
    invariants_demo = path_invariants;
    p_frame_init = linear_frame_initial;
    R_frame_init = angular_frame_initial;

    % constraints
    constraints = struct();
    constraints.start_pose = init_traj(1,:);
    constraints.target_pose = [goal_new, rvec_data(N,:)];

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

    % Constraints on the start and end pose
    opti.subject_to(p_DHB_var{1} == constraints.start_pose(1:3));
    opti.subject_to(rvec_DHB_var{1} == constraints.start_pose(4:6));
    opti.subject_to(p_DHB_var{end} == constraints.target_pose(1:3));
    opti.subject_to(rvec_DHB_var{end} == constraints.target_pose(4:6));

    % Construct objective
    weights = [1, 1, 1, 1, 1, 1]';
    % Construct objective
    objective = 0;
    for k=1:N
        e = U_var(k,:) - invariants_demo(k,:); % invariants error
        e_weighted = sqrt(weights).*e';
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
    optim_result = struct();
    optim_result.invariants = zeros(size(invariants_demo));
    optim_result.invariants = sol.value(U_var);

    %% compare the DHB invariants before and after optimization

    % Plot the invariants
    figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');

    dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2'};
    for i=1:3
        subplot(3,1,i)
        plot(invariants_demo(:,i))
        hold on;
        plot(optim_result.invariants(:,i))
        ylabel(dhbInvNames{i});
        legend('original', 'optimized(casadi)');
        grid on
    end

    %% Plot the trajectory after reconstruction

    % reconstruct the data
    [pos_r_opt_data, rvec_r_opt_data] = reconstructTrajectory(optim_result.invariants, linear_frame_initial, angular_frame_initial, 'pos');

    % get the rotation matrices
    rotm_r_opt_data = zeros(3,3,N);
    for i=1:N
        rvec = rvec_r_opt_data(i,:);
        rotm = rotationVectorToMatrix(rvec)';
        rotm_r_opt_data(:,:,i) = rotm;
    end

    % trajectory 1
    path_first = struct();
    path_first.pos_data = pos_data;
    path_first.rot_data = rotm_data;

    % trajectory 2
    path_second = struct();
    path_second.pos_data = grp_traj;
    path_second.rot_data = rotm_data;

    % trajectory 3
    path_third = struct();
    path_third.pos_data = pos_r_opt_data;
    path_third.rot_data = rotm_r_opt_data;

    % trajectory 4
    path_fourth = struct();
    path_fourth.pos_data = eFSI_pos_traj;
    path_fourth.rot_data = eFSI_rot_traj;

    path_data = {path_first, path_second, path_third, path_fourth};
    color_data = random_color(size(path_data,2),'jet',1212);
    legend_texts = {'The initial Traj', 'Result (GRP)', 'Result (DHB-NLP)', 'Result (eFSI-NLP)'};
    params = struct('auto_calculate_scale', true, 'show_rotation', false);

    plot_se3_trajectories(path_data, color_data, 'Comparison on the results with the demo path', legend_texts, params);
end

%%% Thu Nov 30 05:18:09 AM EST 2023
% 1) refactor the trajectory adaptation with DHB-NLP
% 1) a simple unit test for the result
if (select_program == 6)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    use_randomly_generated_data = false;
    if (use_randomly_generated_data)
        N = 100; % Number of poses
        numWaypoints = 5; % Number of waypoints
        positionRange = [-1, 1]; % Range for x, y, z coordinates
        time = linspace(0, N, N);

        % get the generated se3 path
        measured_pose_coordinates = [time', generateRandomComplexSE3Path(N, numWaypoints, positionRange)];
    else
        % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw

        load data/vive_data.mat
        % load data/pouring_motion_data.mat
        % load data/scooping_motion_data.mat
        % load data/curved_motion_data.mat
    end

    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    rvec_data = zeros(N, 3);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        rvec = rotationMatrixToVector(rotm);

        % store the data
        rotm_data(:,:,i) = rotm';
        rvec_data(i, :) = rvec;
    end
    rotm_data_orig = rotm_data;

    % Get the initial transform
    T_init_orig = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];
    T_final_orig = [rotm_data(:,:,end-3) pos_data(end-3,:)'; 0 0 0 1];

    %- call the refactored method for trajectory adaptation

    close all;

    % Define input parameters
    inputPoseData = struct('pos_data', pos_data, 'rotm_data', rotm_data);
    params = struct('Nframes', 50, 'smoothing', false, ...
        'plot_comparison_invariants', true, 'weights', ones(6,1));

    % initialize unit test
    testCase = matlab.unittest.TestCase.forInteractiveUse;
    tol = 1e-2*ones(4,4); % this can be reduced futher later

    %% each unit test
    disp('Test 1 - zero offsets')

    % construct initial and final pose
    pose_offset_init = [0, 0, 0, 0, 0, 0];  % [m] and [deg]
    pose_offset_final = [0, 0, 0, 0, 0, 0]; % [m] and [deg]
    T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
    T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

    % Generate adapted trajectory
    result = generate_trajectory(inputPoseData, params, T_init, T_final);

    % Check the result
    T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
    T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

    verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
    verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

    %% each unit test
    disp('Test 2 - only positional offsets')

    % construct initial and final pose
    pose_offset_init = [0.3, -0.2, 0.1, 0, 0, 0];  % [m] and [deg]
    pose_offset_final = [-0.2, 0.3, 0.3, 0, 0, 0]; % [m] and [deg]
    T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
    T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

    % Generate adapted trajectory
    result = generate_trajectory(inputPoseData, params, T_init, T_final);

    % Check the result
    T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
    T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

    verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
    verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

    %% each unit test
    disp('Test 3 -  offsets in initial pose')

    % construct initial and final pose
    pose_offset_init = [0.3, -0.2, 0.1, 30, 20, 40];  % [m] and [deg]
    pose_offset_final = [0, 0, 0, 0, 0, 0]; % [m] and [deg]
    T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
    T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

    % Generate adapted trajectory
    result = generate_trajectory(inputPoseData, params, T_init, T_final);

    % Check the result
    T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
    T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

    verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
    verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

    %% each unit test
    disp('Test 4 -  offsets in final pose')

    % construct initial and final pose
    pose_offset_init = [0, 0, 0, 0, 0, 0];  % [m] and [deg]
    pose_offset_final = [0.3, -0.2, 0.3, 45, 60, 30]; % [m] and [deg]
    T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
    T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

    % Generate adapted trajectory
    result = generate_trajectory(inputPoseData, params, T_init, T_final);

    % Check the result
    T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
    T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

    verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
    verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

    %% each unit test
    disp('Test 5 -  offsets in both poses')

    % construct initial and final pose
    pose_offset_init = [0.3, -0.1, 0.2, 30, -20, 45];  % [m] and [deg]
    pose_offset_final = [0.2, 0.3, 0.5, -45, 30, 10]; % [m] and [deg]
    T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
    T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

    % Generate adapted trajectory
    result = generate_trajectory(inputPoseData, params, T_init, T_final);

    % Check the result
    T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
    T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

    verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
    verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

    %% Plot the result trajectory after adaptation

    close all;

    % trajectory 1
    path_first = struct();
    path_first.pos_data = result.pos_data_demo(1:end-3,:);
    path_first.rot_data = result.rotm_data_demo(:,:,1:end-3);

    % trajectory 2
    path_second = struct();
    path_second.pos_data = result.pos_data;
    path_second.rot_data = result.rotm_data;

    path_data = {path_first, path_second};
    color_data = random_color(size(path_data,2),'jet',1232);
    legend_texts = {'The initial Traj', 'The Adapted Traj (DHB-NLP)'};
    params_plot = struct('auto_calculate_scale', true, 'scale', 2, ...
        'show_rotation', true, 'show_coordinates', true);
    handle_paths = plot_se3_trajectories(path_data, color_data, ...
        'Comparison on the results with the demo path', params_plot);

    % overlay T_final
    SpatialRobotModel.drawFrame(T_final, 0.2, 3);

    legend(handle_paths, legend_texts, 'Location', 'best');
end

%%% Sat Dec 2
% 1) try a variety of changes and data
if (select_program == 7)

    %% Data Preparation

    addpath(genpath(pwd));
    ccc;

    use_randomly_generated_data = true;
    if (use_randomly_generated_data)
        N = 100; % Number of time steps
        numWaypoints = 5; % Number of waypoints
        positionRange = [-1, 1]; % Range for x, y, z coordinates
        time = linspace(0, N, N);

        % get the generated se3 path
        measured_pose_coordinates = [time', generateRandomComplexSE3Path(N, numWaypoints, positionRange)];
    else
        % loading pose data, format: rows=samples, columns=x|y|z|qx|qy|qz|qw

        load data/vive_data.mat
    end
    N = length(measured_pose_coordinates);
    dt = 1/60; % timestep

    % Convert quaternion to rotation matrix
    pos_data = measured_pose_coordinates(:,2:4); % position
    rotm_data = zeros(3,3,N);
    rvec_data = zeros(N, 3);
    for i=1:N
        % [qw, qx, qy, qz]
        pos = pos_data(i,:)';
        qt = [measured_pose_coordinates(i,8), measured_pose_coordinates(i,5:7)];
        rotm = quat2rotm(qt); % rotation matrix
        rvec = rotationMatrixToVector(rotm);

        % store the data
        rotm_data(:,:,i) = rotm;
        rvec_data(i, :) = rvec;
    end
    rotm_data_orig = rotm_data;

    % Get the initial transform
    T_init_orig = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];
    T_final_orig = [rotm_data(:,:,end) pos_data(end,:)'; 0 0 0 1];

    %% run trajectory adaptation for multiple transforms

    close all;
    % Define input parameters
    inputPoseData = struct('pos_data', pos_data, 'rotm_data', rotm_data);
    params = struct('Nframes', 50, 'smoothing', true, ...
        'plot_comparison_invariants', false, 'weights', [1 1 1 1 1 1]');

    % loop for different transforms to the init and target poses
    numTests = 10;

    % trajectory original
    path_data = cell(1, numTests + 1);
    path_data{1} = struct();
    path_data{1}.pos_data = pos_data;
    path_data{1}.rot_data = rotm_data;
    legend_texts = {'The initial Traj'};

    % initialize unit test
    testCase = matlab.unittest.TestCase.forInteractiveUse;
    tol = 1e-2*ones(4,4); % this can be reduced futher later
    HALF_RANGE_POS = 0.3;
    HALF_RANGE_ROT = 45;

    % rng(1235);
    result_data = cell(numTests, 1);
    for k = 1:numTests

        fprintf('Generating %d/%d-th trajectory!\n', k, numTests);
        % Generate random translation and rotation values
        posOffset_init = SpatialRobotModel.randRange(3, -HALF_RANGE_POS, HALF_RANGE_POS); % Random translation in the range
        posOffset_final = SpatialRobotModel.randRange(3, -HALF_RANGE_POS, HALF_RANGE_POS); % Random translation in the range
        rotOffset_init = SpatialRobotModel.randRange(3, -HALF_RANGE_ROT, HALF_RANGE_ROT); % Random rotation in the range
        rotOffset_final = SpatialRobotModel.randRange(3, -HALF_RANGE_ROT, HALF_RANGE_ROT); % Random rotation in the range

        % Apply random translation and rotation to initial and final transforms
        pose_offset_init = [posOffset_init, rotOffset_init];
        pose_offset_final = [posOffset_final, rotOffset_final];

        T_init = constructPoseWithOffset(T_init_orig, pose_offset_init);
        T_final = constructPoseWithOffset(T_final_orig, pose_offset_final);

        % Generate adapted trajectory
        result = generate_trajectory(inputPoseData, params, T_init, T_final);

        % Check the result
        T_init_out = [result.rotm_data(:,:,1), result.pos_data(1,:)'; 0 0 0 1];
        T_final_out = [result.rotm_data(:,:,end), result.pos_data(end,:)'; 0 0 0 1];

        verifyEqual(testCase, T_init, T_init_out, 'AbsTol', tol);
        verifyEqual(testCase, T_final, T_final_out, 'AbsTol', tol);

        % store the result data
        path_data{k + 1} = struct();
        path_data{k + 1}.pos_data = result.pos_data;
        path_data{k + 1}.rot_data = result.rotm_data;
        legend_texts{k + 1} = sprintf('%d-th Adaptated Path', k);
        result_data{k} = result;
    end

    %-- Plot the trajectory after reconstruction
    color_data = random_color(size(path_data,2),'jet',1232);

    params = struct('auto_calculate_scale', false, 'scale', 2, ...
        'show_rotation', true, 'show_coordinates', false);
    handle_paths = plot_se3_trajectories(path_data, color_data, ...
        'Comparison on the results with the demo path', params);

    legend(handle_paths, legend_texts, 'Location', 'best');

    % -- Plot the invariants and se3 data
    plot_invariants = true;
    if (plot_invariants)
        % Plot the invariants
        figure('NumberTitle', 'off', 'Name', 'DHB Invariants');

        dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};
        for k=1:numTests
            invariants_pos = result_data{k}.invariants;
            for i=1:6
                subplot(2,3,i)
                plot(invariants_pos(:,i),'LineWidth',2)
                hold on;
                ylabel(dhbInvNames{i});
                grid on
            end
        end

        % plot SE3 coordinate paths
        figure('NumberTitle', 'off', 'Name', 'SE3 Paths');
        coordNames = {'x' 'y' 'z' 'rx' 'ry' 'rz'};
        for k=1:numTests
            for i=1:6
                subplot(2,3,i)
                if (i < 4)
                    plot(result_data{k}.pos_data(:,i),'LineWidth',2)
                else
                    plot(result_data{k}.rvec_data(:,i-3),'LineWidth',2)
                end
                hold on;
                ylabel(coordNames{i});
                grid on
            end
        end
    end
end

%%% Fri Jan 19 10:27:56 AM EST 2024
% 1) visualize the data from Franka Teleop
% 2) compare with the DHB invariants computed in Python-implementation
if (select_program == 8)

%% Data Preparation

addpath(genpath(pwd));
ccc;

% load data/data_pouring_011924.mat
load data/data_fruits_011924.mat

demo_list = fieldnames(trajectories);
numTests = numel(demo_list);

%-- store the path data
legend_texts = cell(numTests,1);
rvec_prev = [1 1 1];
for k = 1:numTests

    pos_data = trajectories.(demo_list{k}).positions;
    quat_data = trajectories.(demo_list{k}).quaternions;
    num_data = size(pos_data,1);

    % quaternion to rotm data
    rotm_data = zeros(3,3,num_data);
    rvec_data = zeros(num_data, 3);
    for idx=1:num_data
        % [qw, qx, qy, qz]
        qt = [quat_data(idx,4), quat_data(idx,1:3)];
        rotm = quat2rotm(qt); % rotation matrix
        rvec = rotationMatrixToVector(rotm);

        % make the signs consistent in rvec
        flipped = false;
        if idx == 1
           flipped = rvec(1) < 0;
        else
           flipped = dot(rvec, rvec_data(idx-1,:)) < 0;
        end

        if flipped
            rvec = -rvec;
        end

        % store the data
        rotm_data(:,:,idx) = rotm;
        rvec_data(idx, :) = rvec;
    end

    % store the result data
    path_data{k} = struct();
    path_data{k}.pos_data = pos_data;
    path_data{k}.rot_data = rotm_data;
    path_data{k}.rvec_data = rvec_data;
    legend_texts{k} = sprintf('%d-th SE3 Path', k);
end

%-- Plot the SE3 trajectories in 3D
plot_se3_trajectories = false;
if (plot_se3_trajectories)
    color_data = random_color(size(path_data,2),'jet',1232);

    params = struct('auto_calculate_scale', false, 'scale', 0.3, ...
        'show_rotation', true, 'show_coordinates', false);
    handle_paths = plot_se3_trajectories(path_data, color_data, ...
        'The SE3 paths from demonstrations', params);

    legend(handle_paths, legend_texts, 'Location', 'best');
end

%% Pre-processing the data
disp('Pre-processing the data...');

% make the path lengths the same
resample_paths = true;
num_samples = 250;
smoothing_enabled = true;
smoothing_factor = 10; % 10 by default

if resample_paths
    for k = 1:numTests
        pos_data_orig = path_data{k}.pos_data;
        rvec_data_orig = path_data{k}.rvec_data;
        rotm_data_orig = path_data{k}.rot_data;

        % resample the path
        num_data = size(pos_data_orig,1);
        timeNew = linspace(1, num_data, num_samples-3);
        timeNew = [timeNew timeNew(end) timeNew(end) timeNew(end)];
        pos_data = SpatialRobotModel.cubicInterp(pos_data_orig, 1:num_data, timeNew);
        rvec_data = SpatialRobotModel.cubicInterp(rvec_data_orig, 1:num_data, timeNew);

        % resize the data
        rotm_data = zeros(3,3,num_samples);
        for i=1:num_samples
            rvec = rvec_data(i,:);
            rotm = rotationVectorToMatrix(rvec);
            rotm_data(:,:,i) = rotm;
        end

        % smooth the data a bit
        if smoothing_enabled
            pos_data = smoothdata(pos_data, "sgolay", smoothing_factor);
            rvec_data = smoothdata(rvec_data, "sgolay", smoothing_factor);
        end

        % store the modified data
        path_data{k}.pos_data = pos_data;
        path_data{k}.rot_data = rotm_data;
        path_data{k}.rvec_data = rvec_data;
    end
end

% compute DHB invariants from MATLAB side
invariants_pos_matlab_data = cell(k,1);
for k=1:numTests
    pos_data = path_data{k}.pos_data;
    rotm_data = path_data{k}.rot_data;
    rvec_data = path_data{k}.rvec_data;

    pos_diff_data = diff(path_data{k}.pos_data);
    T0 = [rotm_data(:,:,1) pos_data(1,:)'; 0 0 0 1];
    [linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = ...
        computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
    invariants_data_matlab{k} = struct();
    invariants_data_matlab{k}.invariants_pos = [linear_motion_invariant, angular_motion_invariant];
    invariants_data_matlab{k}.linear_frame_initial = linear_frame_initial;
    invariants_data_matlab{k}.angular_frame_initial = angular_frame_initial;
end

%% Plot the invariants and SE3 data
close all;
haxes1 = figure('NumberTitle', 'off', 'Name', 'DHB Invariants');
dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};

haxes2 = figure('NumberTitle', 'off', 'Name', 'SE3 Paths');
coordNames = {'x' 'y' 'z' 'rx' 'ry' 'rz'};
use_invariants_from_matlab = true;

method_plot = 1;

% data_skip_list = [5, 8, 11, 15]; % pouring
data_skip_list = []; % fruits

for k=1:numTests
    if ismember(k, data_skip_list)
        continue;
    end

    fprintf('Plotting for the data %d ....\n', k);

    % invariants
    invariants_python = dhb_invariants.(demo_list{k});

    if (use_invariants_from_matlab)
        % use invarians from matlab
        invariants_pos = invariants_data_matlab{k}.invariants_pos;
    else
        invariants_pos = [invariants.linear_motion_invariants, ...
            invariants.angular_motion_invariants];
    end

    set(0,'CurrentFigure',haxes1)
    if method_plot == 1
        for i=1:6
            subplot(2,3,i)
            plot(invariants_pos(:,i),'LineWidth',2)
            hold on;
            ylabel(dhbInvNames{i});
            grid on
        end
    elseif method_plot == 2


    end

    % SE3 coordinates
    set(0,'CurrentFigure',haxes2)
    for i=1:6
        subplot(2,3,i)
        if (i < 4)
            plot(path_data{k}.pos_data(:,i),'LineWidth',2)
        else
            plot(path_data{k}.rvec_data(:,i-3),'LineWidth',2)
        end
        hold on;
        ylabel(coordNames{i});
        grid on;
        drawnow;
    end

    keyboard;
end

%% Reconstructed data from the invariants

haxes3 = figure('NumberTitle', 'off', 'Name', 'SE3 Paths (Reconstructed)');
coordNames = {'x' 'y' 'z' 'rx' 'ry' 'rz'};
for k=1:numTests
    if ismember(k, data_skip_list)
        continue;
    end
    invariants_pos = invariants_data_matlab{k}.invariants_pos;
    linear_frame_initial = invariants_data_matlab{k}.linear_frame_initial;
    angular_frame_initial = invariants_data_matlab{k}.angular_frame_initial;
    [pos_r_data, rvec_r_data] = reconstructTrajectory(invariants_pos, linear_frame_initial, angular_frame_initial, 'pos');

    % SE3 coordinates
    set(0,'CurrentFigure',haxes3)
    for i=1:6
        subplot(2,3,i)
        if (i < 4)
            plot(pos_r_data(:,i),'LineWidth',2)
        else
            plot(rvec_r_data(:,i-3),'LineWidth',2)
        end
        hold on;
        ylabel(coordNames{i});
        grid on;
        drawnow;
    end
end

end

% Construct a pose with an offset in translation and rotation
function poseNew = constructPoseWithOffset(pose, offset)
rot_offset = deg2rad(offset(4:6));
poseNew = SpatialRobotModel.transl(offset(1), offset(2), offset(3)) * pose * ...
    SpatialRobotModel.r2t(SpatialRobotModel.rpy2R(rot_offset(1), rot_offset(2), rot_offset(3)));
end

% a function that converts rotation vector data to quaternion vector data
function qt_r_data = rvecdata2qt(rvec_r_data)
    N = size(rvec_r_data, 1);
    qt_r_data = zeros(N, 4);
    for i = 1:N
        rvec_r = rvec_r_data(i,:);
        rotm_r = rotationVectorToMatrix(rvec_r);
        qt_r = rotm2quat(rotm_r);
        qt_r_data(i, :) = qt_r;
    end
end

% an example cost function for trajectory optimization
% for shape preservation in terms of the shape invariant descriptor
% with normalized cost per parameter distribution
function cost = cost_shape_descriptor(pos_data, rvec_data, T0, pos_invariant_orig)

    % map the current position path into the invariants
    pos_diff_data = diff(pos_data');
    [pos_invariant, rot_invariant, pos_initial, rot_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);

    % compute reconstruction errors
    error = abs(pos_invariant - pos_invariant_orig);
    error_normalized = normalize(error, "range");
    cost = norm(error_normalized);
end

% an example cost function for trajectory optimization
% for shape preservation in terms of the shape invariant descriptor
% with normalized cost per parameter distribution
% with mexified computeDHB method
function cost = cost_shape_descriptor_mex(pos_data, rvec_data, T0, pos_invariant_orig, weights)

    % map the current position path into the invariants
    pos_diff_data = diff(pos_data');
    [pos_invariant, rot_invariant, pos_initial, rot_initial] = computeDHBMex_mex(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);

    % compute reconstruction errors
    error = abs(pos_invariant - pos_invariant_orig);
    error_normalized = normalize(error, "range");

    % apply some weight for importance
    cost = norm(weights .* error_normalized);
end

% an example cost function for trajectory optimization
% for shape preservation in terms of the shape invariant descriptor
% with normalized cost per parameter distribution
% with mexified computeDHB method
% with additional cost of the length
% with additional weights on the initial and final direction
function cost = cost_shape_descriptor_mex2(pos_data, rvec_data, T0, pos_invariant_orig, weights)

    % map the current position path into the invariants
    pos_diff_data = diff(pos_data');
    [pos_invariant, rot_invariant, pos_initial, rot_initial] = computeDHBMex_mex(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);

    % compute reconstruction errors
    error = abs(pos_invariant - pos_invariant_orig);
    error_normalized = normalize(error, "range");

    % compute the length of the path
    cost_path_length = weights(4) * norm(pos_invariant(:,1));

    % construct the cost
    cost_error = norm(weights(1:3) .* error_normalized);

    % apply some weight for importance
    cost = cost_error + cost_path_length;
end

function se3Path = generateRandomComplexSE3Path(N, numWaypoints, positionRange)
    % Generates a complex and smooth SE3 path with random waypoints.
    % N - Number of poses in the path
    % numWaypoints - Number of waypoints to generate
    % positionRange - Range for position coordinates, e.g., [-10 10]

    % Randomly generate position waypoints
    waypoints = positionRange(1) + (positionRange(2) - positionRange(1)) .* rand(numWaypoints, 3);

    % Randomly generate orientation waypoints
    orientationWaypoints = arrayfun(@(x) generateRandomUnitQuaternion(), 1:numWaypoints, 'UniformOutput', false);

    % Interpolate the path
    se3Path = interpolateSE3Path(N, waypoints, orientationWaypoints);
end

function se3Path = interpolateSE3Path(N, waypoints, orientationWaypoints)
    % Time vector for waypoints and path
    tWaypoints = linspace(0, 1, size(waypoints, 1));
    tPath = linspace(0, 1, N);

    % Spline interpolation for positions
    positionPath = interp1(tWaypoints, waypoints, tPath, 'spline');

    % SLERP interpolation for orientations
    quatPath = slerpPath(orientationWaypoints, tPath, tWaypoints);

    % Combine position and quaternion paths
    se3Path = [positionPath, cell2mat(quatPath')];
end

function quat = generateRandomUnitQuaternion()
    % Generates a random unit quaternion
    quat = [rand; randn(3, 1)];
    quat = quat / norm(quat);
end

function quatPath = slerpPath(controlQuats, tPath, tControl)
    % Interpolates a path of quaternions using slerp
    quatPath = cell(1, length(tPath));
    for i = 1:length(tPath)
        % Find the interval in the control points
        t = tPath(i);
        idx = find(tControl <= t, 1, 'last');
        if idx == length(controlQuats)
            qt = controlQuats{end};
            quatPath{i} = [qt(1:3)', qt(4)];
        else
            t0 = tControl(idx);
            t1 = tControl(idx + 1);
            frac = (t - t0) / (t1 - t0);
            qt = SpatialRobotModel.quatInterp(controlQuats{idx}, controlQuats{idx + 1}, frac);
            quatPath{i} = [qt(1:3)', qt(4)];
        end
    end
end