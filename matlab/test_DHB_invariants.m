%%% Testing the concept of DHB invariants
% Andy Park, Nov 2023

addpath(genpath(pwd));
ccc;

select_program = 2;

%%% 11-10-23
% 1) show the invariance over affine transforms
% 2) visualize the reconstructed path with different initial poses
% 3) fit a spline
if (select_program == 1)

    %% Data Preparation

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
    [linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
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
if (select_program == 2)

    %% Data Preparation

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
    toc;

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
function [cost, cost_array] = cost_shape_descriptor(pos_data, rvec_data, T0, pos_invariant_orig)

    % map the current position path into the invariants
    pos_diff_data = diff(pos_data');
    [pos_invariant, rot_invariant, pos_initial, rot_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);

    % compute reconstruction errors
    cost = norm(abs(pos_invariant - pos_invariant_orig));

    % disp(['Reconstruction errors in pose (RMSE): ' num2str(cost)])
end