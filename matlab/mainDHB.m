% Author: Matteo Saveriano - 01.03.18
% Modified by Andy Park in Nov, 2023

clear
close all
clc

%% Construct 6DOF trajectory
dt = 0.1; % sample time
tf = 4;   % total time
N = 1000; % number of samples

t_data = linspace(0,tf,N)';

x = 0.1*exp(t_data);
y = 5+1.5*sin(t_data);
z = cos(t_data);

roll  = 0.5*sin(t_data);
pitch = cos(t_data);
yaw   = 0.1*t_data;

%% construct the input to the DHB invariant calcuation

% generate position and rotation inputs
pos_data = [x, y, z];
rot_data = [roll, pitch, yaw];

% Initial position
x0 = pos_data(1,1);
y0 = pos_data(1,2);
z0 = pos_data(1,3);

% Initial orientation in radians
roll0 = rot_data(1,1);
pitch0 = rot_data(1,2);
yaw0 = rot_data(1,3);

% Construct the rotation matrix from roll, pitch, yaw
Rz = [cos(yaw0) -sin(yaw0) 0; sin(yaw0) cos(yaw0) 0; 0 0 1]; % Rotation around Z axis (yaw)
Ry = [cos(pitch0) 0 sin(pitch0); 0 1 0; -sin(pitch0) 0 cos(pitch0)]; % Rotation around Y axis (pitch)
Rx = [1 0 0; 0 cos(roll0) -sin(roll0); 0 sin(roll0) cos(roll0)]; % Rotation around X axis (roll)

% Combined rotation matrix
R = Rz * Ry * Rx;

% Homogeneous transformation matrix
T0 = [R [x0; y0; z0]; 0 0 0 1];

% rvec_data
rvec_data = zeros(size(rot_data));

% Loop through each set of Euler angles in 'rot_data'
for i = 1:size(rot_data, 1)
    % Extract the Euler angles for the i-th sample
    eul = rot_data(i, :);

    % Convert Euler angles to rotation matrix
    rotm = eul2rotm(eul, 'ZYX'); % 'ZYX' is one of the most common sequences

    % Convert the rotation matrix to axis-angle representation
    axang = rotm2axang(rotm);

    % Convert axis-angle to rotation vector and store it
    rvec_data(i, :) = axang(4) * axang(1:3);
end

% Now 'rvec_data' contains the rotation vectors for all sets of Euler angles

% Compute velocity (twist)
twists(:,4:6) = diff([x y z],1,1)/dt;

orientRate = diff([roll pitch yaw],1,1)/dt;
for i=1:N-1
    Tr = [ 1           0           -sin(pitch(i));...
          0   cos(roll(i))  cos(pitch(i))*sin(roll(i));...
          0  -sin(roll(i))  cos(pitch(i))*cos(roll(i)) ];

    twists(i,1:3) = (Tr * orientRate(i,:)')';

    a(i) = norm(twists(i,1:3));
end


%% Compute DHB invariants

% Compute position based DHB invariants
pos_diff_data = diff(pos_data);
[linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = computeDHB(pos_diff_data, rvec_data(1:end-1,:), 'pos', T0);
invariants_pos = [linear_motion_invariant, angular_motion_invariant];

% Compute velocity based DHB invariants
[linear_motion_invariant_v, angular_motion_invariant_v, linear_frame_initial_v, angular_frame_initial_v] = computeDHB(twists(:,4:6), twists(:,1:3), 'vel');
invariants_vel = [linear_motion_invariant_v, angular_motion_invariant_v];

%% Plot invariant trajectories

figure('NumberTitle', 'off', 'Name', 'Cartesian pose to DHB');
dhbInvNames = {'m_p' '\theta_p^1' '\theta_p^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};
for i=1:6
    subplot(2,3,i)
    plot(t_data(1:end-3),invariants_pos(:,i),'LineWidth',2)
    ylabel(dhbInvNames{i});
    grid on
end

figure('NumberTitle', 'off', 'Name', 'Cartesian velocity to DHB');
dhbInvNames = {'m_v' '\theta_v^1' '\theta_v^2' 'm_{\omega}' '\theta_{\omega}^1' '\theta_{\omega}^1'};
for i=1:6
    subplot(2,3,i)
    plot(t_data(1:end-3),invariants_vel(:,i),'LineWidth',2)
    ylabel(dhbInvNames{i});
    grid on
end

%% Reconstruct original trajectory
% position
[pr, rvec_r] = reconstructTrajectory(invariants_pos, linear_frame_initial, angular_frame_initial, 'pos');

% velocity
[vr, wr] = reconstructTrajectory(invariants_vel, linear_frame_initial_v, angular_frame_initial_v, 'vel');

%% Compute reconstruction errors

%-- error with position
errSP_pos = [(pr - pos_data(1:N-3,:)).^2 (rvec_r - rvec_data(1:N-3,:)).^2];

% Compute rmse error
err_pos = zeros(1,6);
for i=1:6
    err_pos(i) = sum(errSP_pos(:,i));
end

RMSE_pos = sqrt([sum(err_pos(1:3)) sum(err_pos(4:6))]./(N-2));
disp(['Reconstruction errors in pose (RMSE): ' num2str(RMSE_pos)])

%-- error with velocity
errSP_vel = [(vr - twists(1:N-3,4:6)).^2 (wr - twists(1:N-3,1:3)).^2];

% Compute rmse error
err_vel = zeros(1,6);
for i=1:6
    err_vel(i) = sum(errSP_vel(:,i));
end

RMSE_vel = sqrt([sum(err_vel(1:3)) sum(err_vel(4:6))]./(N-3));
disp(['Reconstruction errors in velocity (RMSE): ' num2str(RMSE_vel)])

%% Plot original and reconstructed paths

% Plot original and reconstructed paths -- position
pose_data = [rvec_data, pos_data]; % original
pose_r_data = [rvec_r, pr]; % reconstructed

figure('NumberTitle', 'off', 'Name', 'DHB to Cartesian pose');
for i=1:6
    subplot(2,3,i)
    plot(t_data, pose_data(:,i),'g','LineWidth',4)
    hold on;
    plot(t_data(1:end-3),pose_r_data(:,i),'b','LineWidth',2)
    if(i<4)
        ylabel(['r_' num2str(i)])
    else
        ylabel(['p_' num2str(i-3)])
    end
    grid on
end

% Plot original and reconstructed paths -- velocity
figure('NumberTitle', 'off', 'Name', 'DHB to Cartesian velocity');
for i=1:6
    subplot(2,3,i)
    plot(t_data(1:end-1),twists(:,i),'g','LineWidth',4)
    hold on;
    if(i<4)
        plot(t_data(1:end-3),wr(:,i),'b','LineWidth',2)
        ylabel(['\omega_' num2str(i)])
    else
        plot(t_data(1:end-3),vr(:,i-3),'b','LineWidth',2)
        ylabel(['x_' num2str(i-3)])
    end
    grid on
end
