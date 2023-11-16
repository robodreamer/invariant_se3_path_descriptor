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
                                 
% Compute DHB invariants

% Compute position based DHB invariants
[linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = computeDHBMex(pos_data, rvec_data, 'vel', T0);
invariants_pos = [linear_motion_invariant, angular_motion_invariant];


% Reconstruct original trajectory
% position
[pr, rvec_r] = reconstructTrajectory(invariants_pos, linear_frame_initial, angular_frame_initial, 'vel');

% Compute reconstruction errors

%-- error with position
errSP_pos = [(pr - pos_data(1:N-2,:)).^2 (rvec_r - rvec_data(1:N-2,:)).^2];

% Compute rmse error
err_pos = zeros(1,6);
for i=1:6
    err_pos(i) = sum(errSP_pos(:,i));
end

RMSE_pos = sqrt([sum(err_pos(1:3)) sum(err_pos(4:6))]./(N-2));
disp(['Reconstruction errors in pose (RMSE): ' num2str(RMSE_pos)])


