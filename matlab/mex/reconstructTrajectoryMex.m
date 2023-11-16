%#codegen
%% Functions that reconstruct a Cartesian trajectory from its DHB invariant representation
% This script reconstructs the trajectory of a rigid body in 3D space from
% its Denavit–Hartenberg inspired Bidirectional (DHB) invariant representation.

%% Reconstruct Cartesian trajectory (position or velocity)
% Inputs:
%        dhb_invariants - Matrix of DHB invariants ([N-2]x6 array) 
%        initial_linear_frame - Initial linear frame (4x4 matrix)
%        initial_angular_frame - Initial angular frame (4x4 matrix)
%        method - String, 'pos' for position or otherwise for velocity
%
% Outputs:
%        trajectory_position - Position or linear velocity trajectory ([N-2]x3 array)
%        trajectory_rotation - Relative rotation vector or angular velocity trajectory ([N-2]x3 array)

function [trajectory_position, trajectory_rotation] = reconstructTrajectoryMex(dhb_invariants, initial_linear_frame, initial_angular_frame, method)
    linear_magnitude = dhb_invariants(:,1);
    linear_angle_y = dhb_invariants(:,2);
    linear_angle_x = dhb_invariants(:,3);
    angular_magnitude = dhb_invariants(:,4);
    angular_angle_y = dhb_invariants(:,5);
    angular_angle_x = dhb_invariants(:,6);

    num_samples = size(linear_angle_y,1);

    trajectory_position = zeros(num_samples,3);
    trajectory_rotation = zeros(num_samples,3);

    position_mode = strcmp(method, 'pos');

    initial_linear_frame(4,4) = double(position_mode);

    for i = 1:num_samples
        if position_mode
            trajectory_position(i,:) = initial_linear_frame(1:3,4)';
        end
        
        % Compute rotation matrices for the position or velocity
        rotation_matrix_position = rotateY(linear_angle_y(i)) * rotateX(linear_angle_x(i));
        translation_vector = [linear_magnitude(i) 0 0]';
        transformation_matrix = [rotation_matrix_position translation_vector; 0 0 0 double(position_mode)];
        initial_linear_frame = initial_linear_frame * transformation_matrix;
        
        % Only update the position if we're dealing with velocities
        if ~position_mode
            trajectory_position(i,:) = initial_linear_frame(1:3,4)';
        end

        % Apply the rotation to the angular frame and update the rotation trajectory
        trajectory_rotation(i,:) = (initial_angular_frame * [angular_magnitude(i); 0; 0])';
        rotation_matrix_rotation = rotateY(angular_angle_y(i)) * rotateX(angular_angle_x(i));
        initial_angular_frame = initial_angular_frame * rotation_matrix_rotation;
    end 
end

%% Elementary rotation around the x-axis
% phi - Angle to rotate about the x-axis
function rotation_matrix = rotateX(phi)
    rotation_matrix = [1        0         0; ...
                       0 cos(phi) -sin(phi); ...
                       0 sin(phi)  cos(phi)];
end

%% Elementary rotation around the y-axis
% beta - Angle to rotate about the y-axis
function rotation_matrix = rotateY(beta)
    rotation_matrix = [cos(beta) 0 sin(beta); ...
                       0         1         0; ...
                      -sin(beta) 0 cos(beta)];
end
