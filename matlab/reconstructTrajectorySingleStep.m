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

function [new_linear_frame, new_angular_frame, new_position, new_rotation] = reconstructTrajectorySingleStep(dhb_invariant, current_linear_frame, current_angular_frame, method)
    % Extract individual components from the invariant
    linear_magnitude = dhb_invariant(1);
    linear_angle_y = dhb_invariant(2);
    linear_angle_x = dhb_invariant(3);
    angular_magnitude = dhb_invariant(4);
    angular_angle_y = dhb_invariant(5);
    angular_angle_x = dhb_invariant(6);

    position_mode = strcmp(method, 'pos');

    new_position = casadi.SX.zeros(1,3);
    new_rotation = casadi.SX.zeros(1,3);

    if (position_mode)
        new_position = current_linear_frame(1:3,4)';
    end

    % Compute rotation matrices for the position or velocity
    rotation_matrix_position = rotateY(linear_angle_y) * rotateX(linear_angle_x);
    translation_vector = [linear_magnitude; 0; 0];
    transformation_matrix = [rotation_matrix_position, translation_vector; 0, 0, 0, position_mode];
    new_linear_frame = current_linear_frame * transformation_matrix;

    % Determine the new position
    if (~position_mode)
        new_position = new_linear_frame(1:3,4)';
    end

    % Apply the rotation to the angular frame and update the rotation trajectory
    new_rotation = (current_angular_frame * [angular_magnitude; 0; 0])';
    rotation_matrix_rotation = rotateY(angular_angle_y) * rotateX(angular_angle_x);
    new_angular_frame = current_angular_frame * rotation_matrix_rotation;
end

% Elementary rotation around the x-axis
function rotation_matrix = rotateX(phi)
    rotation_matrix = [1        0         0; ...
                       0 cos(phi) -sin(phi); ...
                       0 sin(phi)  cos(phi)];
end

% Elementary rotation around the y-axis
function rotation_matrix = rotateY(beta)
    rotation_matrix = [cos(beta) 0 sin(beta); ...
                       0         1         0; ...
                      -sin(beta) 0 cos(beta)];
end
