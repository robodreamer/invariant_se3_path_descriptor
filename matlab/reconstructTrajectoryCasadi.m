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

function [trajectory_position, trajectory_rotation] = reconstructTrajectoryCasadi(dhb_invariants, initial_linear_frame, initial_angular_frame, method)

    num_samples = size(dhb_invariants,1);

    % Use CasADi SX or MX for trajectory arrays
    trajectory_position = casadi.SX.zeros(num_samples,3);
    trajectory_rotation = casadi.SX.zeros(num_samples,3);
    % trajectory_position = zeros(num_samples,3);
    % trajectory_rotation = zeros(num_samples,3);

    % Initialize frames as CasADi symbolic expressions
    linear_frame = initial_linear_frame;
    angular_frame = initial_angular_frame;

    for i = 1:num_samples
        % Call the single step reconstruction for each row of invariants
        [linear_frame, angular_frame, new_position, new_rotation] = reconstructTrajectorySingleStep(dhb_invariants(i,:), linear_frame, angular_frame, method);

        % Store the position and rotation in the trajectory arrays
        trajectory_position(i,:) = new_position;
        trajectory_rotation(i,:) = new_rotation;
    end
end