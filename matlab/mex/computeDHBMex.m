%#codegen
%% Functions that implement the Denavit–Hartenberg inspired Bidirectional (DHB) invariant representation
% See: D. Lee, R. Soloperto, and M. Saveriano, "Bidirectional invariant
%      representation of rigid body motions and its application to gesture
%      recognition and reproduction", Auton. Robots, 42(1):125–145, 2018.

%% Compute DHB invariants (position- or velocity-based)
% Input: 
%        position_diff - position difference or linear velocity (Nx3 array) 
%        rotation_diff - rotation vector difference or angular velocity (Nx3 array)
%        method - 'pos': position-based DHB, 'vel': velocity-based DHB
%        initial_pose - initial pose (used only if method='pos')
% 
% Output:
%        linear_motion_invariant - first linear (position or velocity) invariant (N-2 array)
%        angular_motion_invariant - first angular (position or velocity) invariant (N-2 array)
%        linear_frame_initial - Initial linear frame
%        angular_frame_initial - Initial angular frame


function [linear_motion_invariant, angular_motion_invariant, linear_frame_initial, angular_frame_initial] = computeDHB(position_diff, rotation_diff, method, initial_pose)

[num_samples, ~] = size(rotation_diff);

% Compute initial frames for linear and angular components
linear_frame_x = computeFrameAxisX(position_diff(1,:), [1 0 0]);
linear_frame_x2 = computeFrameAxisX(position_diff(2,:), linear_frame_x);
linear_frame_y = computeFrameAxisY(linear_frame_x, linear_frame_x2, ...
    [linear_frame_x(2)-linear_frame_x(3), ...
        linear_frame_x(3)-linear_frame_x(1), ...
        linear_frame_x(1)-linear_frame_x(2)]/norm([linear_frame_x(2)-linear_frame_x(3), ...
        linear_frame_x(3)-linear_frame_x(1), linear_frame_x(1)-linear_frame_x(2)]));
linear_frame_z = cross(linear_frame_x,linear_frame_y);
linear_frame_z = linear_frame_z / norm(linear_frame_z);

% Construct the initial linear frame
linear_frame_initial = eye(4);
linear_frame_initial(1:3,1:3) = ([linear_frame_x', linear_frame_y', linear_frame_z']);
if(strcmp(method, 'pos'))
    linear_frame_initial(1:3,4) = initial_pose(1:3,4);
else
    linear_frame_initial(1:3,4) = position_diff(1,:);
end

% Similar construction for the initial angular frame
angular_frame_x = computeFrameAxisX(rotation_diff(1,:), [1 0 0]);
angular_frame_x2 = computeFrameAxisX(rotation_diff(2,:), angular_frame_x);
angular_frame_y = computeFrameAxisY(angular_frame_x, angular_frame_x2, [0 1 0]);
angular_frame_z = cross(angular_frame_x,angular_frame_y);
angular_frame_z = angular_frame_z / norm(angular_frame_z);

angular_frame_initial = eye(3);

angular_frame_initial(1:3,1:3) = ([angular_frame_x', angular_frame_y', angular_frame_z']);

% Initialize arrays to hold the invariant values
linear_motion_invariant = zeros(num_samples-2,3);
angular_motion_invariant = zeros(num_samples-2,3);

% Compute invariant values for each sample
for i = 1:num_samples-2  
    % Update frames based on new position or rotation
    linear_frame_x3 = computeFrameAxisX(position_diff(i+2,:), linear_frame_x2);
    linear_frame_y2 = computeFrameAxisY(linear_frame_x2, linear_frame_x3, linear_frame_y);
    linear_motion_invariant(i,:) = computeInvariants(position_diff(i,:), linear_frame_x, linear_frame_x2, linear_frame_y, linear_frame_y2);
    
    linear_frame_x = linear_frame_x2;
    linear_frame_x2 = linear_frame_x3;
    linear_frame_y = linear_frame_y2;

    angular_frame_x3 = computeFrameAxisX(rotation_diff(i+2,:), angular_frame_x2);
    angular_frame_y2 = computeFrameAxisY(angular_frame_x2, angular_frame_x3, angular_frame_y);
    angular_motion_invariant(i,:) = computeInvariants(rotation_diff(i,:), angular_frame_x, angular_frame_x2, angular_frame_y, angular_frame_y2);
    
    angular_frame_x = angular_frame_x2;
    angular_frame_x2 = angular_frame_x3;
    angular_frame_y = angular_frame_y2;
end

end

%% Compute the three invariant values given linear or angular frame axes (x,y)
function invariants = computeInvariants(vector_u, frame_x, frame_x2, frame_y, frame_y2)
    magnitude = dot(frame_x, vector_u);
    angle1 = atan2(dot(cross(frame_x, frame_x2), frame_y), dot(frame_x, frame_x2));
    angle2 = atan2(dot(cross(frame_y, frame_y2), frame_x2), dot(frame_y, frame_y2));
    invariants = [magnitude, angle1, angle2];
end

%% Compute the x axis of linear or angular frame
function frame_x = computeFrameAxisX(vector_u, default_x)
    frame_x = vector_u;
    norm_frame_x = norm(frame_x);
    if norm_frame_x > 1e-10
        frame_x = frame_x / norm_frame_x;
    else
        frame_x = default_x;
    end
end

%% Compute the y axis of linear or angular frame
function frame_y = computeFrameAxisY(frame_x, frame_x2, default_y)
    frame_y = cross(frame_x, frame_x2);
    frame_y = normalizeVector(frame_y, default_y);
end

%% Normalize a vector, replace with default if below threshold
function normalized_vector = normalizeVector(vector, default_vector)
    norm_vector = norm(vector);
    if(norm_vector > 1e-10)
        normalized_vector = vector / norm_vector;
    else
        normalized_vector = default_vector;
    end
end

