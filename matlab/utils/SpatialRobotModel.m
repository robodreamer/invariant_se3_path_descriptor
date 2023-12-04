classdef SpatialRobotModel < handle
    % SpatialRobotModel Robotics Library Class
    %   This is a superclass for robot models in AP Toolbox. All the
    %   generic computations are defined in this class. Spatial Vector is
    %   used in kinematics and dynamics computations.
    %
    %   @Author: Andy Park (andy.park@kuka.com)

    properties (Abstract = true, SetAccess=protected)
        name; % robot name
    end

    properties (SetAccess = protected)

        %-- numbers of bodies, joints, frames

        NB; % number of bodies:
        N_joint; % number of joints (including floating-body)
        N_ext; % total number of bodies + dummy bodies (end-effector frames)
        NJ; % number of actual joints

        URDF_name; % URDF filename
        URDF_par; % URDF object
        baseName; % base link name
        linkNames; % link names
        frameNames; % frame names (end-effectors)
        collisionShapeFile; % collision shape prestored file

        jointNames; % joint names

        revoluteJointNames = {}; % revolute joint list

        %-- visualiation related

        Nmeshes; % no. of mesh files
        fv; % patch variables (faces + vertices) for meshes
        h_fv; % handles for patches
        T_fv; % transform for patches
        h_frame; % handle for coordinate frames
        vdata; % data for vertices

        fv_com; % patches for CoM sphere meshes
        h_fv_com; % handles for patches of CoM
        p_fv_com; % position vectors for patches of CoM

        p_mesh_offset; % position offsets for mesh visualization
        R_mesh_offset; % rotation offsets for mesh visualization

        fv_col; % collision shape patches
%         h_fv_col; % handles for patches
        p_fv_col; % position vectors for patches of collision shapes
        T_fv_col; % transform for patches
        collision; % collision shape related field

        param_backup; % model parameters backup

        %-- joint limits

        joint_pos_min; % joint limit -- minimum position
        joint_pos_max; % joint limit -- maximum position
        joint_vel_max; % joint limit -- maximum velocity
        joint_acc_max; % joint limit -- maximum acceleration
        joint_tq_max; % joint limit -- maximum torque
    end

    properties (Constant = true)
        ag_0 = [0 0 0 0 0 -9.80665]'; % gravity vector
    end

    properties (Access = public)
        %-- base motion

        base_motion = true; % enable/disable setting a custom pose for base
        R_base = eye(3); % base orientation (zero by default)

        %-- scale factors for model kinematics and inertial parameters

        scale_length = 1e0; % scale factor for length
        scale_inertia = 1e0; % scale factor for inertia tensors
        scale_ratio = 1e0; % scale ratio for 3d mesh patches

        %-- properties related to model kinematics/dynamics compuatation

        parent; % parent links
        phi; % link axis of rotation
        i_p_com; % link CoM in link coordinate frame
        m; % link mass
        I_com; % link inertia tensor
        pi_p_i; % link coordinate frame position w.r.t. joint frame
        pi_R_i; % link coordinate frame rotation w.r.t. joint frame
        joints; % corresponding joint indices for the current link
        X_J; % 6x6 transforms for joint motion
        i_XL_pi; % 6x6 transforms for link coordinate frame w.r.t. joint coordinate frame (excluding joint motion)
        i_X_pi; % 6x6 transforms for link coordinate frame w.r.t. parent coordinate frame (after joint motion)
        T0; % 4x4 transforms for link coordinate frames
        R0; % 3x4 rotation matrices for link coordinate frames (i.e., 3x3 portion of T0)
        I; % Inertia matrices
        I_C; % 6x6 composite inertia matrices
        i_X_0; % 6x6 transforms for link coordinate frames
        pi_T_i; % 4x4 transform between link coordinate frames and joint coordinate frames
        R0_6; % 6x6 rotation matrices for link coordinate frames
        v; % spatial velocity vectors
        a; % spatial acceleration vectors
        chi; % velocity product vectors
        f; % force vectors for each joint
        a0; % acceleration vector for the base
        q; % joint position variable
        q_a; % joint position for actuated joints
        q_p; % joint position for passive joints (6DOF for base)
        bodyname; % names of bodies
        f_ext_zero; % external force set to zero by default
        type; % type of the joint

        % inertia index ordering indices
        % These default indices [1 4 6 2 5 3] indicate that
        % [ixx iyy izz ixy iyz ixz] is what SpatialRobotModel needs,
        % but URDF provides them as [ixx, ixy, ixz, iyy, iyz, izz]
        inertia_index = [1 4 6 2 5 3];

        model_param; % struct of model parameters used in mex functions

        %-- visualization related

        frame_len = 3e-2; % coordinate frame line lengths
        frame_width = 2; % coordinate frame line widths

        display_animation = false; % animation flag
        display_mass = 0; % display CoM
        display_frames = 1; % display frames
        display_inertia = 1; % display inertia and CoM
        mesh_alpha = 0.6; % mesh transparency
        mesh_color = {}; % color of the body links
        display_collision_shapes = 0; % display collision shapes (e.g., 1: cube, 2: capsule)
        load_stored_collision_shapes = 0; % load stored collision shapes (capsules)
        collision_shape_color = [103, 156, 239]/255; % collision shape colors
        collision_shape_alpha = 0.3; % collision shape alpha
%         collision_capsules_param;

        color_segment = [115 159 238]/255; % color of the segment
        width_segment = 2; % width of segment
        color_inertia = [36 97 255]/255; % color for inertia ellipsoid
        width_inertia = 1; % width of the inertia ellipsoids
        frame_scale = 0.03; % scale factor for frames

        mesh_name; % link mesh file names

        reduce_patch_ratio = 0.2 % patch reduced ratio for faster rendering of 3D mesh files

        %-- custom coordinate frame indices

        ind_H; % end-effector frame

        %-- miscellaneous

        use_mex = true; % use mex functions (set to true by default)
        T_im; % transform matrix used for interactive marker

        %-- for collision shapes
        h_fv_col; % handles for patches

    end

    %%% useful static functions %%%

    methods (Static = true)

        %% interpolation related functions

        function dataNew = cubicInterp(data, time, timeNew, method)
            % CUBICINTERP cubic interpolation for a N-dim data
            %   dataNew = CUBICINTERP(data, time, timeNew) interpolates
            %   data with timeNew and generates dataNew.
            %
            %   Example:
            %   time = 0:10:360;
            %   data = sind(time)';
            %   timeNew = 0:5:360;
            %   out = SpatialRobotModel.CUBICINTERP(data,time,timeNew)
            
            if ~exist('method','var')
                method = 'pchip';
            end

            dataNew = zeros(length(timeNew),size(data,2));
            for i=1:size(data,2)
                % one dim. data
                data_tmp = data(:,i);
                % apply cubic interpolation to get a new data
                dataNew(:,i) = interp1(time, data_tmp, timeNew, method);
            end
        end

        function [qt,qdt,qddt] = jtraj(q0, q1, tv, qd0, qd1)
            % JTRAJ compute interpolated trajectory between two configurations
            %   [q,qd,qdd] = JTRAJ(q0, qf, m)
            %
            %   This function returns a joint space trajectory q (MxN)
            %   where the joint coordinates vary from q0 (1xN) to qf (1xN).
            %   A quintic (5th order) polynomial is used with default zero
            %   boundary conditions for velocity and acceleration. Time is
            %   assumed to vary from 0 to 1 in m steps. Joint velocity and
            %   acceleration can be optionally returned as qd (MxN) and qdd
            %   (MxN) respectively. The trajectory q, qd and qdd are MxN
            %   matrices, with one row per time step, and one column per
            %   joint.
            %
            %   [q,qd,qdd] = JTRAJ(q0, qf, m, qd0, qdf) as above but also
            %   specifies initial and final joint velocity for the
            %   trajectory.
            %
            %   [q,qd,qdd] = JTRAJ(q0, qf, T) as above but the trajectory
            %   length is defined by the length of the time vector T (Mx1).
            %
            %   [q,qd,qdd] = JTRAJ(q0, qf, T, qd0, qdf) as above but
            %   specifies initial and final joint velocity for the
            %   trajectory and a time vector.
            %
            %   Copyright (C) 1993-2008, by Peter I. Corke

            if length(tv) > 1
                tscal = max(tv);
                t = tv(:)/tscal;
            else
                tscal = 1;
                t = [0:(tv-1)]'/(tv-1);	% normalized time from 0 -> 1
            end

            q0 = q0(:);
            q1 = q1(:);

            if nargin == 3,
                qd0 = zeros(size(q0));
                qd1 = qd0;
            elseif nargin == 5,
                qd0 = qd0(:);
                qd1 = qd1(:);
            else
                error('incorrect number of arguments')
            end

            % compute the polynomial coefficients
            A = 6*(q1 - q0) - 3*(qd1+qd0)*tscal;
            B = -15*(q1 - q0) + (8*qd0 + 7*qd1)*tscal;
            C = 10*(q1 - q0) - (6*qd0 + 4*qd1)*tscal;
            E = qd0*tscal; % as the t vector has been normalized
            F = q0;

            tt = [t.^5 t.^4 t.^3 t.^2 t ones(size(t))];
            c = [A B C zeros(size(A)) E F]';

            qt = tt*c;

            % compute optional velocity
            if nargout >= 2,
                c = [ zeros(size(A)) 5*A 4*B 3*C  zeros(size(A)) E ]';
                qdt = tt*c/tscal;
            end

            % compute optional acceleration
            if nargout == 3,
                c = [ zeros(size(A))  zeros(size(A)) 20*A 12*B 6*C  zeros(size(A))]';
                qddt = tt*c/tscal^2;
            end
        end

        function d = T2diff(t1,t2)
            % T2DIFF Convert homogeneous transform to differential motion
            %
            %   d = T2DIFF(T0, T1) is the differential motion (6 � 1)
            %   corresponding to infinitessimal motion from pose T0 to T1
            %   which are homogeneous transformations (4�4).
            %   d = (dRx, dRy, dRz, dx, dy, dz) and is an approximation to
            %   the average spatial velocity multiplied by time.
            %
            %   d = T2DIFF(T) is the differential motion corresponding to
            %   the infinitessimal relative pose T expressed as a
            %   homogeneous transformation.
            %
            %   Copyright (C) 1993-2008, by Peter I. Corke

            if nargin == 1,
                d = [0.5*[t1(3,2)-t1(2,3); t1(1,3)-t1(3,1); t1(2,1)-t1(1,2)];t1(1:3,4)];
            else

                d = [0.5*(cross(t1(1:3,1), t2(1:3,1)) + ...
                    cross(t1(1:3,2), t2(1:3,2)) + ...
                    cross(t1(1:3,3), t2(1:3,3)) ...
                    ); t2(1:3,4)-t1(1:3,4)];
            end
        end

        %% transform related functions

        function out = passive(in)
            % PASSIVE take the passive part of the vector/matrix
            %   q_p = PASSIVE(q) is the first 6 dim vector in q that
            %   corresponds to the passive DOFs.

            if(size(in,2) == 1) % only one column (e.g. q,qd,qdd)
                out = in(1:6,1);
            else % if no. column > 1 (e.g. J)
                out = in(:,1:6);
            end
        end

        function out = active(in)
            % ACTIVE take the actuated part of the vector/matrix
            %   q_a = ACTIVE(q) is the q(7:end) vector that corresponds to
            %   the active DOFs in q. This works for matrices as well
            %   (i.e., J_a = ACTIVE(J) where J_a is Jacobian for active
            %   DOF components).

            if(size(in,2) == 1) % only one column (e.g. q,qd,qdd)
                out = in(7:end,1);
            else % if no. column > 1 (e.g. J)
                out = in(:,7:end);
            end
        end

        function out = takePos(in)
            % TAKEPOS take the positional part of the vector/matrix
            %   v = TAKEPOS(twist) is the 3x1 linear velocity vector in
            %   twist (=[w^T v^T]^T) vector.
            out = in(4:6,:);
        end

        function out = takeOri(in)
            % TAKEORI take the rotational part of the vector/matrix
            %   w = TAKEORI(twist) is the 3x1 angular velocity vector in
            %   twist (=[w^T v^T]^T) vector.
            out = in(1:3,:);
        end

        function out = crossmany(a, b)
            % CROSSMANY cross product with multiple vectors
            %   OUT = CROSSMANY(A,B) where a is a single 3x1 vector, and b
            %   is a 3xn matrix (n 3x1 column vectors). out is a 3xn
            %   matrix where each column corresponds to the cross product
            %   of a with each column of b.

            n = size(b,2);
            out = zeros(3,n);
            for i=1:n
                out(:,i) = cross(a,b(:,i));
            end
        end

        function [out, varargout] = isrot(r)
            % ISROT check if it is a rotation matrix
            %   [out, varargout] = ISROT(r)
            %
            %   ISROT  returns true if the matrix is a rotation matrix
            %
            %	T = ISROT(R)
            %	[T REASON] = ISTROT(R)
            %   Copyright (C) 2005, by Brad Kratochvil

            global DebugLevel;

            out = true;
            varargout = {'rot'};

            if ~isequal([3 3], size(r))
                out = false;
                varargout = {'matrix must be 3x3'};
                return;
            end

            if (isempty(DebugLevel)) || (DebugLevel > 1)

                if ~isequalf(eye(3), r' * r, 1e-12)
                    out = false;
                    varargout = {'transpose(r) * r = I violated'};
                    return;
                end

                if ~isequalf(1, det(r))
                    out = false;
                    varargout = {'det(r) == 1 violated'};
                    return;
                end

            end

        end

        function [roll pitch yaw] = rpy(R)
            % RPY convert a rotation matrix R to roll-pitch-yaw angles
            %	[ROLL PITCH YAW] = RPY(R)
            %
            %   Copyright (C) 2005, by Brad Kratochvil

            if ~isrot(R)
                error('R is not a rotation matrix');
            end

            beta = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
            if isequalf(beta, pi/2)
                alpha = 0;
                gamma = atan2(R(1,2), R(2,2));
            elseif isequalf(beta, -pi/2),
                alpha = 0;
                gamma = -atan2(R(1,2), R(2,2));
            else
                alpha = atan2(R(2,1)/cos(beta), R(1,1)/cos(beta));
                gamma = atan2(R(3,2)/cos(beta), R(3,3)/cos(beta));
            end

            roll = gamma;
            pitch = beta;
            yaw = alpha;
        end

        function M = crossSkew(x)
            % CROSSSKEW get a skew matrix from a 3x1 vector
            %   M = CROSSSKEW(X) returns a 3x3 skew matrix constructed from
            %   a 3x1 vector x.
            M = [0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0];
        end

        function R = t2r(T)
            % T2R convert a 4x4 transformation matrix to a rotation matrix
            %   R = T2R(T)

            R = T(1:3,1:3);
        end

        function T = r2t(R)
            % R2T convert a rotation matrix to a 4x4 transformation matrix
            %   T = r2t(R)
            T = [R(1:3,1:3) [0;0;0]; 0 0 0 1];
        end

        function r = rotx(t)
            % ROTX get a rotation matrix for a rotation of t (rad) about x-axis
            %   R = ROTX(t)
            ct = cos(t);
            st = sin(t);
            r =    [1	0	0
                0	ct	-st
                0	st	ct];
        end

        function r = roty(t)
            % ROTY get a 3x3 rotation matrix for a rotation of t (rad) about y-axis
            %   R = ROTY(t)
            ct = cos(t);
            st = sin(t);
            r =    [ct	0	st
                0	1	0
                -st	0	ct];
        end

        function r = rotz(t)
            % rotz get a 3x3 rotation matrix for a rotation of t (rad) about z-axis
            %   R = ROTZ(t)
            ct = cos(t);
            st = sin(t);
            r =    [ct	-st	0
                st	ct	0
                0	0	1];
        end

        function r = transl(x, y, z)
            % TRANSL Create or unpack an SE3 translational transform
            %
            %   This function can create a translational transform matrix
            %   for a translation of (x, y, z), or extract a 3x1 position
            %   vector from a 4x4 transform matrix.
            %
            %   T = TRANSL(x,y,z) or TRANSL(p), where p is a 3x1 position
            %   vector and T is a translational 4x4 transform matrix
            %   (rotational part is an identity matrix).
            %
            %   p = TRANSL(T), where T is a 4x4 transform matrix and p is a
            %   3x1 position vector.
            %
            %   Copyright (C) 1993-2008, by Peter I. Corke

            if nargin == 1,
                if SpatialRobotModel.ishomog(x),
                    r = x(1:3,4);
                elseif ndims(x) == 3,
                    r = squeeze(x(1:3,4,:))';
                else
                    t = x(:);
                    r =    [eye(3)			t;
                        0	0	0	1];
                end
            elseif nargin == 3,
                t = [x; y; z];
                r =    [eye(3)			t;
                    0	0	0	1];
            end
        end

        function h = ishomog(tr)
            % ISHOMOG check if the matrix is a homogeneous tranform matrix
            %   h = ISHOMOG(tr) is 1 when tr is a 4x4 matrix.

            if ndims(tr) == 2,
                h =  all(size(tr) == [4 4]);
            else
                h = 0;
            end
        end

        function u = unit(v)
            % UNIT get a unit vector of v
            %   u = UNIT(v) is in the same direction as v with unit size
            u = v / norm(v,'fro');
        end

        function X = xlt( r )
            % XLT  spatial coordinate transform (translation of origin).
            %   XLT(r)  calculates the coordinate transform matrix from A
            %   to B coordinates for spatial motion vectors, in which frame
            %   B is translated by an amount r (3D vector) relative to
            %   frame A.  r can be a row or column vector.

            X = [  1     0     0    0  0  0 ;
                0     1     0    0  0  0 ;
                0     0     1    0  0  0 ;
                0     r(3) -r(2) 1  0  0 ;
                -r(3)  0     r(1) 0  1  0 ;
                r(2) -r(1)  0    0  0  1];
        end

        function R = angvec2r(theta, k)
            % ANGVEC2R Convert angle and vector orientation to a 3x3 rotation matrix
            %   R = ANGVEC2R(theta, v) Return an orthonormal rotation
            %   matrix, R, equivalent to a rotation of theta about the
            %   vector v.
            %
            %   Copyright (C) 1993-2008, by Peter I. Corke
            cth = cos(theta);
            sth = sin(theta);
            vth = (1 - cth);
            kx = k(1); ky = k(2); kz = k(3);

            R = [kx*kx*vth+cth      ky*kx*vth-kz*sth   kz*kx*vth+ky*sth
                kx*ky*vth+kz*sth   ky*ky*vth+cth          kz*ky*vth-kx*sth
                kx*kz*vth-ky*sth   ky*kz*vth+kx*sth   kz*kz*vth+cth];
        end

        function X = r2X(R)
            % R2X get a 6x6 rotation matrix X from a 3x3 rotation matrix R
            %   X = r2X(R)
            X = [R eye(3);eye(3) R];
        end

        function M = crossSpv(spV)
            % CROSSSPV get a 6x6 skew matrix from 6x1 spatial motion
            %   M = CROSSSPV(spV)
            w = spV(1:3);
            v = spV(4:6);

            M = [SpatialRobotModel.crossSkew(w), zeros(3); ...
                SpatialRobotModel.crossSkew(v), SpatialRobotModel.crossSkew(w)];
        end

        function M = crossSpf(spF)
            % CROSSSPF get a 6x6 skew matrix from 6x1 spatial force
            %   M = CROSSSPF(spF)
            M = -SpatialRobotModel.crossSpv(spF)';
        end

        function  out = skew( in )
            % SKEW  convert 3D vector <--> 3x3 skew-symmetric matrix
            %   S=SKEW(v) and v=SKEW(A) calculate the 3x3 skew-symmetric
            %   matrix S corresponding to the given 3D vector v, and the 3D
            %   vector corresponding to the skew-symmetric component of the
            %   given arbitrary 3x3 matrix A.  For vectors a and b,
            %   SKEW(a)*b is the cross product of a and b.  If the argument
            %   is a 3x3 matrix then it is assumed to be A, otherwise it is
            %   assumed to be v.  SKEW(A) produces a column-vector result,
            %   but SKEW(v) will accept a row or column vector argument.

            if all(size(in)==[3 3])			% do v = skew(A)
                out = 0.5 * [ in(3,2) - in(2,3);
                    in(1,3) - in(3,1);
                    in(2,1) - in(1,2) ];
            else					% do S = skew(v)
                out = [  0,    -in(3),  in(2);
                    in(3),  0,    -in(1);
                    -in(2),  in(1),  0 ];
            end
        end

        function Tout = SRT(Tori, Tref)
            % SRT similarity transform of transform matrices
            %   Tout = SRT(Tori, Tref)
            %
            %   Tori: origininal transformation matrix.
            %   Tref: reference transformation matrix that we want to express
            %         with respect to (a rotational 4x4 transform matrix).
            %   Tout: Tori expressed w.r.t Tref

            if(size(Tori) == [3 3])
                Tout = Tref'*Tori*Tref;
            elseif(size(Tori) == [4 4])
                Tout = SpatialRobotModel.r2t(Tref)'*Tori*SpatialRobotModel.r2t(Tref);
            else
                error('unrecognized transform')
            end
        end

        function Iout = inertia(I, list_index)
            % INERTIA generate an inertia matrix from 1x6 row vector definition.
            %   Iout = INERTIA(I, list_index)
            %
            %   list_index defines the ordering so that 1x6 row vector
            %   definition maps to [ixx iyy izz ixy iyz ixz] for
            %   SpatialRobotModel computations.

            if ~exist('list_index','var')
                list_index = 1:6; %in order by default
            end

            Ixx = I(list_index(1));
            Iyy = I(list_index(2));
            Izz = I(list_index(3));
            Ixy = I(list_index(4));
            Iyz = I(list_index(5));
            Ixz = I(list_index(6));

            Iout = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz];
        end

        function Iout = invInertia(I)
            % INVINERTIA generate 1x6 row vector definition from an inertia matrix.
            %   Iout = INVINERTIA(I)
            %
            %   This function returns a 1x6 row inertia definition from a
            %   3x3 inertia tensor matrix in the order of [Ixx, Iyy, Izz,
            %   Ixy, Iyz, Ixz].

            Ixx = I(1,1);
            Iyy = I(2,2);
            Izz = I(3,3);
            Ixy = I(1,2);
            Iyz = I(2,3);
            Ixz = I(3,1);

            Iout = [Ixx, Iyy, Izz, Ixy, Iyz, Ixz];
        end

        function R = rpy2R(ang_r,ang_p,ang_y)
            % rpy2R convert RPY angles to a 3x3 rotation matrix R
            %   R = rpy2R(ang_r,ang_p,ang_y)
            %
            %   This function returns a rotation matrix R which is a
            %   rotation matrix as a result of rotation of ang_r about x
            %   axis followed by rotation of ang_p about y axis followed by
            %   rotation of ang_y about z axis and all the rotations are
            %   made with respect to base frame.
            R = SpatialRobotModel.rotz(ang_y)*SpatialRobotModel.roty(ang_p)*SpatialRobotModel.rotx(ang_r);
        end

        function [theta, v] = tr2angvec(R)
            % TR2ANGVEC obtain ang and vector from a 4x4 TR or 3x3 R matrix
            %   [theta, v] = TR2ANGVEC(R)

            if ~SpatialRobotModel.isrot(R)
                R = SpatialRobotModel.t2r(R);
            end

            if(norm(R-eye(3))<1e-10) % if R is an identity matrix

                theta = 0;
                v = [0 0 1];

            else

                % Compute theta
                theta = real(acos((1/2)*(R(1,1,:)+R(2,2,:)+R(3,3,:)-1)));

                % Determine initial axis vectors from theta
                v = [ R(3,2,:)-R(2,3,:),...
                    R(1,3,:)-R(3,1,:),...
                    R(2,1,:)-R(1,2,:)] ./ (repmat(2*sin(theta),[1,3]));

                % Handle the degenerate cases where theta is divisible by pi
                singularLogical = mod(theta, cast(pi,'like',R)) == 0;
                if any(singularLogical)
                    vspecial = zeros(3,sum(singularLogical),'like',R);

                    inds = find(singularLogical);
                    for i = 1:sum(singularLogical)
                        [~,~,V] = svd(eye(3)-R(:,:,inds(i)));
                        vspecial(:,i) = V(:,end);
                    end
                    v(1,:,singularLogical) = vspecial;
                end

                % Extract final values
                theta = reshape(theta,[numel(theta) 1]);
                v = reshape(v,[3, numel(v)/3]).';

                if nargout == 0
                    fprintf('Rotation: %f rad x [%f %f %f]\n', theta, v(1), v(2), v(3));
                end
            end
        end

        function I = parallelAxis(I0,mass,R,p)
            % PARALLELAXIS perform parallel axis theorem for an inertia matrix
            %   I = PARALLELAXIS(I0,mass,R,p)
            %
            %   This function translates the moment of inertia of a body
            %   about its CM to the origin. For a point mass, use I0 = 0.
            %   p (= [0;0;0]-[x;y;z]) is a column vector from the
            %   center-of-mass to the origin.

            I0_new = R*I0*R';
            I = I0_new + mass.*(dot(p,p)*eye(3,3) - p*p');
        end

        %% plotting related functions

        function drawFrame(T, len, width)
            % DRAWFRAME draw a simple coordinate frame
            %   DRAWFRAME(T, len, width)
            %
            %   This function draws a frame at a coordinate defined by a
            %   transform matrix T with desired lenth (len) and width
            %   properties.

            if nargin == 1,
                len = 1;
            end

            if ~exist('width', 'var'),
                width = 1;
            end

            colr1 = [1 0 0];
            colr2 = [0 1 0];
            colr3 = [0 0 1];

            origin = T(1:3, 4);             % 1st three elements of 4th column
            X = origin + len*T(1:3, 1);     % point 'len' units out along x axis
            Y = origin + len*T(1:3, 2);     % point 'len' units out along y axis
            Z = origin + len*T(1:3, 3);     % point 'len' units out along z axis

            line([origin(1),X(1)], [origin(2), X(2)], [origin(3), X(3)], 'color', colr1, 'LineWidth', width);
            line([origin(1),Y(1)], [origin(2), Y(2)], [origin(3), Y(3)], 'color', colr2, 'LineWidth', width);
            line([origin(1),Z(1)], [origin(2), Z(2)], [origin(3), Z(3)], 'color', colr3, 'LineWidth', width);

        end

        function handle = drawFrameHandle(T, len, width, color)
            % DRAWFRAMEHANDLE draw a simple coordinate frame and outputs a handle
            %   handle = DRAWFRAMEHANDLE(T, len, width, color)
            %
            %   This function draws a coordinate frame with desired
            %   properties and returns the handle for the frame.

            if nargin == 1,
                len = 1;
            end

            if ~exist('width', 'var'),
                width = 1;
            end

            if ~exist('color', 'var'),
                color.colr1 = [1 0 0];
                color.colr2 = [0 1 0];
                color.colr3 = [0 0 1];
            end

            origin = T(1:3, 4);             % 1st three elements of 4th column
            X = origin + len*T(1:3, 1);     % point 'len' units out along x axis
            Y = origin + len*T(1:3, 2);     % point 'len' units out along y axis
            Z = origin + len*T(1:3, 3);     % point 'len' units out along z axis

            handle.bodies(1) = line([origin(1),X(1)], [origin(2), X(2)], [origin(3), X(3)], 'color', color.colr1, 'LineWidth', width);
            handle.bodies(2) = line([origin(1),Y(1)], [origin(2), Y(2)], [origin(3), Y(3)], 'color', color.colr2, 'LineWidth', width);
            handle.bodies(3) = line([origin(1),Z(1)], [origin(2), Z(2)], [origin(3), Z(3)], 'color', color.colr3, 'LineWidth', width);

        end

        function drawFrameHandleSet(handle, T, len)
            % DRAWFRAMEHANDLESET update simple frame drawing through a handle
            %   DRAWFRAMEHANDLESET(handle, T, len)
            %
            %   This function updates a frame via handle with desired T and
            %   length (len).

            if ~exist('len', 'var'),
                len = 1;
            end

            % colr1 = [1 0 0];
            % colr2 = [0 1 0];
            % colr3 = [0 0 1];

            origin = T(1:3, 4);             % 1st three elements of 4th column
            X = origin + len*T(1:3, 1);     % point 'len' units out along x axis
            Y = origin + len*T(1:3, 2);     % point 'len' units out along y axis
            Z = origin + len*T(1:3, 3);     % point 'len' units out along z axis

            set(handle.bodies(1),'XData',[origin(1),X(1)],'YData', [origin(2), X(2)],'ZData', [origin(3), X(3)]);
            set(handle.bodies(2),'XData',[origin(1),Y(1)],'YData', [origin(2), Y(2)],'ZData', [origin(3), Y(3)]);
            set(handle.bodies(3),'XData',[origin(1),Z(1)],'YData', [origin(2), Z(2)],'ZData', [origin(3), Z(3)]);
        end

        function plotRows(A,color,line_width)
            % plotRows plot lines between links
            %   plotRows(A,color,line_width)
            %
            %   Copyright (C) 2014, by Pat Wensing

            if ~exist('color','var')
                color = 'b'; %world coordinate by default
            end
            if ~exist('line_width','var')
                line_width = 1; %world coordinate by default
            end
            plot3(A(1,:), A(2,:), A(3,:),'color',color,'LineWidth', line_width);
        end

        function scatterRows(A)
            % scatterRows plotting mass at CoM with circles
            %   scatterRows(A)
            %
            %   Copyright (C) 2014, by Pat Wensing

            scatter3(A(1,:), A(2,:), A(3,:));
        end

        function handles = plotEllipse(A, scale, centre, varargin)
            % PLOTELLIPSE Draw an ellipse on the current plot
            %
            %   PLOT_ELLIPSE(A, LS) draws an ellipse defined by X'AX = 0 on the
            %   current plot, centred at the origin, with Matlab line style LS.
            %
            %   PLOT_ELLIPSE(A, C, LS) as above but centred at C=[X,Y].
            %   current plot.  If C=[X,Y,Z] the ellipse is parallel to the XY plane
            %   but at height Z.
            %
            %   H = PLOT_CIRCLE(C, R, options) as above but return handles. For multiple
            %   circles H is a vector of handles, one per circle.
            %
            %   Options::
            %   'edgecolor'   the color of the circle's edge, Matlab color spec
            %   'fillcolor'   the color of the circle's interior, Matlab color spec
            %   'alpha'       transparency of the filled circle: 0=transparent, 1=solid
            %   'alter',H     alter existing circles with handle H
            %
            %   See also PLOT_CIRCLE.
            %   Copyright (C) 1993-2014, by Peter I. Corke

            if size(A,1) ~= size(A,2)
                error('ellipse is defined by a square matrix');
            end

            if size(A,1) > 3
                error('can only plot ellipsoid for 2 or 3 dimenions');
            end

            if nargin < 2
                scale = 1;
            end

            if nargin < 3
                centre = zeros(1, size(A,1));
            end

            if nargin < 4
                varargin = {};
            end

            opt.fillcolor = [];
            opt.alpha = 1;
            opt.edgecolor = 'k';
            opt.alter = [];

            [opt,arglist] = tb_optparse(opt, varargin);

            if ~isempty(opt.alter) & ~ishandle(opt.alter)
                error('RTB:plot_circle:badarg', 'argument to alter must be a valid graphic object handle');
            end

            holdon = ishold();
            hold on

            if size(A,1) == 3
                %= plot an ellipsoid

                % define mesh points on the surface of a unit sphere
                [Xs,Ys,Zs] = sphere();
                ps = scale*[Xs(:) Ys(:) Zs(:)]';

                % warp it into the ellipsoid
                pe = sqrtm(A) * ps;

                % offset it to optional non-zero centre point
                if nargin > 1
                    pe = bsxfun(@plus, centre(:), pe);
                end

                % put back to mesh format
                Xe = reshape(pe(1,:), size(Xs));
                Ye = reshape(pe(2,:), size(Ys));
                Ze = reshape(pe(3,:), size(Zs));

                % plot it
                if isempty(opt.alter)
                    h = mesh(Xe, Ye, Ze, arglist{:});
                else
                    set(opt.alter, 'xdata', Xe, 'ydata', Ye, 'zdata', Ze, arglist{:});

                end

            else
                %= plot an ellipse


                [V,D] = eig(A);

                % define points on a unit circle
                th = linspace(0, 2*pi, 50);
                pc = [cos(th);sin(th)];

                % warp it into the ellipse
                pe = sqrtm(A)*pc;

                % offset it to optional non-zero centre point
                centre = centre(:);
                if nargin > 1
                    pe = bsxfun(@plus, centre(1:2), pe);
                end
                x = pe(1,:); y = pe(2,:);


                if length(centre) > 2
                    % plot 3D data
                    z = ones(size(x))*centre(3);
                    if isempty(opt.alter)
                        h = plot3(x, y, z, varargin{:});
                    else
                        set(opt.alter, 'xdata', x, 'ydata', y, 'zdata', z, arglist{:});
                    end
                else
                    % plot 2D data
                    if isempty(opt.fillcolor)
                        if isempty(opt.alter)
                            h = plot(x, y, arglist{:});
                        else
                            set(opt.alter, 'xdata', x, 'ydata', y, arglist{:});
                        end
                    else
                        if isempty(opt.alter)
                            h = patch(x, y, 0*y, 'FaceColor', opt.fillcolor, ...
                                'FaceAlpha', opt.alpha, 'EdgeColor', opt.edgecolor, arglist{:});
                        else
                            set(opt.alter, 'xdata', x, 'ydata', y, arglist{:});
                        end

                    end
                end
            end
            holdon = ishold;
            hold on

            if nargout > 0
                handles = h;
            end
        end

        %% quaternion-related functions

        function qt = r2q(t)
            % R2Q get a unit quaternion from a homogeneous/rotation matrix
            %   qt = R2Q(t)

            % for sufficiently large q0, this function formulates 2*q0 times the
            % correct return value; otherwise, it formulates 4*|q1| or 4*|q2| or 4*|q3|
            % times the correct value.  The final normalization step yields the correct
            % return value.
            E = t(1:3,1:3);
            tr = trace(E);				% trace is 4*q0^2-1
            v_ = -SpatialRobotModel.skew(E);				% v is 2*q0 * [q1;q2;q3]

            if tr > 0
                qt = [ (tr+1)/2; v_ ];
            else
                E = E - (tr-1)/2 * eye(3);
                E = E + E';
                if E(1,1) >= E(2,2) && E(1,1) >= E(3,3)
                    qt = [ 2*v_(1); E(:,1) ];
                elseif E(2,2) >= E(3,3)
                    qt = [ 2*v_(2); E(:,2) ];
                else
                    qt = [ 2*v_(3); E(:,3) ];
                end
                if qt(1) < 0
                    qt = -qt;
                end
            end

            qt = qt / norm(qt);

            % reverse the sign
            qt(2:4) = -qt(2:4);

        end

        function R = q2R(q)
            % Q2R converts a quaternion into a 3 x 3 Rotation Matrix according to the Euler-Rodrigues formula.
            %   R = Q2R(q) where q = [q0;qv] is a 4x1 quaternion vector.
            %   R = I + 2*q0*hat(qv) + 2*hat(qv)^2

            R = eye(3)+2*q(1)*SpatialRobotModel.crossSkew(q(2:4))+2*SpatialRobotModel.crossSkew(q(2:4))*SpatialRobotModel.crossSkew(q(2:4));
        end

        function Q = qt2Q(qt)
            % QT2Q convert a quaternion into the matrix that converts w into quaternion rate
            %   Q = QT2Q(qt) returns a 4x4 matrix

            q0 = qt(1);
            q1 = qt(2:4);
            Q = [-q1'; q0*eye(3) + SpatialRobotModel.skew(q1)];
        end

        function qt = angvec2qt(theta, u)
            % ANGVEC2QT construct a quaternion from an angle and axis
            %   qt = ANGVEC2QT(theta, u) returns a 4x1 quaternion vector

            qt = zeros(4,1);
            qt(1) = cos(theta/2);
            qt(2:end) = u*sin(theta/2);
        end

        function qout = quatLog(qt)
            % QUATLOG compute the log of a quaternion
            %   qout = QUATLOG(qt) returns a 4x1 quaternion

            for index = length(qt):-1:1
                q(index,:) = norm(qt(index),2);
            end

            % Calculate half the rotation angle
            len = size(q,1);
            normv = arrayfun(@(k) norm(q(k,2:4)),1:len,'UniformOutput',true)';
            th = atan2(normv,q(:,1));

            % Initialize outputs
            qout = zeros(size(q));

            % Calculate logarithm
            tmp = arrayfun(@(k) th(k)*q(k,2:4)/normv(k),1:len,'UniformOutput',false);
            tmp = reshape(cell2mat(tmp'),length(modq),3);
            qout(normv~=0,2:4)=tmp(normv~=0,:);

            qout = qout';
        end

        function qExp = quatExp(qt)
            % QUATEXP quaternion exponential map
            %   qExp = QUATEXP(qt)

            q0 = qt(1);
            q1 = qt(2:4);
            nq1 = norm(q1);
            cosq1 = cos(nq1);
            if(norm(q1)<1e-2)
                sinq1 = 1 - nq1^2/6 + nq1^4/120 - nq1^6/5040;
            else
                sinq1 = sin(nq1)/nq1;
            end
            qExp = exp(q0)*[cosq1; q1*sinq1];
        end

        function qv = quatV(qt)
            % QUATV get the vector part from a quaternion
            %   qv = QUATV(qt)
            qv = qt(2:end);
        end

        function Q = quatMat(qt)
            % QUATMAT compute the conversion matrix (q->Q) from a quaternion
            %   Q = QUATMAT(qt)
            %
            %   This function allows matrix multiplication to perform
            %   quaternion multiplication operation in the same order,
            %   i.e., q*p = Qp.
            %
            %   See also quatMat2.

            q0 = qt(1);
            q1 = qt(2:4);
            Q = [q0 -q1';q1 q0*eye(3) + SpatialRobotModel.skew(q1)];
        end

        function P = quatMat2(qt)
            % QUATMAT2 compute the conversion matrix with negative sign (q->Q)
            %   P = QUATMAT2(qt)
            %
            %   This function allows matrix multiplication to perform
            %   quaternion multiplication operation in the opposite order,
            %   i.e., q*p = P(-)q
            %
            %   See also quatMat.

            q0 = qt(1);
            q1 = qt(2:4);
            P = [q0 -q1';q1 q0*eye(3) - SpatialRobotModel.skew(q1)];
        end

        function [theta, u] = qt2angvec(qt)
            % QT2ANGVEC find the angle and axis of rotation from a quaternion
            %   [theta, u] = QT2ANGVEC(qt)

            q0 = qt(1);
            q1 = qt(2:4);

            th = acos(q0);
            theta = th*2;
            u = q1/norm(q1);
        end

        function qtOut = quatPow(qt, p)
            % QUATPOW quaternion power operation
            %   qtOut = QUATPOW(qt, p) returns quaternion to the power of p

            q0 = qt(1);
            q1 = qt(2:4);

            th = acos(q0);
            u = q1/norm(q1);
            qtOut = norm(qt)^p*[cos(p*th); u*sin(p*th)];
        end

        function qt_t = quatInterp(qt_0,qt_1,t)
            % QUATINTERP apply SLERP for two consequtive quaternions with t
            %   qt_t = QUATINTERP(qt_0, qt_1, t) where t is a scalar in [0,1].
            %   qt_t = qt_0*(qt_0'*qt_1)^t

            qt_delta = SpatialRobotModel.quatPow((SpatialRobotModel.quatMat(qt_0)'*qt_1),t);
            qt_t = SpatialRobotModel.quatMat(qt_0)*qt_delta;
        end

        function qConj = quatConj(q)
            % QUATCONJ quaternion conjugate
            %   qConj = QUATCONJ(q)

            qConj = [q(1);-q(2:4)];
        end

        function qtUnit = quatNorm(qt)
            % QuatNorm nomalize a quaternion into a unit quaternion
            %   qtUnit = QUATNORM(qt)
            qtUnit = qt / norm(qt);
        end

        function E = quatE(qt)
            % QUATE compute the E matrix (mu*I - S(eps))
            %   E = QUATE(qt)

            E = qt(1)*eye(3) - SpatialRobotModel.skew(qt(2:4));
        end

        function qt = quatRand(bounds)
            % QUATRAND random quaternion rotation
            %   qt = QUATRAND(bounds)
            %
            %   bounds: [thetaMin, thetaMax] -- lower and upper bounds

            if ~exist('bounds','var')
                bounds = [-pi pi];
            end

            % create a random rotation in terms of angle and axis
            theta = bounds(1) + (bounds(2)-bounds(1))*rand();
            u_axis = SpatialRobotModel.unit(rand(3,1));
            qt = SpatialRobotModel.angaxis2qt(theta, u_axis);
        end

        function [theta, u_axis] = angleAxisRand(bounds)
            % ANGLEAXISRAND get random angle and axis
            %   [theta, u_axis] = ANGLEAXISRAND(bounds)
            %
            %   bounds: [thetaMin, thetaMax] -- lower and upper bounds

            if ~exist('bounds','var')
                bounds = [-pi pi];
            end

            % create a random rotation in terms of angle and axis
            theta = bounds(1) + (bounds(2)-bounds(1))*rand();
            u_axis = SpatialRobotModel.unit(rand(1,3));
        end

        function qtOut = quatSRT(qt_ori, qt_ref)
            % QUATSRT apply similarity transform to quaternions
            %   qtOut = QUATSRT(qt_ori, qt_ref)
            %
            %   qt_out expresses qt_ori with respect to the qt_ref frame.

            % Rout = Rref'*Rori*Rref
            qtOut = SpatialRobotModel.r2q(SpatialRobotModel.SRT(SpatialRobotModel.q2R(qt_ori), SpatialRobotModel.q2R(qt_ref)));
        end

        function [qtTwist, qtSwing] = quatST(qt, v)
            % QUATST Swing-twist decomposition
            %   [qtTwist, qtSwing] = QUATST(qt, v)
            %
            %   This function decompose a quaternion into twist and swing
            %   quaternions (swing after twist). The original quaternion
            %   can be reconstructed as q = q_swing*q_twist (so q_swing =
            %   q*q_twist').
            %
            %   INPUT
            %       qt: 4x1 quaternion vector
            %       v: 3x1 desired axis of rotation (twist axis)
            %   OUTPUT
            %       qtTwist: a twist quaternion vector
            %       qtSwing: a swing quaternion vector

            % compute twist quaternion
            qtTwist = SpatialRobotModel.unit([qt(1); SpatialRobotModel.proj(qt(2:4), v)]);
            % compute swing quaternion
            qtSwing = SpatialRobotModel.quatMat2(qt_twist)'*qt;
        end

        %% kinematics related functions

        function J_inv = SRINV(J, method, tolerance, lambda_max)
            % SRINV singularity-robust inverse implementation
            %   J_inv = SRINV(J, method, tolerance, lambda_max)
            %
            %   This function computes a singularity-robust inverses for a
            %   Jacobian matrix.
            %
            %   J: a Jacobian matrix
            %   method: an integer (0-4) that selects the desired
            %           inverse implementation
            %           0: matlab's pinv
            %           1: damped least-square inverse
            %           2: remove near-singular values
            %           3: user-defined accuracy
            %           4: damped least-square inverse with numerical filtering
            %           5: damped least-square inverse with numerical filtering (updated version)
            %   tolerance: a tolerance
            %   lamda_max: a maximum lambda value used in method 3

            if ~exist('method','var')
                method = 1;
            end

            if ~exist('tolerance','var')
                tolerance = 1e-4;
            end

            if ~exist('lambda_max','var')
                lambda_max = 1e-1;
            end

            if(method == 0)
                % method 0: matlab's pinv
                J_inv = pinv(J);

            elseif(method == 1)
                % method 1: damped least-square inverse
                lambda = tolerance;

                J_inv = J'/(J*J' + lambda^2*eye(size(J,1)));

            elseif(method == 2)
                % method 2: remove near-singular values

                [U,S,V] = svd(J);
                S(abs(S)<(sum(sum(S))*tolerance)) = 0; % removes near-singular values.
                Stemp = 1./S;
                % This take the inverse of non-zero terms... I'm sure there
                % is faster way
                Stemp(isinf(Stemp)) = 0;
                J_inv = V * Stemp' * U';

            elseif(method == 3)
                % method 3: user-defined accuracy

                thres = tolerance;
                %     lambda_max = 1e-1;
                %     [u,s,v] = svd(J);
                s = svd(J);
                s_min = min(s);
                if(s_min >= thres)
                    lambda = 0;
                else
                    lambda = (1-(s_min/thres)^2)*lambda_max;
                end

                J_inv = J'/(J*J' + lambda*eye(size(J,1)));

            elseif(method == 4)
                % method 4: damped least-square inverse with numerical filtering
                [u,s,v] = svd(J);
                u_min = u(:,end);

                lambda = tolerance;
                J_inv = J'/(J*J' + lambda*u_min*u_min');

            elseif(method == 5)
                % method 5: damped least-square inverse with numerical
                % filtering (multiple singular output vectors)
                tolerance = 1e-1;

                [u,s,v] = svd(J);
                s_min = min(diag(s));

                % select the singular output vectors that correspond to
                % near zero singular values
                u_sum = zeros(size(J,1),1);
                if(s_min < tolerance)
                    indices_u = find(diag(s) < tolerance);
                    for i=1:numel(indices_u)
                        u_sum = u_sum + u(:,indices_u(i))*u(:,indices_u(i))';
                    end
                end

%                 beta = 1e-1;
%                 lambda = 1e-4;
%
%                 J_inv = J'/(J*J' + lambda*eye(size(J,1)) + beta*u_sum);

                beta_max = 2e-1; %1e-1
                lambda = 1e-6;
                thres = tolerance;
                if(s_min > thres)
                    beta = 0;
                else
                    beta = (1-(s_min/thres)^2)*beta_max;
                end
                J_inv = J'/(J*J' + lambda^2*eye(size(J,1)) + beta*u_sum);

            elseif(method == 6)
                % method 6: Jacobian Transpose
                J_inv = J';
            end

        end

        function Nproj = nullProj(J, method, tolerance)
            % NullProj compute nullspace projection matrix
            %   J_inv = SRINV(J, method, tolerance, lambda_max)
            %
            %   This function computes nullspace projection matrix for a
            %   Jacobian matrix.
            %
            %   J: a Jacobian matrix
            %   method: an integer (0-1) that selects the desired
            %           inverse implementation
            %           0: using matlab's pinv
            %           1: damped least-square inverse
            %           2: remove near-singular values
            %           3: user-defined accuracy
            %           4: damped least-square inverse with numerical filtering
            %           5: damped least-square inverse with numerical filtering (updated version)
            %   tolerance: a tolerance
            %   lamda_max: a maximum lambda value used in method 3

            if ~exist('method','var')
                method = 0;
            end

            if ~exist('tolerance','var')
                tolerance = 1e-5;
            end

            if (method == 0)
                %- regular nullspace projection
                Nproj = (eye(size(J,2))-SpatialRobotModel.SRINV(J,0)*J);

            elseif (method == 1)
                %- compute continuous nullspace projection

                % get SVD of Jacobian
                [U,S,V] = svd(J);

                m = size(S,1);
                n = size(S,2);
                A = zeros(n, n);

                % for the first m-1 directions, set it to 1
                a_i = ones(1, m);

                % compute the activation matrix components based on
                % singular values and the threshold
                s_thres = 0.05;
                for i=1:m
                    if(S(i,i)^2 < s_thres)
                        a_i(i) = S(i,i)^2/s_thres;
                    else
                        a_i(i) = 1;
                    end
                end

                % m x m block square matrix
                A(1:m,1:m) = diag(a_i);

                % compute nullspace projection
                Nproj = (eye(size(J,2)) - V*A*V')';
            end

        end

        function TR = decompMul(T1, T2, method)
            % DECOMPMUL decompositional multiplication for two transform matrices
            %   TR = DECOMPMUL(T1, T2, method)
            %
            %   This implements decompositional multiplication for two 4x4 HT matrices
            %   it will perform both T1 and T2 tranformations in sequential order
            %   w.r.t base coordinate system and returns the output pose.
            %
            %   usage: TR = DECOMPMUL(T1, T2, method)
            %           - method = 1 (the description below)
            %           - method = 2 (use similarity tranformation)
            %
            %   TR = | R2*R1 | p2+p1 |
            %        | 0 0 0 |   1   |

            if ~exist('method','var')
                method = 1;
            end

            if (method == 1)
                TR = [SpatialRobotModel.t2r(T2)*SpatialRobotModel.t2r(T1), SpatialRobotModel.transl(T2)+SpatialRobotModel.transl(T1); 0 0 0 1];
            elseif (method == 2)
                TR = T1*SpatialRobotModel.r2t(SpatialRobotModel.t2r(T1)')*T2*SpatialRobotModel.r2t(SpatialRobotModel.t2r(T1));
            end

        end

        function out = computePerformanceIndex(J, select_method, u)
            % COMPUTEPERFORMANCEINDEX compute performance indices of jacobian
            %   out = COMPUTEPERFORMANCEINDEX(J, select_method, u)
            %
            %   Based on select_method value, a different performance index
            %   is calculated:
            %       1) velocity manipulability
            %       2) force manipulability
            %       3) condition number
            %       4) velocity transmission ratio along u
            %       5) force transmission ratio along u

            if ~exist('u','var')
                u = [1; 0; 0]; % unit vector
            end

            if(select_method == 1)
                %-- compute velocity manipulability
                out = sqrt(det(J * J'));

            elseif(select_method == 2)
                %-- compute force manipulability
                out = 1/sqrt(det(J * J'));

            elseif(select_method == 3)
                %-- compute condition number
                out = cond(J * J');

            elseif(select_method == 4)
                %-- compute the velocity transmission ratio along u
                out = 1/sqrt(u'/(J * J')*u);

            elseif(select_method == 5)
                %-- compute the force transmission ratio along u
                out = sqrt(u'/(J * J')*u);
            end
        end

        function numgrad = ngrad(Pind, theta)
            % NGRAD compute a numerical gradient
            %   numgrad = NGRAD(Pind, theta)
            %
            %   This function computes the gradient using "finite
            %   differences" and gives us a numerical estimate of the
            %   gradient. numgrad = NGRAD(J, theta) computes the numerical
            %   gradient of the function J around theta. Calling y =
            %   J(theta) should return the function value at theta.
            %
            %   Notes: The following code implements numerical gradient
            %   checking, and returns the numerical gradient.It sets
            %   numgrad(i) to (a numerical approximation of) the partial
            %   derivative of J with respect to the i-th input argument,
            %   evaluated at theta. (i.e., numgrad(i) should be the
            %   (approximately) the partial derivative of J with respect to
            %   theta(i).)

            numgrad = zeros(size(theta));
            perturb = zeros(size(theta));
            e = 1e-4;

            for p = 1:numel(theta)
                % Set perturbation vector
                perturb(p) = e;
                loss1 = Pind(theta - perturb);
                loss2 = Pind(theta + perturb);
                % Compute Numerical Gradient
                numgrad(p) = (loss2 - loss1) / (2*e);
                perturb(p) = 0;
            end
        end

        %% other useful functions

        function r = randRange(N,lb,ub)
            % RANDRANGE create a random number (N x 1) in a range
            %   r = RANDRANGE(N,lb,ub)

            r = (ub-lb).*rand(N,1) + lb;
        end

        function out = struct2num(s)
            % STRUCT2NUM turn struct of string items to numbers (recursively)
            %   out = STRUCT2NUM(s)

            F = fieldnames(s);
            nsub = length(F);
            out = s;
            for i=1:nsub
                % if the subitem is struct, do a recursive operation
                if(isstruct(s.(F{i})))
                    out.(F{i}) = SpatialRobotModel.struct2num(out.(F{i}));
                else
                    out.(F{i}) = str2num(out.(F{i}));
                    % if string put it back to string
                    if isempty( out.(F{i}) )
                        out.(F{i}) = s.(F{i});
                    end
                end
            end
        end

        function out = struct2str(s)
            % STRUCT2STR turn struct of number items to strings (recursively)
            %   out = STRUCT2STR(s)

            F = fieldnames(s);
            nsub = length(F);
            out = s;
            for i=1:nsub
                % if the subitem is struct, do a recursive operation
                if(isstruct(s.(F{i})))
                    out.(F{i}) = SpatialRobotModel.struct2num(out.(F{i}));
                else
                    out.(F{i}) = num2str(out.(F{i}));
                    % if string put it back to string
                    if isempty( out.(F{i}) )
                        out.(F{i}) = s.(F{i});
                    end
                end
            end
        end

        function c = proj(a,b)
            % PROJ  Projection of a vector A onto vector B.
            %   C = PROJ(A,B) returns the vector of the projection of A
            %   onto B. A and B must be vectors of the same length.

            if nargin < 2
                error('proj:Input', 'Requires exactly two input arguments.');
            end

            % Check dimensions
            if any(size(a) ~= size(b))
                error('MATLAB:proj:InputSizeMismatch', 'A and B must be same size.');
            end

            c = dot(a,b)/norm(b)^2*b;
        end

    end

    %%% functions for class object %%%%

    methods

        %% model related functions

        function model_param = model2struct(model)
            % MODEL2STRUCT copy the model parameters into a struct
            %   model_param = MODEL2STRUCT(model)
            %
            %   This function creates a struct that holds all the model
            %   parameters used by mex functions.

            model_param.NB = model.NB;
            model_param.NJ = model.NJ;
            model_param.N_joint = model.N_joint;
            model_param.N_ext = model.N_ext;
            model_param.ag_0 = model.ag_0;

            for j = 1:model.N_ext
                model_param.parent(j) = model.parent{j};
            end

            model_param.pi_R_i = zeros(3, 3, model.N_ext);
            for j = 1:model.N_ext
                model_param.pi_R_i(:, :, j) = model.pi_R_i{j};
            end

            model_param.phi = zeros(model.N_ext, 6);
            model_param.phi_base = model.phi{1};
            for j = 2:model.NB
                model_param.phi(j, :) = model.phi{j};
            end

            model_param.pi_p_i = zeros(model.N_ext, 3);
            for j = 1:model.N_ext
                model_param.pi_p_i(j, :) = model.pi_p_i{j}(1:3);
            end

            for j = 1:model.N_ext
                model_param.joints(j,:) = zeros(1,6);
                if ~isempty(model.joints{j})
                    for k=1:numel(model.joints{j})
                        joints_tmp = model.joints{j};
                        model_param.joints(j,k) = joints_tmp(k);
                    end
                end
            end

            model_param.m = zeros(model.N_ext, 1);
            for j = 1:model.N_ext
                model_param.m(j) = model.m{j};
            end

            model_param.I_com = zeros(3, 3, model.N_ext, 1);
            for j = 1:model.N_ext
                model_param.I_com(:, :, j) = model.I_com{j};
            end

            model_param.i_p_com = zeros(model.N_ext, 3);
            for j = 1:model.N_ext
                model_param.i_p_com(j, :) = model.i_p_com{j};
            end

            model_param.T0 = zeros(4, 4, model_param.N_ext);
            model_param.R0 = zeros(3, 3, model_param.N_ext);
            model_param.i_X_pi = zeros(6, 6, model_param.N_ext);
            model_param.I = zeros(6, 6, model_param.N_ext);
        end

        function struct2model(model, model_param)
            % STRUCT2MODEL copy the model parameters from struct to model
            %   STRUCT2MODEL(model, model_param)
            %
            %   This function can be used to update the model object's
            %   parameters once mex functions update the struct variable's
            %   parameters.

            if ~exist('model_param','var')
                model_param = model.model_param;
            end

            for j = 1:model.N_ext
                model.T0{j} = model_param.T0(:,:,j);
                model.R0{j} = model_param.R0(:,:,j);
                model.i_X_pi{j} = model_param.i_X_pi(:,:,j);
                model.I{j} = model_param.I(:,:,j);
            end
        end

        function result = isInSubtree(model, leaf, candidate)
            % ISINSUBTREE check if the candidate is in the tree given a leaf
            %   result = ISINSUBTREE(model, leaf, candidate)
            %
            %   This function checks if the candidate node (index) is in
            %   the tree following the tree structure starting from a leaf
            %   node provided. Returns 1 or 0.

            parent = model.parent{leaf};
            while parent ~= 0
                if parent == candidate
                    result = 1;
                    return
                else
                    parent = model.parent{parent};
                end
            end
            result = 0;
        end

        function index = FindBody(model,str)
            % FINDBODY find the index of the body by name
            %   index = FINDBODY(model,str)
            %
            %   This function returns the index of the body that has the
            %   provided string name. Returns -1 if the body does not exist.

            a = model.bodyname;
            ind = ind2sub(size(a),find(cellfun(@(x)strcmp(x,str),a)));

            % error if the body name is incorrect
            if(isempty(ind))
                error('link name is incorrect or the body does not exist!');
                index = -1;
            else
                index = ind;
            end
        end

        function pos = posBody(model, link, pt)
            % POSBODY compute the position of a point in a body
            %   pos = POSBODY(model, link, pt)
            %
            %   This function returns a position vector of the origin of a
            %   chosen link's link coordinate frame w.r.t base frame. If a
            %   3x1 point (pt) is provided, it will add the point vector to
            %   the origin position vector. That is, a point is defined
            %   w.r.t link coordinate frame and the returned position
            %   vector is the point expressed w.r.t base frame.

            if ~exist('pt','var')
                pt = zeros(3,1); %world coordinate by default
            end
            pt = [pt;1];
            T0 = model.T0{link};
            pos = T0*pt;
            pos = pos(1:3);
        end

        function constructModelFromURDF(model, URDF_filename)
            % CONSTRUCTMODELFROMURDF construct model from URDF file
            %   CONSTRUCTMODELFROMURDF(model, URDF_filename)
            %
            %   This function creates a model object from a URDF file
            %   (specified by URDF_filename). It parses the URDF file, and
            %   fills all the model parameters necessary for visualization
            %   and computations. NOTE that this function was not meant to
            %   work for any URDF files, and it needs some cleanup work on
            %   URDF file to some extent.

            %% load URDF file
            model.URDF_name = URDF_filename;
            URDF_data = xml2struct(URDF_filename);
            robot = URDF_data.robot;

            % define collision shape prestored filename
            model.collisionShapeFile = strcat(model.URDF_name(1:end-9),'_caps_param.mat');

            %% get link and joint numbers
            %-- get number of rigid links
            model.NB = numel(robot.link);
            model.N_ext = model.NB;

            %-- get number of revolute joints and their names
            % remove fixed joints
            model.jointNames = {};
            cntJnt = 0;
            for i=1:numel(robot.joint)
                if(strcmp(robot.joint{i}.Attributes.type,'revolute'))
                    cntJnt = cntJnt + 1;
                    model.jointNames{cntJnt} = robot.joint{i}.Attributes.name;
                else
                    % remove virtual links
                    if (i > 1)
                        % diminish NB only if the joint type is fixed
                        model.NB = model.NB-1;
                        fprintf('%s Link is a virtual link!\n', robot.joint{i}.child.Attributes.link);
                    end
                end
            end
            model.NJ = cntJnt;

            % set the number of joints including floating base 6DOF
            model.N_joint = model.NJ + 6;

            %% set link and joint names
            q = zeros(1,model.N_joint);

            model.linkNames = cell(model.NB,1);
            model.mesh_name = cell(model.NB,1);

            % set link and frame names
            for i=1:model.N_ext
                if (i <= model.NB)
                    model.linkNames{i} = robot.link{i}.Attributes.name;
                    model.mesh_name{i} = model.linkNames{i};
                else
                    % The names of the virtual link becomes end-effectors' frameNames
                    model.frameNames{i-model.NB} = robot.link{i}.Attributes.name;
                end
            end

            %% get parameters from links and joints
            for i = 1:numel(robot.joint)
                robot_joint_i = model.struct2num(robot.joint{i});
                model.URDF_par(i).jointname = robot_joint_i.Attributes.name;
                model.URDF_par(i).type = robot_joint_i.Attributes.type;
                model.URDF_par(i).parentlink = robot_joint_i.parent.Attributes.link;
                model.URDF_par(i).childlink = robot_joint_i.child.Attributes.link;
                % joint coordinate
                model.URDF_par(i).joint_rpy = robot_joint_i.origin.Attributes.rpy;
                model.URDF_par(i).joint_xyz = robot_joint_i.origin.Attributes.xyz;
                % some joints do not have field: axis
                try
                    model.URDF_par(i).joint_axis = robot_joint_i.axis.Attributes.xyz;
                end
                % some joints do not have field: limit
                try
                    model.URDF_par(i).joint_limit = robot_joint_i.limit.Attributes;
                end
            end

            for i = 1:numel(robot.joint)
                for j = 1:numel(robot.link)
                    % find child link for the joint
                    if strcmp(model.URDF_par(i).childlink, robot.link{j}.Attributes.name)
                        robot_link_j = model.struct2num(robot.link{j});

                        % some link do not have field: collision
                        try
                            model.URDF_par(i).collision.origin = robot_link_j.collision.origin.Attributes;
                        end
                        try
                            model.URDF_par(i).collision.geometry.cylinder = robot_link_j.collision.geometry.cylinder.Attributes;
                        end
                        try
                            model.URDF_par(i).collision.geometry.sphere = robot_link_j.collision.geometry.sphere.Attributes;
                        end

                        % get the properties for mesh files
                        try
                            model.URDF_par(i).visual.origin = robot_link_j.visual.origin.Attributes;
                        end

                        % some link do not have field: inertial
                        if(isfield(robot_link_j, 'inertial'))
                            model.URDF_par(i).mass = robot_link_j.inertial.mass.Attributes.value;
                            model.URDF_par(i).mass_xyz = robot_link_j.inertial.origin.Attributes.xyz;
                            model.URDF_par(i).mass_rpy = robot_link_j.inertial.origin.Attributes.rpy;
                            % inertia in urdf frame
                            I = robot_link_j.inertial.inertia.Attributes;
                            % our model uses [ixx iyy izz ixy iyz ixz]
                            model.URDF_par(i).inertia = [I.ixx, I.iyy, I.izz, I.ixy, I.iyz, I.ixz];
                        else
                            model.URDF_par(i).mass = 0;
                            model.URDF_par(i).mass_xyz = zeros(1,3);
                            model.URDF_par(i).mass_rpy = zeros(1,3);
                            model.URDF_par(i).inertia = [1, 1, 1, 0, 0, 0];
                        end
                    end
                end
            end

            %% construct a spatial model from URDF struct

            %== base body (floating) - to be modified to fixed later
            KBody = 1;
            model.bodyname{KBody} = model.linkNames{1};
            model.parent{KBody} = 0; % parent body (0: none)
            model.phi{KBody}    = [zeros(3) eye(3); eye(3) zeros(3)]; % joint moving this body
            index = find(cellfun( @(x) isequal(x,model.bodyname{KBody}), {model.URDF_par.childlink}));
            if isempty(index)
                index = 1;
            end
            model.i_p_com{KBody}  = model.URDF_par(index).mass_xyz';
            model.m{KBody}      = model.URDF_par(index).mass;
            model.I_com{KBody}  = model.inertia(model.URDF_par(index).inertia, model.inertia_index);
            model.pi_p_i{KBody}   = model.URDF_par(index).joint_xyz'; % position of the joint
            model.pi_R_i{KBody}   = eye(3); % rotation of the joint
            model.type{KBody} = model.URDF_par(index).type; % type of the joint: revolute or fixed
            model.joints{KBody}  = [1:KBody + 5]; % joints (moving this body)
            model.q{KBody}      = q(1:KBody + 5)'; % initial joint angles
            model.p_mesh_offset{KBody} = zeros(3,1);
            model.R_mesh_offset{KBody} = zeros(3,1);

            %=== body
            kJoint = 7;
            for i = 2 : model.NB
                KBody = KBody + 1;
                model.bodyname{KBody} = model.linkNames{i}; % body: upper-shoulder, joint: s0
                index = find(cellfun( @(x) isequal(x,model.bodyname{KBody}), {model.URDF_par.childlink}));
                model.parent{KBody} = model.FindBody(model.URDF_par(index).parentlink); % parent body name
                model.phi{KBody}    = [ model.URDF_par(index).joint_axis' ; zeros(3,1)] ; % axis (xyz)
                model.i_p_com{KBody}  = model.URDF_par(index).mass_xyz';
                model.m{KBody}      = model.URDF_par(index).mass;
                model.I_com{KBody}  = model.inertia(model.URDF_par(index).inertia, model.inertia_index);
                model.pi_p_i{KBody}   = model.URDF_par(index).joint_xyz'; % origin xyz
                model.pi_R_i{KBody}   = model.rpy2R(model.URDF_par(index).joint_rpy(1),model.URDF_par(index).joint_rpy(2),model.URDF_par(index).joint_rpy(3)); % origin rpy
                model.type{KBody} = model.URDF_par(index).type; % type of the joint: revolute or fixed
%                 model.revoluteJointNames{numel(model.revoluteJointNames)+1} = model.jointNames{i-1};
                if(strcmp(model.type{KBody},'revolute'))
                    model.joints{KBody} = kJoint;
                    kJoint = kJoint + 1;

                    % add the current joint name to the revolute joint names
                    model.revoluteJointNames{numel(model.revoluteJointNames)+1} = model.jointNames{i-1};
                elseif(strcmp(model.type{KBody},'fixed'))
                    % update model parameters
                    model.NJ = model.NJ - 1;
                    model.N_joint = model.N_joint - 1;
                    model.joints{KBody} = []; % empty
                end
                model.q{KBody}      = q(KBody + 5)';

                % get joint limit parameters
                if(isfield(model.URDF_par(index), 'joint_limit'))
                    if(~isempty(model.URDF_par(index).joint_limit))
                        model.joint_pos_min{KBody} = model.URDF_par(index).joint_limit.lower;
                        model.joint_pos_max{KBody} = model.URDF_par(index).joint_limit.upper;
                        model.joint_vel_max{KBody} = model.URDF_par(index).joint_limit.velocity;
                        model.joint_tq_max{KBody}  = model.URDF_par(index).joint_limit.effort;
                    else
                        model.joint_pos_min{KBody} = -2*pi;
                        model.joint_pos_max{KBody} = 2*pi;
                        model.joint_vel_max{KBody} = 10;
                        model.joint_tq_max{KBody}  = 300;
                    end
                else
                    % if the joint limit field does not exist, set default values
                    model.joint_pos_min{KBody} = -2*pi;
                    model.joint_pos_max{KBody} = 2*pi;
                    model.joint_vel_max{KBody} = 10;
                    model.joint_tq_max{KBody}  = 300;
                end

                if(isfield(model.URDF_par(index),'collision'))
                    % get collision shape
                    model.collision{KBody} = model.URDF_par(index).collision;
                end

                % get visual properties
                if(isfield(model.URDF_par(index),'visual'))
                    if(numel(model.URDF_par(index).visual)==0)
                        model.p_mesh_offset{KBody} = zeros(3,1);
                        model.R_mesh_offset{KBody} = zeros(3,1);
                    else
                        model.p_mesh_offset{KBody} = model.URDF_par(index).visual.origin.xyz';
                        model.R_mesh_offset{KBody} = model.URDF_par(index).visual.origin.rpy';
                    end
                else
                    model.p_mesh_offset{KBody} = zeros(3,1);
                    model.R_mesh_offset{KBody} = zeros(3,1);
                end
            end

            %=== end-effector body (only for frames)
            for i = 1 : numel(model.frameNames)
                KBody = KBody + 1;
                model.bodyname{KBody} = model.frameNames{i};
                index = find(cellfun( @(x) isequal(x,model.bodyname{KBody}), {model.URDF_par.childlink}));
                model.parent{KBody} = model.FindBody(model.URDF_par(index).parentlink);
                model.phi{KBody}    = [ model.URDF_par(index).joint_axis' ; zeros(3,1)] ; % axis (xyz)
                model.i_p_com{KBody}  = model.URDF_par(index).mass_xyz';
                model.m{KBody}      = model.URDF_par(index).mass;
                model.I_com{KBody}  = model.inertia(model.URDF_par(index).inertia);
                model.pi_p_i{KBody}   = model.URDF_par(index).joint_xyz';
                model.pi_R_i{KBody}   = model.rpy2R(model.URDF_par(index).joint_rpy(1),model.URDF_par(index).joint_rpy(2),model.URDF_par(index).joint_rpy(3)); % origin rpy
                model.type{KBody} = model.URDF_par(index).type; % type of the joint: revolute or fixed
                model.joints{KBody} = [];
                model.q{KBody}       = [];
                model.p_mesh_offset{KBody} = zeros(3,1);
            end

            %= set the mesh color
            for i=1:model.NB+1
                % grey
                model.mesh_color{i} = [196 205 213]/255; % color of the body links
            end

            % copy model parameters into a struct variable
            model.model_param = model.model2struct();

            %= backup model parameters
            model.param_backup.i_p_com = model.i_p_com;
            model.param_backup.m = model.m;
            model.param_backup.I_com = model.I_com;
            model.param_backup.pi_p_i = model.pi_p_i;
            model.param_backup.pi_R_i = model.pi_R_i;
            model.param_backup.revoluteJointNames = model.revoluteJointNames;

        end

        %% kinematics/dynamics compuations

        function jac = jacob(model, link, coord)
            % JACOB compute jacobian for the provided link frame
            %   jac = JACOB(model, link, coord)
            %
            %   This function computes Jacobian matrix for the provided
            %   link index. If coord is 0, Jacobian is expressed w.r.t base
            %   frame, and if coord is 1, it is expressed w.r.t the link
            %   coordinate frame.
            %
            %   See also jacob_dot.

            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end

            if(model.use_mex)
                jac = jacobMex_mex(model.model_param, link, coord);
            else
                jac = jacob_pt(model,link,[0 0 0]',coord);
            end
        end

        function jac = jacob_pt(model,link,pt,coord)
            % JACOB_PT compute jacobian of a certain position of the body i
            %   jac = JACOB_PT(model,link,pt,coord)
            %
            %   This function computes a Jacobian for a certain point w.r.t
            %   a link. jacob function uses this function with zero
            %   position vector.

            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jac = zeros(6,model.N_joint);

            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                c_X_j = SpatialRobotModel.xlt(pt);
                parent = link;
            end

            while parent ~= 0
                j = parent;
                if(j==1)
                    c_X_j = c_X_j*model.i_X_pi{j};
                    if(~isempty(model.joints{j}))
                        jac(:,model.joints{j}) = c_X_j*model.phi{j};
                    end
                else
                    if(~isempty(model.joints{j}))
                        jac(:,model.joints{j}) = c_X_j*model.phi{j};
                    end
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end

            if(coord == 0)
                %Transform to global
                jac = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jac;
            else
            end
        end

        function Jd = jacob_dot(model, ind, qd_a)
            % JACOB_DOT computes jacobian dot using 3D vector representation
            %   Jd = JACOB_DOT(model, ind, qd_a)
            %
            %   This function computes the derivative of the Jacobian for a
            %   link coordinate frame based on my analytical derivation
            %   using 3D vector notation. It has been validated, but it
            %   needs to be further optimized for speed later.
            %   TODO: cleanup and optimize the code.

            % compute jacobian
            J = model.active(model.jacob(ind));
            Jv = model.takePos(J);
            Jw = model.takeOri(J);

            %%% initialization
            p_0_k = model.transl(model.T0{ind});
            p_0_i_1_i = zeros(3, model.NJ);
            pdot_0_i_1_i = zeros(3, model.NJ);
            p_0_i_1 = zeros(3, model.NJ);
            p_0_i_1_k = zeros(3, model.NJ);
            pdot_0_i_1_k = zeros(3, model.NJ);

            w = zeros(3, model.NJ);

            Jdw = zeros(3, model.NJ);
            Jdv = zeros(3, model.NJ);

            Njoints = model.joints{model.parent{ind}}-6;

            %%% compute Jdw
            for i = 1:Njoints
                if (i > 1)
                    w_i_1 = w(:,i-1);
                else
                    w_i_1 = zeros(3,1);
                end

                %-- compute Jw dot
                Jdw(:,i) = cross(w_i_1, Jw(:,i));

                %-- compute quantities for Jdv
                % compute w_0_0_i
                w(:,i) = qd_a(i)*Jw(:,i) + w_i_1;

                % compute link position vector w.r.t base frame
                p_0_i_1_i(:,i) = model.transl(model.T0{i+2}) - model.transl(model.T0{i+1});

                % sum of w_i x p_0_i_1_i
                pdot_0_i_1_i(:,i) = cross(w(:,i), p_0_i_1_i(:,i));

                % p_0_i_1
                p_0_i_1(:,i) = model.transl(model.T0{i+1});

                % p_i_1_k
                p_0_i_1_k(:,i) = p_0_k - p_0_i_1(:,i);
            end

            %%% compute Jdv

            for i=1:Njoints
                % sum of w_i and p_0_i_1_i crossproduct
                pdot_0_i_1_k(:,i) = sum(pdot_0_i_1_i(:,i:end),2);

                % compute Jvdot
                Jdv(:,i) = cross(Jdw(:,i), p_0_i_1_k(:,i)) + cross(Jw(:,i), pdot_0_i_1_k(:,i));
            end

            Jd = [Jdw; Jdv];
        end

        function jacdot = jacob_dot_spatial(model,link,qd,coord)
            % JACOB_DOT compute jacobian dot for the provided link frame
            %   jacdot = JACOB_DOT_SPATIAL(model,link,qd,coord)
            %
            %   This function computes jacobian dot in terms of spatial
            %   vector computation, and it is currently incorrect and it
            %   needs to be revisited at some point. PLEASE USE JACOB_DOT
            %   FUNCTION INSTEAD.
            %   TODO: revisit the implementation and fix the problem.

            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jacdot = zeros(6,model.N_joint);

            % compute the velocity of each body
            %% RNEA outward recursion
            for i=1:model.NB
                p = model.parent{i};
                if p == 0
                    v{i} = model.phi{i}*qd(model.joints{i});
                else
                    v{i} = model.i_X_pi{i}*v{p} + model.phi{i}*qd(model.joints{i});
                end
            end

            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                c_X_j = eye(6);
                parent = link;
            end

            while parent ~= 0
                j = parent;
                % v_i x si*qd_i
                if(j == 1)
                    c_X_j = c_X_j*model.i_X_pi{j};
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                else
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end

            if(coord == 0)
                %Transform to global
                jacdot = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jacdot;
            else
            end
        end

        function jacdot = jacob_pt_dot_spatial(model,link,pt,qd,coord)
            % JACOB_PT_DOT_SPATIAL compute jacobian dot for a point in the body
            %   jacdot = JACOB_PT_DOT_SPATIAL(model,link,pt,qd,coord)
            %
            %   This function computes the derivative of the Jacobian for a
            %   point w.r.t a link coordinate frame based on spatial vector
            %   representation. It is currently incorrect and this function
            %   will be revisited once jacob_dot_spatial is fixed. TODO:
            %   revisit the implementation and fix it with
            %   jacob_dot_spatial function.

            if ~exist('coord','var')
                coord = 0; %world coordinate by default
            end
            % Initalize Jacobian
            jacdot = zeros(6,model.N_joint);

            % compute the velocity of each body
            %% RNEA outward recursion
            for i=1:model.NB
                p = model.parent{i};
                if p == 0
                    v{i} = model.phi{i}*qd(model.joints{i});
                else
                    v{i} = model.i_X_pi{i}*v{p} + model.phi{i}*qd(model.joints{i});
                end
            end

            % Recursively build jacobian
            if link > model.NB % for virtual body (frames)
                parent = model.parent{link};
                c_X_j = model.i_X_pi{link};
            else
                %                 c_X_j = eye(6);
                c_X_j = SpatialRobotModel.xlt(pt);
                parent = link;
            end

            while parent ~= 0
                j = parent;
                % v_i x si*qd_i
                if(j==1)
                    c_X_j = c_X_j*model.i_X_pi{j};
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                else
                    jacdot(:,model.joints{j}) = c_X_j*model.crossSpv(v{j})*model.phi{j};
                    c_X_j = c_X_j*model.i_X_pi{j};
                end
                parent = model.parent{j};
            end

            if(coord == 0)
                %Transform to global
                jacdot = [model.R0{link} zeros(3,3); zeros(3,3) model.R0{link}]*jacdot;
            else
            end
        end

        function updateModel(model, q, display)
            % UPDATEMODEL update the model with the joint configurations
            %   UPDATEMODEL(model, q, display)
            %
            %   This function updates model kinematics and dynamics with
            %   the provided joint configuration vector. display is set to
            %   0 by default and it enables the default visualization for
            %   frames and inertia without using 3d mesh files. Note that
            %   is use_mex field is true, it calls the mex function for
            %   speed.
            %
            %   See also updateModelFast.

            if ~exist('display','var')
                display = 0; %off by default
            end

            if ~exist('q','var')
                model.q_p = {zeros(6,1)};
                model.q = [model.q_p; num2cell(model.q_a)];
                q = [zeros(6,1); model.q_a];
            else % if q exists
                if numel(q) == model.NJ
                    model.q_p = {zeros(6,1)};
                    model.q_a = q;
                    model.q = [model.q_p; num2cell(model.q_a)];
                    q = [zeros(6,1); model.q_a];
                elseif numel(q) == model.N_joint
                    model.q_p = q(1:6);
                    model.q_a = q(7:end);
                    model.q = [model.q_p; num2cell(model.q_a)];
                else
                    error('Invalid size of q!')
                end
            end

            if(model.use_mex) % using mex function
                model.model_param = updateModelMex_mex(model.model_param, q);
                model.struct2model();
            else

                % base body (floating)
                KBody = 1;
                q_base = model.q{KBody};
                if(model.base_motion)
                    model.pi_p_i{KBody}   = q_base(1:3); % position of the joint
                    q_ = q_base(4:6);
                    if(model.use_mex)
                        R_base = rotz_mex(q_(3))*roty_mex(q_(2))*rotx_mex(q_(1));
                    else
                        R_base = SpatialRobotModel.rotz(q_(3))*SpatialRobotModel.roty(q_(2))*SpatialRobotModel.rotx(q_(1));
                    end
                    model.pi_R_i{KBody} = R_base; % rotation of the joint
                end

                %% Setup all neccecaries for model structs
                for i=1:model.N_ext
                    %                 j  = model;
                    pred = model.parent{i};

                    % Setup Coordinates

                    if i==1 % base bodies
                        % SpatialRobotModel.translation and rotation
                        model.i_XL_pi{i} = eye(6);
                    else % non-base body
                        model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -model.pi_R_i{i}'*SpatialRobotModel.skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                    end

                    if i == 1 % floating-base
                        R = eye(3);
                        model.X_J{i} = [model.pi_R_i{i}' zeros(3); zeros(3) model.pi_R_i{i}'];

                    elseif i > model.NB %end-effector frame (no transform due to joint rotation)
                        R = eye(3);
                        model.X_J{i} = eye(6);

                    else %1 < i < model.NB
                        if ~isempty(model.joints{i})
                            q_ = model.q{model.joints{i}-5};
                        else
                            q_ = 0;
                        end

                        R = SpatialRobotModel.angvec2r(-q_, model.phi{i}(1:3)); % equivalent expression
                        model.X_J{i} = [R zeros(3); zeros(3) R];
                    end

                    % 6x6 Adj transform
                    model.i_X_pi{i} = model.X_J{i} * model.i_XL_pi{i};

                    % 4x4 HT
                    if i == 1 % at floating-joint
                        model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1];
                    else
                        model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1]*[ R' [0;0;0] ; 0 0 0 1];
                    end

                    % Initialize inertia
                    I_com_bar = model.I_com{i} + model.m{i} * SpatialRobotModel.crossSkew(model.i_p_com{i}) * SpatialRobotModel.crossSkew(model.i_p_com{i})';

                    % inertia w.r.t the i-th link coordinate
                    model.I{i} = [I_com_bar model.m{i}*SpatialRobotModel.crossSkew(model.i_p_com{i}) ; model.m{i} * SpatialRobotModel.crossSkew(model.i_p_com{i})' model.m{i} * eye(3);];

                    % Initialize global transforms
                    if model.parent{i} == 0
                        model.T0{i} = real(model.pi_T_i{i});
                        model.R0{i} = model.pi_T_i{i}(1:3,1:3);
                        model.i_X_0{i} = model.i_X_pi{i};

                        if(display == 1)
                            % display coordinate frames at the base
                            if(model.display_frames ==1 )
                                SpatialRobotModel.drawFrame(model.T0{i}, model.frame_scale, 2);
                            end
                        end
                    else

                        model.T0{i} = model.T0{pred} * model.pi_T_i{i};
                        model.R0{i} = model.R0{pred} * model.pi_T_i{i}(1:3,1:3);
                        model.i_X_0{i} = model.i_X_pi{i}*model.i_X_0{pred};

                        if(display == 1)

                            % display segments (linkages) - by default
                            SpatialRobotModel.plotRows([model.T0{i}*[0 0 0 1]' model.T0{pred}*[0 0 0 1]'],model.color_segment, model.width_segment);

                            % display coordinate frames
                            if(model.display_frames ==1 )
                                SpatialRobotModel.drawFrame(model.T0{i}, model.frame_scale, 2);
                            end
                        end
                    end
                    % for jacobian transform
                    model.R0_6{i} = [model.R0{i} zeros(3) ; zeros(3) model.R0{i}];

                    if(display == 1)

                        if i<=model.NB

                            if(model.display_mass == 1)
                                % Plot Masses (only for the bodies that are moved by joints)
                                SpatialRobotModel.scatterRows(model.T0{i}*[model.i_p_com{i} ;1]);
                            end

                            if(model.display_inertia == 1)

                                [V, D] = eig(model.I_com{i});

                                S = (model.m{i}/5*(ones(3)-eye(3)))\diag(D);
                                S = sqrt(abs(S))/5*(model.m{i})^.25;
                                v1 = V(:,1)*S(1); v2 = V(:,2)*S(2); v3 = V(:,3)*S(3);

                                theta = [0:.1:2*pi 0];
                                ct = cos(theta); st = sin(theta); onest = 0*ct+1;

                                % Plot inertias
                                SpatialRobotModel.plotRows(model.T0{i}*[(v1*ct+v2*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                                SpatialRobotModel.plotRows(model.T0{i}*[(v1*ct+v3*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                                SpatialRobotModel.plotRows(model.T0{i}*[(v3*ct+v2*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                            end
                        end
                    end
                end

            end
        end

        function updateModelFast(model, q, display)
            % UPDATEMODELFAST update the model with the joint configurations only for FK purposes
            %   UPDATEMODELFAST(model, q, display)
            %
            %   As opposed to updateModel This function only updates
            %   kinematic model for faster FK and IK purposes.
            %
            %   See also updateModel.

            if ~exist('display','var')
                display = 0; %off by default
            end

            if ~exist('q','var')
                model.q_p = {zeros(6,1)};
                model.q = [model.q_p; num2cell(model.q_a)];
                q = [zeros(6,1); model.q_a];
            else % if q exists
                if numel(q) == model.NJ
                    model.q_p = {zeros(6,1)};
                    model.q_a = q;
                    model.q = [model.q_p; num2cell(model.q_a)];
                    q = [zeros(6,1); model.q_a];
                elseif numel(q) == model.N_joint
                    model.q_p = q(1:6);
                    model.q_a = q(7:end);
                    model.q = [model.q_p; num2cell(model.q_a)];
                else
                    error('Invalid size of q!')
                end
            end

            if(model.use_mex) % using mex function
                model.model_param = updateModelMex_mex(model.model_param, q);
                model.struct2model();
            else

                % base body (floating)
                KBody = 1;
                q_base = model.q{KBody};
                if(model.base_motion)
                    model.pi_p_i{KBody}   = q_base(1:3); % position of the joint
                    q_ = q_base(4:6);
                    if(model.use_mex)
                        R_base = rotz_mex(q_(3))*roty_mex(q_(2))*rotx_mex(q_(1));
                    else
                        R_base = SpatialRobotModel.rotz(q_(3))*SpatialRobotModel.roty(q_(2))*SpatialRobotModel.rotx(q_(1));
                    end
                    model.pi_R_i{KBody} = R_base; % rotation of the joint
                end

                if(display == 1)
                    if(0)
                        figure(1);
                        clf;
                        hold on
                        axis square

                        x = model.pi_p_i{1}(1);
                        y = model.pi_p_i{1}(2);
                        z = model.pi_p_i{1}(3);

                        xlabel('x');
                        ylabel('y');
                        zlabel('z');

                        hold on;
                        grid on;

                    else
                        figure(1);
                        if(refresh)
                            clf;
                        end

                        % older version than 2014b
                        if((strcmp(version, '8.4.0.150421 (R2014b)') || strcmp(version, '8.5.0.197613 (R2015a)')) || strcmp(version, '8.6.0.267246 (R2015b)') == 0)
                            set(gcf,'NumberTitle', 'off', ...
                                'BackingStore','off', ...
                                'MenuBar','default', ...
                                'Color',[1 1 1], ...
                                'Renderer','OpenGL');

                            axis equal;

                            xlabel('x');
                            ylabel('y');
                            zlabel('z');

                            hold on;
                            grid on;

                            % change some properties
                            set(gca, 'DrawMode',   'fast', ...
                                'Color',     [1 1 1], ...
                                'XColor',    [0 0 0], ...
                                'YColor',    [0 0 0]);
                        else  % for 2014b
                            set(gcf,'NumberTitle', 'off', ...
                                'BackingStore','off', ...
                                'MenuBar','default', ...
                                'Color',[1 1 1]);

                            axis equal;

                            xlabel('x');
                            ylabel('y');
                            zlabel('z');

                            hold on;
                            grid on;

                        end

                    end
                end

                %% Setup all neccecaries for model structs
                for i=1:model.N_ext
                    pred = model.parent{i};

                    % Setup Coordinates

                    if i==1 % base bodies
                        % SpatialRobotModel.translation and rotation
                        model.i_XL_pi{i} = eye(6);
                    else % non-base body
                        model.i_XL_pi{i} = [model.pi_R_i{i}' zeros(3); -model.pi_R_i{i}'*SpatialRobotModel.skew(model.pi_p_i{i}) model.pi_R_i{i}'];
                    end

                    if i == 1 % floating-base
                        R = eye(3);
                        model.X_J{i} = [model.pi_R_i{i}' zeros(3); zeros(3) model.pi_R_i{i}'];

                    elseif i > model.NB %end-effector frame (no transform due to joint rotation)
                        R = eye(3);
                        model.X_J{i} = eye(6);

                    else %1 < i < model.NB
                        if ~isempty(model.joints{i})
                            q_ = model.q{model.joints{i}-5};
                        else
                            q_ = 0;
                        end

                        R = SpatialRobotModel.angvec2r(-q_, model.phi{i}(1:3)); % equivalent expression
                        model.X_J{i} = [R zeros(3); zeros(3) R];
                    end

                    % 6x6 Adj transform
                    model.i_X_pi{i} = model.X_J{i} * model.i_XL_pi{i};

                    % 4x4 HT
                    if i == 1 % at floating-joint
                        model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1];
                    else
                        model.pi_T_i{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1]*[ R' [0;0;0] ; 0 0 0 1];
                    end

                    % Initialize global transforms
                    if model.parent{i} == 0
                        model.T0{i} = real(model.pi_T_i{i});
                        model.R0{i} = model.pi_T_i{i}(1:3,1:3);
                        %                     model.R0{i}
                        model.i_X_0{i} = model.i_X_pi{i};

                        if(display == 1)
                            % display coordinate frames at the base
                            if(model.display_frames ==1 )
                                SpatialRobotModel.drawFrame(model.T0{i}, model.frame_scale, 2);
                            end
                        end
                    else

                        model.T0{i} = model.T0{pred} * model.pi_T_i{i};
                        model.R0{i} = model.R0{pred} * model.pi_T_i{i}(1:3,1:3);
                        model.i_X_0{i} = model.i_X_pi{i}*model.i_X_0{pred};

                        if(display == 1)

                            % display segments (linkages) - by default
                            SpatialRobotModel.plotRows([model.T0{i}*[0 0 0 1]' model.T0{pred}*[0 0 0 1]'],model.color_segment, model.width_segment);

                            % display coordinate frames
                            if(model.display_frames ==1 )
                                %                             drawframe(model.T0{i}, model.frame_scale);
                                SpatialRobotModel.drawFrame(model.T0{i}, model.frame_scale, 2);
                            end
                        end
                    end
                    % for jacobian transform
                    model.R0_6{i} = [model.R0{i} zeros(3) ; zeros(3) model.R0{i}];

                    if(display == 1)

                        if i<=model.NB

                            if(model.display_mass == 1)
                                % Plot Masses (only for the bodies that are moved by joints)
                                SpatialRobotModel.scatterRows(model.T0{i}*[model.i_p_com{i} ;1]);
                            end

                            if(model.display_inertia == 1)

                                [V, D] = eig(model.I_com{i});
                                %D=sqrt(diag(D));
                                %D=D/10;
                                S = (model.m{i}/5*(ones(3)-eye(3)))\diag(D);
                                %                             S = sqrt(S)/5*(model.m{i})^.25;
                                S = sqrt(abs(S))/5*(model.m{i})^.25;
                                v1 = V(:,1)*S(1); v2 = V(:,2)*S(2); v3 = V(:,3)*S(3);

                                theta = [0:.1:2*pi 0];
                                ct = cos(theta); st = sin(theta); onest = 0*ct+1;

                                % Plot inertias
                                SpatialRobotModel.plotRows(model.T0{i}*[(v1*ct+v2*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                                SpatialRobotModel.plotRows(model.T0{i}*[(v1*ct+v3*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                                SpatialRobotModel.plotRows(model.T0{i}*[(v3*ct+v2*st)+ model.i_p_com{i}*onest;onest], model.color_inertia, model.width_inertia)
                            end
                        end
                    end
                end
            end
        end

        function T0 = FKin(model,q)
            % FKIN compute forward kinematics of the model and return transform matrices
            %   T0 = FKin(model,q)
            %
            %   This function does forward kinematics computation which is
            %   only a subset of what updateModel does and returns a cell
            %   array of tranform matrices for all the link coordinate
            %   frames.
            %
            %   See also updateModel, updateModelFast.

            % cell array initialization
            T0 = cell(model.N_ext,1);

            %= Setup all neccecaries for model structs
            for i=1:model.N_ext
                pred = model.parent{i};

                if i == 1 % floating-base
                    R = eye(3);

                elseif i > model.NB %end-effector frame (no transform due to joint rotation)
                    R = eye(3);

                else %1 < i < model.NB
                    R = SpatialRobotModel.angvec2r(-q(i-1), model.phi{i}(1:3)); % equivalent expression
                end

                % 4x4 HT
                if i == 1 % at floating-joint
                    pi_T_i_tmp{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1];
                else
                    pi_T_i_tmp{i} =  [ model.pi_R_i{i} model.pi_p_i{i} ; 0 0 0 1]*[ R' [0;0;0] ; 0 0 0 1];
                end

                % Initialize global transforms
                if model.parent{i} == 0
                    T0{i} = pi_T_i_tmp{i};
                else
                    T0{i} = T0{pred} * pi_T_i_tmp{i};
                end
            end
        end

        function tau = inverseDynamics(model,qdd,qd,f_ext,grav0)
            % INVERSEDYNAMICS compute inverse dynamics for a fixed-base robot
            %   tauTotal = INVERSEDYNAMICS(model,qdd,qd,f_ext,grav0)
            %
            %   This function computes inverse dynamics for a fixed base
            %   robot model. It returns a torque vector (NJ x 1).
            %
            %   See also inverseDynamicsFB, inverseDynamicsAll.

            if ~exist('grav','var')
                grav0 = model.ag_0; % default gravity constant from model
            end

            if ~exist('f_ext','var')
                % if external force not given, we assume zero external
                % force
                f_ext = cell(model.NB,1);
                for i = 1:model.NB
                    f_ext{i} = zeros(6,1);
                end
            end

            tauTotal = zeros(model.N_joint,1);

            %% RNEA outward recursion
            for i=1:model.NB
                %                 r = model;
                p = model.parent{i};
                if p == 0 % at base body
                    model.v{i}   = model.phi{i}*qd(model.joints{i});
                    model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                    model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                    model.f{i}   = model.I{i}*model.a{i} +...
                        SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                else
                    if ~isempty(model.joints{i})
                        qd_joint = qd(model.joints{i});
                        qdd_joint = qdd(model.joints{i});
                    else % for fixed joints
                        qd_joint = 0;
                        qdd_joint = 0;
                    end
                    model.v{i} = model.i_X_pi{i}*model.v{p} + model.phi{i}*qd_joint;
                    model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd_joint;
                    model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i} + model.phi{i}*qdd_joint;
                    model.f{i}   = model.I{i}*model.a{i} + ...
                        SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                end
            end
            %% rnea inward recursion
            for i=model.NB:-1:1
                p = model.parent{i};
                if ~isempty(model.joints{i}) % only for revolute joints
                    tauTotal(model.joints{i}) = model.phi{i}'*model.f{i};
                end
                if p>0
                    model.f{p} = model.f{p}+model.i_X_pi{i}'*model.f{i};
                end
            end

            % extract actuated joint torque vector
            tau = tauTotal(7:end);

        end

        function [tau, tauAllAxes] = inverseDynamicsAll(model,qdd,qd,payload_only,f_ext,grav0)
            % INVERSEDYNAMICSALL compute inverse dynamics for a fixed-base robot (consider payload and offaxis loading)
            %   [tau, tauAllAxes] = INVERSEDYNAMICSALL(model,qdd,qd,payload_only,f_ext,grav0)
            %
            %   This function returns joint torques as well as off-axis
            %   torques from inverse dynamics computation for a fixed-base
            %   robot. It also considers payload defined at end-effector
            %   frames.
            %
            %   See also inverseDynamics.

            if ~exist('grav','var')
                grav0 = model.ag_0; % default gravity constant from model
            end

            if ~exist('f_ext','var')
                % if external force not given, we assume zero external
                % force
                f_ext = cell(model.N_ext,1);
                for i = 1:model.N_ext
                    f_ext{i} = zeros(6,1);
                end
            end

            if ~exist('payload_only','var')
                % by default, the ID considers the inertia&mass in all the
                % bodies including payload
                payload_only = false;
            end

            % torques about axis of rotation
            tauTotal = zeros(model.N_joint,1);

            % torques about all axes
            tau = zeros(model.NJ,1);
            tauAllAxes = zeros(model.NJ,3);

            %% RNEA outward recursion
            for i=1:model.N_ext

                if (i <= model.NB) % bodies with joints

                    if(payload_only == false)
                        p = model.parent{i};
                        if p == 0 % at base body
                            model.v{i}   = model.phi{i}*qd(model.joints{i});
                            model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                            model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                            model.f{i}   = model.I{i}*model.a{i} +...
                                SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                        else
                            if ~isempty(model.joints{i}) % only for revolute joints
                                qd_joint = qd(model.joints{i});
                                qdd_joint = qdd(model.joints{i});
                            else % for fixed joints
                                qd_joint = 0;
                                qdd_joint = 0;
                            end
                            model.v{i} = model.i_X_pi{i}*model.v{p} + model.phi{i}*qd_joint;
                            model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd_joint;
                            model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i} + model.phi{i}*qdd_joint;
                            model.f{i}   = model.I{i}*model.a{i} + ...
                                SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                        end

                        % zero out the mass and inertia in the robot links
                    elseif(payload_only == true)
                        I_zero = 0*model.I{i};

                        p = model.parent{i};
                        if p == 0 % at base body
                            model.v{i}   = model.phi{i}*qd(model.joints{i});
                            model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                            model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i} + model.phi{i}*qdd(model.joints{i});
                            model.f{i}   = I_zero*model.a{i} +...
                                SpatialRobotModel.crossSpf(model.v{i})*I_zero*model.v{i} - f_ext{i};
                        else
                            if ~isempty(model.joints{i}) % only for revolute joints
                                qd_joint = qd(model.joints{i});
                                qdd_joint = qdd(model.joints{i});
                            else % for fixed joints
                                qd_joint = 0;
                                qdd_joint = 0;
                            end
                            model.v{i} = model.i_X_pi{i}*model.v{p} + model.phi{i}*qd_joint;
                            model.chi{i} = SpatialRobotModel.crossSpv(model.v{i})*model.phi{i}*qd_joint;
                            model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i} + model.phi{i}*qdd_joint;
                            model.f{i}   = I_zero*model.a{i} + ...
                                SpatialRobotModel.crossSpf(model.v{i})*I_zero*model.v{i} - f_ext{i};
                        end
                    end

                else % frames with payload

                    p = model.parent{i};
                    if p == 0 % at base body
                        model.v{i}   = 0;
                        model.chi{i} = 0;
                        model.a{i}   = model.i_X_pi{i}*(-grav0) + model.chi{i};
                        model.f{i}   = model.I{i}*model.a{i} - f_ext{i};
                    else
                        model.v{i} = model.i_X_pi{i}*model.v{p};
                        model.chi{i} = 0;
                        model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i};
                        model.f{i}   = model.I{i}*model.a{i} + ...
                            SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
                    end

                end

            end

            %% rnea inward recursion
            for i=model.N_ext:-1:1
                if (i <= model.NB)
                    f_i_all = model.f{i}; % wrench transmitted between i and p(i) bodies
                    if ~isempty(model.joints{i}) % only for revolute joints
                        tauTotal(model.joints{i}) = model.phi{i}'*model.f{i};
                        if (i > 1)
                            % take the torques for all three axes
                            tauAllAxes(model.joints{i}-6, :) = f_i_all(1:3)';
                        end
                    end
                end

                p = model.parent{i};
                if p>0
                    model.f{p} = model.f{p}+model.i_X_pi{i}'*model.f{i};
                end
            end

            % get the actuated joint torque
            tau = tauTotal(7:end);
        end

        function [a0, tau] = inverseDynamicsFB(model,qdd,qd,f_ext)
            % INVERSEDYNAMICSFB compute inverse dynamics for a floating-base robot
            %   [a0_method3, tau_method3] = INVERSEDYNAMICSFB(model,qdd,qd,f_ext,grav0)
            %
            %   This function computes inverse dynamics for a floating-base
            %   robot model. It returns the spatial acceleration at the
            %   base as well as the actuated joint torques.
            %
            %   See also inverseDynamics, inverseDynamicsAll.

            % initialization for the floating base
            model.a{1} = -model.ag_0;
            model.I_C{1} = model.I{1};
            model.v{1} = model.phi{1}*qd(1:6);
            %             model.v{1} = qd(1:6); % because of this multiplication qd is reversed again!
            model.f{1} = model.I{1}* model.a{1} + ...
                SpatialRobotModel.crossSpf(model.v{1}) * model.I{1}*model.v{1} -...
                f_ext{1};

            % inverse dynamics with floating base accel ==0
            % outward rnea
            %             qdd = [zeros(6,1) ; qddJoints];
            for i = 2:model.NB
                p = model.parent{i};
                model.v{i} = model.i_X_pi{i}*model.v{p}+model.phi{i}*qd(model.joints{i});
                model.chi{i} = model.crossSpv(model.v{i})*model.phi{i}*qd(model.joints{i});
                model.a{i}   = model.i_X_pi{i}*model.a{p} + model.chi{i}+model.phi{i}*qdd(model.joints{i});
                model.f{i}   = model.I{i}*model.a{i} + ...
                    SpatialRobotModel.crossSpf(model.v{i})*model.I{i}*model.v{i} - f_ext{i};
            end
            % inward rnea
            tauTotal = zeros(model.N_joint,1);
            for i=model.NB:-1:2
                p = model.parent{i};
                tauTotal(model.joints{i}) =  model.phi{i}'*model.f{i};
                model.f{p} = model.f{p}+model.i_X_pi{i}'*model.f{i};
            end

            % Composite rigid body recursion for I_C of floating base
            for i = 1:model.NB
                model.I_C{i} = model.I{i};
            end
            for i = model.NB:-1:2
                p = model.parent{i};
                %             r = model(i);
                model.I_C{p} = model.I_C{p} + model.i_X_pi{i}' * model.I_C{i} * model.i_X_pi{i};
            end

            a0 = (-model.I_C{1})\model.f{1};
            model.a0{1} = a0;
            % Add in torques due to floating base acceleration
            for i = 2:model.NB
                p = model.parent{i};
                model.a0{i} = model.i_X_pi{i} * model.a0{p};
                % Add in additional torques
                tauTotal(model.joints{i}) = tauTotal(model.joints{i}) + model.phi{i}'*model.I_C{i}*model.a0{i};
            end
            tau = tauTotal(7:end);
        end

        function Hout = computeH (model)
            % COMPUTEH compute the joint space mass matrix
            %   Hout = COMPUTEH (model)
            %
            %   This function computes an inertia matrix (NJ x NJ) for
            %   actuated joints.
            %
            %   See also computeHAll.

            if(model.use_mex)
                H = computeHMex_mex(model.model_param);
            else
                %= Calculate H
                H = zeros(model.N_joint, model.N_joint);

                for i = 1:model.NB
                    I_C{i} = model.I{i};
                end

                for i = model.NB:-1:1
                    F = I_C{i}*model.phi{i};
                    %model(i).I_C
                    %fprintf('Link (%d,%d), (%f,%f,%f,%f,%f,%f)\NB',i,i,F(1),F(2),F(3),F(4),F(5),F(6));
                    if model.parent{i} ~= 0
                        pred = model.parent{i};
                        % compute composite rigid body inertia
                        I_C{pred} = I_C{pred} + model.i_X_pi{i}'*I_C{i}*model.i_X_pi{i};
                    end

                    if ~isempty(model.joints{i})
                        H(model.joints{i},model.joints{i}) = model.phi{i}'*F; % diagonal components
                        j = i;
                        while model.parent{j} ~= 0 % coupling components
                            F = model.i_X_pi{j}'*F;
                            j = model.parent{j};
                            if ~isempty(model.joints{j})
                                %fprintf('Link (%d,%d), (%f,%f,%f,%f,%f,%f)\NB',i,j,F(1),F(2),F(3),F(4),F(5),F(6));
                                H(model.joints{i},model.joints{j}) = F' * model.phi{j};
                                H(model.joints{j},model.joints{i}) = H(model.joints{i},model.joints{j})';
                            end
                        end
                    end
                end
            end

            % extract inertia matrix for actuated joints
            Hout = H(7:end, 7:end);

        end

        function Hout = computeHAll(model)
            % COMPUTEHALL compute the joint inertia matrix including payload
            %   Hout = COMPUTEHALL(model)
            %
            %   This function computes a joint inertia matrix considering
            %   the payload defined at the end-effector frame.
            %
            %   See also computeH.

            if(model.use_mex)
                H = computeHAllMex_mex(model.model_param);
            else
                H = zeros(model.NB,model.NB);
                I_C = cell(model.N_ext, 1);

                for i = 1:model.N_ext
                    I_C{i} = model.I{i};
                end

                for i = model.N_ext:-1:1

                    % this update now includes the payloads
                    if model.parent{i} ~= 0
                        pred = model.parent{i};
                        I_C{pred} = I_C{pred} + model.i_X_pi{i}'*I_C{i}*model.i_X_pi{i};
                    end

                    % update only for the joint-associated bodies
                    if(i <= model.NB && ~isempty(model.joints{i}))
                        F = I_C{i}*model.phi{i};
                        H(model.joints{i},model.joints{i}) = model.phi{i}'*F;

                        j = i;
                        while model.parent{j} ~= 0
                            F = model.i_X_pi{j}'*F;
                            j = model.parent{j};
                            if ~isempty(model.joints{j})
                                %fprintf('Link (%d,%d), (%f,%f,%f,%f,%f,%f)\NB',i,j,F(1),F(2),F(3),F(4),F(5),F(6));
                                H(model.joints{i},model.joints{j}) = F' * model.phi{j};
                                H(model.joints{j},model.joints{i}) = H(model.joints{i},model.joints{j})';
                            end
                        end
                    end
                end
            end

            % extract inertia matrix for actuated joints
            Hout = H(7:end, 7:end);
        end

        function out = computeCandG(model,qd,grav)
            % COMPUTECANDG compute C and G component in the rigid body system
            %   out = COMPUTECANDG(model,qd,grav)
            %
            %   This function computes coriolis, centrifugal and gravity
            %   related torques.
            %
            %   See also computeCandGAll.

            if ~exist('grav','var')
                grav = model.ag_0; % compute C and G both by default
            end

            qdd = zeros(model.N_joint,1);

            if(model.use_mex)
                f_ext_tmp = zeros(model.N_ext, 6);
                out = inverseDynamicsMex_mex(model.model_param, qdd, qd, f_ext_tmp, grav);
                out = out(7:end);
            else
                f_ext_tmp = cell(model.NB,1);
                for i = 1:model.NB
                    f_ext_tmp{i} = zeros(6,1);
                end
                out = model.inverseDynamics(qdd, qd, f_ext_tmp, grav);
            end
        end

        function [tq, tqAll] = computeCandGAll(model, qd, payload_only, grav)
            % COMPUTECANDGALL compute C and G component in the rigid body system (including payload)
            %   [tq, tqAll] = COMPUTECANDGALL(model, qd, payload_only, grav)
            %
            %   This function computes coriolis, centrifugal and gravity
            %   related torques considering the payload defined at the
            %   end-effector frame.
            %
            %   See also computeCandG.

            if ~exist('payload_only','var')
                % by default, the ID considers the inertia&mass in all the
                % bodies including payload
                payload_only = false;
            end

            if ~exist('grav','var')
                grav = model.ag_0; % compute C and G both by default
            end

            qdd = zeros(model.N_joint,1);

            if(model.use_mex)
                f_ext_tmp = zeros(model.N_ext, 6);
                [tq, tqAll] = inverseDynamicsAllMex_mex(model.model_param, qdd, qd, payload_only, f_ext_tmp, grav);
            else
                f_ext_tmp = cell(model.N_ext,1);
                for i = 1:model.N_ext
                    f_ext_tmp{i} = zeros(6,1);
                end
                [tq, tqAll] = model.inverseDynamicsAll(qdd, qd, payload_only, f_ext_tmp, grav);
            end
        end

        function out = computeG(model,grav)
            % COMPUTEG compute Gravity torque component in the rigid body system
            %   out = COMPUTEG(model,grav)
            %
            %   This function computes gravity compensation torques.
            %
            %   See also computeGAll.

            if ~exist('grav','var')
                grav = model.ag_0; % use gravity from the model definition
            end

            qdd = zeros(model.N_joint,1);
            qd = zeros(model.N_joint,1);

            if(model.use_mex)
                f_ext_zero = zeros(model.N_ext, 6);
                out = inverseDynamicsMex_mex(model.model_param, qdd, qd, f_ext_zero, grav);
                out = out(7:end);
            else
                f_ext_zero = cell(model.NB,1);
                for i = 1:model.NB
                    f_ext_zero{i} = zeros(6,1);
                end
                out = model.inverseDynamics(qdd, qd, f_ext_zero, grav);
            end
        end

        function [tq, tqAll] = computeGAll(model, qdd, qd, payload_only, grav)
            % COMPUTEGALL compute Gravity torque component in the rigid body system (including payload)
            %   [tq, tqAll] = COMPUTEGALL(model, qdd, qd, payload_only, grav)
            %
            %   This function computes gravity compensation torques
            %   considering the payload defined at the end-effector frame.
            %   With payload_only set to true, it computes gravity
            %   compensation torques required only for payload.
            %
            %   See also computeG.

            if ~exist('payload_only','var')
                % by default, the ID considers the inertia&mass in all the
                % bodies including payload
                payload_only = false;
            end

            if ~exist('grav','var')
                grav = model.ag_0; % use gravity from the model definition
            end

            if ~exist('qdd','var')
                qdd = zeros(model.N_joint,1);
            end

            if ~exist('qd','var')
                qd = zeros(model.N_joint,1);
            end


            if(model.use_mex)
                f_ext_zero = zeros(model.N_ext, 6);
                [tq, tqAll] = inverseDynamicsAllMex_mex(model.model_param, qdd, qd, payload_only, f_ext_zero, grav);
            else

                f_ext_zero = cell(model.N_ext,1);
                for i = 1:model.N_ext
                    f_ext_zero{i} = zeros(6,1);
                end
                [tq, tqAll] = model.inverseDynamicsAll(qdd, qd, payload_only, f_ext_zero, grav);
            end
        end

        %% inertia related computations

        function Pcom = getCoM(model)
            % GETCOM compute the CoM
            %   Pcom = GETCOM(model)
            %
            %   This function returns a 3x1 position vector for
            %   the Center-of-Mass (CoM) in the robot system.

            c = 0;
            mb = 0;
            for i = 1:model.NB
                mb = mb + model.m{i};
                c_i = model.T0{i}*[model.i_p_com{i} ;1];
                c = c + model.m{i}*c_i;
            end
            Pcom = c(1:3)/mb;
        end

        function Jcom = getJcom(model)
            % GETJCOM compute COM jacobian
            %   Jcom = GETJCOM(model)
            %
            %   This function computes the Jacobian at the CoM in the robot
            %   system.
            %
            %   See also getJcom_dot.

            % Initalize Jacobian
            jac_com = zeros(6,model.N_joint);
            mb = 0;

            for link=1:model.NB
                jac_mi = model.jacob_pt(link, model.i_p_com{link});
                mb = mb + model.m{link};
                jac_com = jac_com + model.m{link}*jac_mi;
            end
            Jcom = jac_com/mb;
        end

        function Jcom_dot = getJcom_dot(model,qd)
            % GETJCOM_DOT compute COM jacobian dot
            %   Jcom_dot = GETJCOM_DOT(model,qd)
            %
            %   This function returns the derivative of Jacobian at the CoM
            %   in the robot system.
            %
            %   See also getJcom.

            % Initalize Jacobian
            jac_com_dot = zeros(6,model.N_joint);

            mb = 0;
            for link=1:model.NB
                jac_mi_dot = model.jacob_spatial_pt_dot(link, model.i_p_com{link}, qd);
                mb = mb + model.m{link};
                jac_com_dot = jac_com_dot + model.m{link}*jac_mi_dot;
            end
            Jcom_dot = jac_com_dot/mb;
        end

        function M_total = getMass(model)
            % GETMASS compute the whole-body mass
            %   M_total = GETMASS(model)
            %
            %   See also getMassAll.

            % sum up individual link mass
            mb = 0;
            for i = 1:model.NB
                mb = mb + model.m{i};
            end

            M_total = mb;
        end

        function M_total = getMassAll(model)
            % GETMASSALL compute the whole-body mass (including payload)
            %   M_total = GETMASSALL(model)
            %
            %   See also getMass.

            % sum up individual link mass as well as payload mass
            mb = 0;
            for i = 1:model.N_ext
                mb = mb + model.m{i};
            end

            M_total = mb;
        end

        function I_C_0 = getCompositeInertia(model)
            % CompositeInertia compute the composite rigid body inertia for the floating-base
            %   I_C_0 = GETCOMPOSITEINERTIA(model)

            % Composite rigid body recursion for I_C of floating base
            for i = 1:model.NB
                I_C{i} = model.I{i};
            end

            I_C_0 = 0;
            for i = model.NB:-1:1
                p = model.parent{i};
                if p > 0
                    I_C{p} = I_C{p} + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                else %p=0 at floating base
                    I_C_0 = I_C_0 + model.i_X_pi{i}' * I_C{i} * model.i_X_pi{i};
                end
            end
        end

        %% pose metrics related computations

        function out = optimize_performance_index(model, q_a, select_method, u)
            % optimize_performance_index compute performance indices of jacobian

            % 1) manipulability
            % 2) condition number
            % 3) directional velocity transmission ratio
            if ~exist('u','var')
                u = [1; 0; 0]; % unit vector
            end

            % get a temporary model object
            modeltemp = model;
            modeltemp.updateModelFast([zeros(6,1); q_a], 0);

            % compute the jacobian at end-point
            J =  model.active(model.jacob(model.ind_H));

            if(select_method == 1)
                %-- compute velocity manipulability
                out = sqrt(det(J * J'));

            elseif(select_method == 2)
                %-- compute force manipulability
                out = 1/sqrt(det(J * J'));

            elseif(select_method == 3)
                %-- compute condition number
                out = cond(J * J');

            elseif(select_method == 4)
                %-- compute the velocity transmission ratio along u
                out = 1/sqrt(u'/(J * J')*u);

            elseif(select_method == 5)
                %-- compute the force transmission ratio along u
                out = sqrt(u'/(J * J')*u);
            end
        end

        %% visualization related functions

        function initAnimation(model, isModelGhost)
            % INITANIMATION initialize a robot graphical model for visualization
            %   INITANIMATION(model)
            %
            %   This function loads STL mesh files and initializes the
            %   figure with a robot animation model displayed.
            %
            %   See also updateAnimationModel.

            if ~exist('isModelGhost','var')
                isModelGhost = false;
            end

            model.Nmeshes = length(model.mesh_name); % # of meshes

            % initialize variables
            model.fv = cell(model.Nmeshes,1);
            model.h_fv = cell(model.Nmeshes,1);
            model.T_fv = cell(model.Nmeshes,1);
            model.T_fv_col = cell(model.Nmeshes,1);
            model.vdata.vmins = cell(model.Nmeshes,1);
            model.vdata.vmaxs = cell(model.Nmeshes,1);
            model.vdata.vmlens = cell(model.Nmeshes,1);

            if(model.display_frames == 1)
                model.h_frame = cell(size(model.N_ext,2),1);
            end

            %= set display settings for ghost model
            if (isModelGhost == true)
                model.display_frames = 0; % turn on/off display coordinate frames
                model.display_inertia = 0; % turn on/off display CoM locations
                model.mesh_alpha = 0.4; % transparency of the links
                for i=1:model.Nmeshes
                    model.mesh_color{i} = [226 209 24]/255; % color of the body links
                end
            end

            % load stl files into patch structure
            for i=1:model.Nmeshes
                try
                    model.fv{i} = stlread(strcat(model.mesh_name{i},'.stl'));

                    % reduce number of faces in patch object
                    model.fv{i} = reducepatch(model.fv{i}, model.reduce_patch_ratio);
                    fv_tmp = model.fv{i};
                    fv_tmp.vertices = model.scale_ratio*fv_tmp.vertices;

                    mesh_color = model.mesh_color{i};
                    mesh_alpha = model.mesh_alpha;

                    model.h_fv{i} = patch(fv_tmp,'FaceColor', mesh_color, ...
                        'EdgeColor', 'none',        ...
                        'FaceLighting', 'gouraud',     ...
                        'AmbientStrength', 0.15, 'FaceAlpha', mesh_alpha);

                    % mesh scaling
                    model.h_fv{i}.Vertices = model.h_fv{i}.Vertices;

                    % get vertices data
                    v_min = min(fv_tmp.vertices);
                    v_max = max(fv_tmp.vertices);
                    v_len = v_max - v_min;

                    model.vdata.vmins{i} = v_min;
                    model.vdata.vmaxs{i} = v_max;
                    model.vdata.vlens{i} = v_len;
                end
            end

            if(model.display_inertia == 1)
                model.fv_com = cell(size(model.N_ext,2),1);
                model.h_fv_com = cell(size(model.N_ext,2),1);
                model.p_fv_com = cell(size(model.N_ext,2),1);

                for i=1:model.N_ext
                    %                     h = drawSphere([0 0 0 0.1]);
                    [X, Y, Z] = drawSphere([0 0 0 0.02]);
                    C = ones(size(X));
                    model.h_fv_com{i} = patch(X,Y,Z,C, 'FaceColor', [0 0 0.8], ...
                        'EdgeColor', 'none', 'FaceAlpha', 1);
                end
            end

            %= apply the transformation
            % get the position and orientation from the model definition
            for i=1:model.N_ext

                if(i <= model.Nmeshes)
                    T = model.T0{i};

                    model.T_fv{i} = T; % store the current T
                    t = model.transl(T);
                    R = model.t2r(T);
                    rpy_offset = model.R_mesh_offset{i};
                    R_offset = model.rpy2R(rpy_offset(1),rpy_offset(2),rpy_offset(3));

                    V = get(model.h_fv{i},'Vertices');
                    V = V*(R*R_offset)' + ones(size(V,1),1)*(t + R*model.p_mesh_offset{i})';
                    set(model.h_fv{i},'Vertices',V);
                    model.fv{i}.vertices = V;

                    if(model.display_inertia >= 1 && i<=model.NB)
                        t = T*[model.i_p_com{i} ;1]; % store the current T
                        model.p_fv_com{i} = t(1:3);
                        V = get(model.h_fv_com{i},'Vertices');
                        model.fv_com{i} = V;
                        V = V + ones(size(V,1),1)*t(1:3)';
                        set(model.h_fv_com{i},'Vertices',V);
                    end

                else % for end-effector frame
                    if(model.display_inertia == 1)
                        t = T*[model.i_p_com{i};1]; % store the current T
                        model.p_fv_com{i} = t(1:3);
                        V = get(model.h_fv_com{i},'Vertices');
                        V = V + ones(size(V,1),1)*t(1:3)';
                        set(model.h_fv_com{i},'Vertices',V);
                    end
                end

                %- display frames and get handles
                if(model.display_frames == 1)
                    model.h_frame{i} = model.drawFrameHandle(T, model.frame_len , model.frame_width);
                end
            end

            % initialize collision shapes
            model.initCollisionShapes();

            %= set display settings
            if (isModelGhost == true)
%                 model.display_frames = 0; % turn on/off display coordinate frames
%                 model.display_inertia = 0; % turn on/off display CoM locations
%                 model.mesh_alpha = 0.4; % transparency of the links
%                 for i=1:model.Nmeshes
%                     model.mesh_color{i} = [226 209 24]/255; % color of the body links
%                 end

%                 material('dull');
            else
                % Add a camera light, and tone down the specular highlighting
                camlight('headlight');

                material('dull');

                % Fix the axes scaling, and set a nice view angle
                axis('image');
            end

            % set the animation flag
            model.display_animation = true;
%             grid on;
        end

        function initCollisionShapes(model, displayFlag)
            % INITCOLLISIONSHAPES initialize collision shapes for visualization
            %   INITCOLLISIONSHAPES(model, displayFlag)
            %
            %   This function initializes collision shapes for the current
            %   robot model.
            %
            %   See also initAnimation.

            if ~exist('displayFlag','var')
                displayFlag = true;
            end

            % if collision shape prestored file exists, set the flag;
            if (exist(model.collisionShapeFile,'file'))
                model.load_stored_collision_shapes = 1;
            else
                model.load_stored_collision_shapes = 0;
            end

            for i=1: model.Nmeshes
                % get transforms for each link
                T = model.T0{i};
                t = model.transl(T);
                R = model.t2r(T);

                % display bounding box
                if(model.display_collision_shapes)
                    model.T_fv_col{i} = T;

                    if (model.display_collision_shapes == 1) % cube
                        % get data for vertices
                        v_min = model.vdata.vmins{i};
                        v_len = model.vdata.vlens{i};

                        % plot cube
                        opt = struct('color',model.collision_shape_color,'alpha',model.collision_shape_alpha,'sub_idx',i,...
                            'line_style','-','line_color',model.collision_shape_color);
                        model.h_fv_col{i} = plot_cube(t',v_min,v_len,R',opt);
                        [V_col, F_col] = get_cube_FV(t',v_min,v_len,R');
                        model.fv_col{i}.Faces = F_col;
                        model.fv_col{i}.Vertices = V_col;

                    elseif (model.display_collision_shapes == 2) % capsule
                        reduce_patch_ratio = 0.1;
                        fv_tmp = reducepatch(model.fv{i}, reduce_patch_ratio);

                        if(model.load_stored_collision_shapes)
                            % load stored capsule collision shapes
                            load(model.collisionShapeFile,'caps_param');
                        else
                            %= find a bounding capsule and draw it
                            % get bounding capsule parameters
                            caps_param{i} = findBoundingCapsule(fv_tmp.vertices);

                            % store transforms w.r.t link frame
                            caps_param{i}.t_link = R'*(caps_param{i}.t - t);
                            caps_param{i}.e1_link = R'*(caps_param{i}.param.e1 - t);
                            caps_param{i}.e2_link = R'*(caps_param{i}.param.e2 - t);
                        end

                        % transform the local transforms
                        caps_param{i}.t = R*caps_param{i}.t_link + t;
                        caps_param{i}.param.e1 = R*caps_param{i}.e1_link + t;
                        caps_param{i}.param.e2 = R*caps_param{i}.e2_link + t;

                        % get rotation matrix
                        caps_param{i}.R = vecRotMat([0 0 1], SpatialRobotModel.unit(caps_param{i}.param.e2-caps_param{i}.param.e1)');

                        % draw the capsule if requested
                        if (displayFlag)
                            model.h_fv_col{i} = drawCapsuleSimple(caps_param{i}.R, caps_param{i}.t, caps_param{i}.param, ...
                                'FaceColor', model.collision_shape_color, 'FaceAlpha', model.collision_shape_alpha, ...
                                'EdgeColor', model.collision_shape_color, 'EdgeAlpha', 0);
                        end

                        % set the type
                        model.h_fv_col{i}.type = 1; % robot links

                        % load local transforms
                        model.h_fv_col{i}.t_link = caps_param{i}.t_link;
                        model.h_fv_col{i}.e1_link = caps_param{i}.e1_link;
                        model.h_fv_col{i}.e2_link = caps_param{i}.e2_link;
                    end
                end
            end

            % store collision shape parameters
            if (model.display_collision_shapes == 2 && ~model.load_stored_collision_shapes)
                save(model.collisionShapeFile, 'caps_param');
            end
        end

        function updateAnimation(model)
            % UPDATEANIMATION update the graphical robot model.
            %   updateAnimation(model)
            %
            %   This function updates the tranforms of the graphics handles
            %   in the animation model.
            %
            %   See also initAnimation.

            % only update link mesh files
            Nmeshes = model.Nmeshes;

            %= apply the transformation
            % get the position and orientation from the model definition
            for i=1:model.N_ext
                T_new = model.T0{i};

                if(i <= Nmeshes)
                    T_old = model.T_fv{i};
                    t_new = model.transl(T_new);
                    R_new = model.t2r(T_new);
                    t_old = model.transl(T_old);
                    R_old = model.t2r(T_old);

                    % determine local transformation
                    R = R_new*R_old';
                    t = t_new - R*t_old;

                    V = get(model.h_fv{i},'Vertices');
                    V = V*R' + ones(size(V,1),1)*t';
                    set(model.h_fv{i},'Vertices',V);

                    model.T_fv{i} = T_new;

                    if(model.display_inertia == 1)
                        t_old_com = model.p_fv_com{i};
                        t_new_com = model.T0{i}*[model.i_p_com{i} ;1]; % store the current T
                        t_com = t_new_com(1:3) - t_old_com;

                        V = get(model.h_fv_com{i},'Vertices');
                        V = V + ones(size(V,1),1)*t_com(1:3)';
                        set(model.h_fv_com{i},'Vertices',V);

                        model.p_fv_com{i} = t_new_com(1:3);
                    end

                    % update verticies in the collision shapes
                    if(model.display_collision_shapes)
                        %- get relative transform to last updated transform

                        % get last updated transform
                        T_old2 = model.T_fv_col{i};
                        t_old2 = model.transl(T_old2);
                        R_old2 = model.t2r(T_old2);

                        % determine relative transformation
                        R_col = R_new*R_old2';
                        t_col = t_new - R_col*t_old2;

                        if(model.display_collision_shapes == 1)
                            V_col = get(model.h_fv_col{i},'Vertices');
                            V_col = V_col*R_col' + ones(size(V_col,1),1)*t_col';
                            set(model.h_fv_col{i},'Vertices',V_col);

                        elseif(model.display_collision_shapes == 2) % capsule
                            V_col = get(model.h_fv_col{i}.bodies(1),'Vertices');
                            V_col = V_col*R_col' + ones(size(V_col,1),1)*t_col';
                            set(model.h_fv_col{i}.bodies(1),'Vertices',V_col);

                            % transform center line segments of capsules
                            model.h_fv_col{i}.t = R_new*model.h_fv_col{i}.t + t_new;
                            model.h_fv_col{i}.param.e1 = R_new*model.h_fv_col{i}.param.e1 + t_new;
                            model.h_fv_col{i}.param.e2 = R_new*model.h_fv_col{i}.param.e2 + t_new;
                        end

                        model.fv_col{i}.Vertices = V_col;
                        model.T_fv_col{i} = T_new;
                    end
                else

                    if(model.display_inertia == 1)
                        t_old_com = model.p_fv_com{i};
                        t_new_com = model.T0{i}*[model.i_p_com{i}; 1]; % store the current T
                        t_com = t_new_com(1:3) - t_old_com;

                        V = get(model.h_fv_com{i},'Vertices');
                        V = V + ones(size(V,1),1)*t_com(1:3)';
                        set(model.h_fv_com{i},'Vertices',V);

                        model.p_fv_com{i} = t_new_com(1:3);
                    end
                end

                %- display frames and get handles
                if(model.display_frames == 1)
                    model.drawFrameHandleSet(model.h_frame{i}, T_new, model.frame_len);
                end
            end

        end

        % update capsule transforms without visualization
        function updateCollisionCapsules(model)
            plotForDebugging = 0;
            for i=1:model.Nmeshes
                if (model.display_collision_shapes == 2) % capsule
                    T_i = model.T0{i};
                    t_i = model.transl(T_i);
                    R_i = model.t2r(T_i);

                    model.h_fv_col{i}.t = R_i*model.h_fv_col{i}.t_link + t_i;
                    model.h_fv_col{i}.param.e1 = R_i*model.h_fv_col{i}.e1_link + t_i;
                    model.h_fv_col{i}.param.e2 = R_i*model.h_fv_col{i}.e2_link + t_i;

                    if (plotForDebugging)
                        hold on;
                        plot3(model.h_fv_col{i}.t(1), model.h_fv_col{i}.t(2), model.h_fv_col{i}.t(3), 'r.', 'MarkerSize', 10);
                        plot3(model.h_fv_col{i}.param.e1(1), model.h_fv_col{i}.param.e1(2), model.h_fv_col{i}.param.e1(3), 'b.', 'MarkerSize', 10);
                        plot3(model.h_fv_col{i}.param.e2(1), model.h_fv_col{i}.param.e2(2), model.h_fv_col{i}.param.e2(3), 'b.', 'MarkerSize', 10);
                    end
                end
            end
        end

        % toggle display collision shapes
        function toggleDisplayCollisionShapes(model)
            if (model.display_collision_shapes == 0)

                % use the capsule shape by default
                model.display_collision_shapes = 2;

                % check if the patches exist. If not, initialize them
                if(~isempty(model.h_fv_col))

                    % turn on the collision shapes visiblity settings
                    for i=1:model.Nmeshes
                        if(model.display_collision_shapes == 1)
                            set(model.h_fv_col{i}, 'Visible', 'on');
                        elseif(model.display_collision_shapes == 2)
                            set(model.h_fv_col{i}.bodies, 'Visible', 'on');
                        end
                    end
                    % if it exists, update the collision shapes
                    model.updateAnimation();
                else
                    model.initCollisionShapes();
                end

            elseif (model.display_collision_shapes > 0) % if the collision shape exists

                % set the current collision shapes' visible settings off
                model.hideCollisionShapes();

                % set the flag to disable collision shape update
                model.display_collision_shapes = 0;

            end
        end

        function hideCollisionShapes(model)
            % HIDECOLLISIONSHAPES hide the collision shapes
            %   HIDECOLLISIONSHAPES(model)
            %
            %   This function sets the visible settings off for the
            %   collision shape handles in the animation model.

            if(~isempty(model.h_fv_col))
                for i=1:model.Nmeshes
                    if(model.display_collision_shapes == 1)
                        set(model.h_fv_col{i}, 'Visible', 'off');
                    elseif(model.display_collision_shapes == 2)
                        set(model.h_fv_col{i}.bodies, 'Visible', 'off');
                    end
                end
            end
        end

        function hideAnimation(model)
            % HIDEANIMATION hide the robot graphical model
            %   HIDEANIMATION(model)
            %
            %   This function sets the visible settings off for the
            %   graphics handles in the animation model.

            if(~isempty(model.h_fv))
                if (model.display_animation)
                    for i=1:model.Nmeshes
                        set(model.h_fv{i}, 'Visible', 'off');

                        if(model.display_collision_shapes == 1)
                            set(model.h_fv_col{i}, 'Visible', 'off');
                        elseif(model.display_collision_shapes == 2)
                            set(model.h_fv_col{i}.bodies(1), 'Visible', 'off');
                        end
                    end

                    if(model.display_inertia == 1)
                        for i=1:model.NB
                            set(model.h_fv_com{i}, 'Visible', 'off');
                        end
                    end

                    %- display frames and get handles
                    for i=1:model.N_ext
                        if(model.display_frames == 1)
                            set(model.h_frame{i}, 'Visible', 'off');
                        end
                    end
                    % turn off the flag
                    model.display_animation = false;
                else
                    for i=1:model.Nmeshes
                        set(model.h_fv{i}, 'Visible', 'on');

                        if(model.display_collision_shapes == 1)
                            set(model.h_fv_col{i}, 'Visible', 'on');
                        elseif(model.display_collision_shapes == 2)
                            set(model.h_fv_col{i}.bodies(1), 'Visible', 'on');
                        end
                    end

                    if(model.display_inertia == 1)
                        for i=1:model.NB
                            set(model.h_fv_com{i}, 'Visible', 'on');
                        end
                    end

                    %- display frames and get handles
                    for i=1:model.N_ext
                        if(model.display_frames == 1)
                            set(model.h_frame{i}, 'Visible', 'on');
                        end
                    end
                    % turn on the flag
                    model.display_animation = true;
                end
            end
        end

        %% other useful functions

        function isLimitViolated = checkJointLimits(model, q)
            % CHECKJOINTLIMITS check joint limits
            %   isLimitViolated = checkJointLimits(model, q)
            %
            %   This function returns true if joint limits are violated

            isLimitViolated = false;

            iJoint = 1;
            while (isLimitViolated == false && iJoint <= model.NJ)
                if(~isempty(model.joints{iJoint}))
                    q_min = model.joint_pos_min{iJoint+1};
                    q_max = model.joint_pos_max{iJoint+1};

                    if(q(iJoint) >= q_max)
                        isLimitViolated = true;
                    elseif(q(iJoint) <= q_min)
                        isLimitViolated = true;
                    end
                end

                iJoint = iJoint+1;
            end
        end

        function clampJointsToLimits(model)
            % CLAMPJOINTSTOLIMITS clamp joint angle to the limit
            %   CLAMPJOINTSTOLIMITS(model)
            %
            %   This function clamps joint positions to its limits, and
            %   updates actuated joint position for the model.
            %
            %   See also checkJointLimits.

            for i=1:model.NJ
                q_min = model.joint_pos_min{i+1};
                q_max = model.joint_pos_max{i+1};

                if(model.q_a(i) > q_max)
                    model.q_a(i) = q_max;
                elseif(model.q_a(i) < q_min)
                    model.q_a(i) = q_min;
                end
            end
        end

        function updatePayload(model, dist, weight, com, TcpInd)
            % UPDATEPAYLOAD update endeffector information (payload and offset)
            %   UPDATEPAYLOAD(model, dist, weight, com)
            %
            %   This function updates end-effector information with payload
            %   mass (in kg) value and offset (dist in meter) and
            %   Center-of-Mass (com in meter) vectors (1x3).

            if ~exist('TcpInd','var')
                TcpInd = model.ind_H;
            end

            % transform the specified dist and com
            dist_pi = model.pi_R_i{TcpInd}*dist';

            model.pi_p_i{TcpInd} = model.param_backup.pi_p_i{TcpInd} + dist_pi;
            model.m{TcpInd} = model.param_backup.m{TcpInd} + weight;
            model.i_p_com{TcpInd} = com';

            % print out the updated values
            fprintf('\n\n End-effector has been updated.\n -- Distance: (%.3f, %.3f, %.3f)\n -- Mass: (%.3f)\n -- CoM: (%.3f, %.3f, %.3f)\n', ...
                dist(1), dist(2), dist(3), model.m{model.ind_H}, com(1), com(2), com(3));

            % copy model parameters into a struct variable
            model.model_param = model.model2struct();
        end

        % normalize the joint positions (joint positions are mapped to [0,1] range)
        function q_out = p2u(model, q_in)
            q_out = zeros(model.NJ,1);
            for p_idx = 1:model.NJ
                q_out(p_idx) = (q_in(p_idx)-model.joint_pos_min{p_idx+1}) / ...
                    (model.joint_pos_max{p_idx+1}-model.joint_pos_min{p_idx+1});
            end
        end

        % map normalized positions back to actual joint positions
        function q_out = u2p(model, q_in)
            q_out = zeros(model.NJ,1);
            for p_idx = 1:model.NJ
                q_out(p_idx) = model.joint_pos_min{p_idx+1}+...
                    (model.joint_pos_max{p_idx+1}-model.joint_pos_min{p_idx+1})*q_in(p_idx);
            end
        end

    end


end