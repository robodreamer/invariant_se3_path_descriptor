function plot_trajectory(trajectory1,trajectory2,titletext,show_rotation)

p_obj_demo = trajectory1.pos_data;
R_demo = trajectory1.rot_data;
if ~isempty(trajectory2)
    p_obj_recon = trajectory2.pos_data;
    R_recon = trajectory2.rot_data;
end

% Parameters for plotting
inc = 5; % how many rigid bodies drawn along trajectory
% scale calculation to get an idea of how big figure is
scale1 = (max(p_obj_demo(:,1))-min(p_obj_demo(:,1))) + (max(p_obj_demo(:,2))-min(p_obj_demo(:,2))) + (max(p_obj_demo(:,3))-min(p_obj_demo(:,3)));
if ~isempty(trajectory2)
    scale2 = (max(p_obj_recon(:,1))-min(p_obj_recon(:,1))) + (max(p_obj_recon(:,2))-min(p_obj_recon(:,2))) + (max(p_obj_recon(:,3))-min(p_obj_recon(:,3)));
else
    scale2 = 0;
end
scale = max(scale1,max(scale2,1));

lencx = 0.1*scale; % length of cube
lency = 0.060*scale; % length of cube
lencz = 0.060*scale; % length of cube
len = 1.2*lencx; % length arrow
t = 1; % thickness arrowhead
linewidth = '1'; % width arrow

figure; clf;
set(gcf,'units','normalized','outerposition',[0 0.2 0.5 0.7]);
hold on; axis equal; view(-161,25); grid on; box on;
xlabel('$x$[m]','Interpreter','LaTex','FontSize',18)
ylabel('$y$[m]','Interpreter','LaTex','FontSize',18)
zlabel('$z$[m]','Interpreter','LaTex','FontSize',18)
title(titletext)

%% Trajectory 1 (demonstration)
[N1,M1] = size(p_obj_demo);

% plot trajectory reference point
for j=1:M1/3 % for multiple markers
    plot3(p_obj_demo(:,3*j-2),p_obj_demo(:,3*j-1),p_obj_demo(:,3*j),'b.-','linewidth',1.5);
    plot3(p_obj_demo(1,3*j-2),p_obj_demo(1,3*j-1),p_obj_demo(1,3*j),'bo','MarkerSize',7,'LineWidth',2);
    plot3(p_obj_demo(N1,3*j-2),p_obj_demo(N1,3*j-1),p_obj_demo(N1,3*j),'b*','MarkerSize',7,'LineWidth',2);
end

% Draw orientation rigid body demonstration
if show_rotation

    % Draw arrows
    Rx_demo = []; Ry_demo = []; Rz_demo = []; p_demo = [];
    for j=round(linspace(1,N1,inc))
        Rx_demo = [Rx_demo ; R_demo(1:3,1,j)'];
        Ry_demo = [Ry_demo ; R_demo(1:3,2,j)'];
        Rz_demo = [Rz_demo ; R_demo(1:3,3,j)'];
        p_demo = [p_demo ; p_obj_demo(j,:)];
    end
    arrow3(p_demo,p_demo+len*Rx_demo,['_r' linewidth],t,2*t)
    arrow3(p_demo,p_demo+len*Ry_demo,['_e' linewidth],t,2*t)
    arrow3(p_demo,p_demo+len*Rz_demo,['_b' linewidth],t,2*t)
    axis equal;

    % Draw cubes
    for j=round(linspace(1,N1,inc))
        fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
        vm = repmat(p_obj_demo(j,:),8,1) + [ zeros(3,1) ...
            R_demo(1:3,1:3,j)*[lencx; 0; 0] R_demo(1:3,1:3,j)*[lencx; lency; 0] ...
            R_demo(1:3,1:3,j)*[0; lency; 0] R_demo(1:3,1:3,j)*[0; 0; lencz] ...
            R_demo(1:3,1:3,j)*[lencx; 0; lencz] R_demo(1:3,1:3,j)*[lencx; lency; lencz] ...
            R_demo(1:3,1:3,j)*[0; lency; lencz] ]';
        patch('Vertices',vm,'Faces',fm, 'EdgeAlpha',0.8,'FaceColor',[0.25 0.80 0.80],'FaceAlpha',0.30,'EdgeColor',[0.15 0.15 0.90],'LineWidth',1.00);
    end
end

%% Trajectory 2 (generated)
if ~isempty(trajectory2)

    [N2,M2] = size(p_obj_recon);

    % plot trajectory reference point
    for j=1:M2/3 % for multiple markers
        plot3(p_obj_recon(:,3*j-2),p_obj_recon(:,3*j-1),p_obj_recon(:,3*j),'r.-','linewidth',1.5);
        plot3(p_obj_recon(1,3*j-2),p_obj_recon(1,3*j-1),p_obj_recon(1,3*j),'ro','MarkerSize',7,'LineWidth',2);
        plot3(p_obj_recon(N2,3*j-2),p_obj_recon(N2,3*j-1),p_obj_recon(N2,3*j),'r*','MarkerSize',7,'LineWidth',2);
    end

    % Draw orientation rigid body reconstruction
    if show_rotation
        Rx_recon = []; Ry_recon = []; Rz_recon = []; p_recon = [];
        Rx_FS_recon = []; Ry_FS_recon = []; Rz_FS_recon = []; p_FS_recon = [];
        for j=round(linspace(1,length(p_obj_recon),inc))
            Rx_recon = [Rx_recon ; R_recon(1:3,1,j)'];
            Ry_recon = [Ry_recon ; R_recon(1:3,2,j)'];
            Rz_recon = [Rz_recon ; R_recon(1:3,3,j)'];
            p_recon = [p_recon ; p_obj_recon(j,:)];
        end
        arrow3(p_recon,p_recon+len*Rx_recon,['_r' linewidth],t,2*t)
        arrow3(p_recon,p_recon+len*Ry_recon,['_e' linewidth],t,2*t)
        arrow3(p_recon,p_recon+len*Rz_recon,['_b' linewidth],t,2*t)
        axis equal
        % Draw cube reconstruction
        for j=round(linspace(1,N2,inc))
            fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
            vm = repmat(p_obj_recon(j,:),8,1) + [ zeros(3,1) ...
                R_recon(1:3,1:3,j)*[lencx; 0; 0] R_recon(1:3,1:3,j)*[lencx; lency; 0] ...
                R_recon(1:3,1:3,j)*[0; lency; 0] R_recon(1:3,1:3,j)*[0; 0; lencz] ...
                R_recon(1:3,1:3,j)*[lencx; 0; lencz] R_recon(1:3,1:3,j)*[lencx; lency; lencz] ...
                R_recon(1:3,1:3,j)*[0; lency; lencz] ]';
            patch('Vertices',vm,'Faces',fm, 'EdgeAlpha',0.8,'FaceColor',[0.80 0.25 0.25],'FaceAlpha',0.10,'EdgeColor',[0.90 0.15 0.15],'LineWidth',1.00);
        end
    end
end
%zoom(0.75)
