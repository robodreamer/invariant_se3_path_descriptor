function handle_plot = plot_se3_trajectories(trajectories, colors, titletext, params)

% Parameters for plotting (same as before)
inc = 5;
lineLength = 0.1;
lineWidth = 1.5;

% Define start and end marker colors and styles
startMarkerColor = [0, 1, 0]; % Green for start
endMarkerColor = [1, 0, 0]; % Red for end
startMarkerStyle = 'o'; % Circle marker
endMarkerStyle = 'x'; % Cross marker

% Initialize axis limits
minX = Inf; maxX = -Inf;
minY = Inf; maxY = -Inf;
minZ = Inf; maxZ = -Inf;

% Calculate scale based on all trajectories
if (params.auto_calculate_scale)
    scale = 1;
    for i = 1:length(trajectories)
        p_obj = trajectories{i}.pos_data;
        scale = max(scale, (max(p_obj(:,1))-min(p_obj(:,1))) + (max(p_obj(:,2))-min(p_obj(:,2))) + (max(p_obj(:,3))-min(p_obj(:,3))));
    end
else
    scale = params.scale;
end

lencx = 0.1*scale;
lency = 0.060*scale;
lencz = 0.060*scale;
len = 1.2*lencx;

% Setup figure
hFig = figure; clf;
set(gcf,'units','normalized','outerposition',[0 0.2 0.5 0.7]);
hold on; axis equal; view(-161,25); grid on; box on;
xlabel('$x$[m]','Interpreter','LaTex','FontSize',18)
ylabel('$y$[m]','Interpreter','LaTex','FontSize',18)
zlabel('$z$[m]','Interpreter','LaTex','FontSize',18)
title(titletext)

% Loop through each trajectory
num_traj = length(trajectories);
handle_plot = zeros(1, num_traj);
for idx = 1:num_traj

    trajectory = trajectories{idx};
    p_obj = trajectory.pos_data;
    R = trajectory.rot_data;
    color = colors(idx,:);

    [N,M] = size(p_obj);

    % Update axis limits
    minX = min(minX, min(p_obj(:,1)));
    maxX = max(maxX, max(p_obj(:,1)));
    minY = min(minY, min(p_obj(:,2)));
    maxY = max(maxY, max(p_obj(:,2)));
    minZ = min(minZ, min(p_obj(:,3)));
    maxZ = max(maxZ, max(p_obj(:,3)));

    % Plot trajectory
    handle_plot(idx) = plot3(p_obj(:,1), p_obj(:,2), p_obj(:,3), 'Color', color, 'linewidth', 2);
    % Mark start and end points
    plot3(p_obj(1,1), p_obj(1,2), p_obj(1,3), ...
          startMarkerStyle, 'Color', startMarkerColor, 'MarkerSize', 10, 'LineWidth', 2);
    plot3(p_obj(end,1), p_obj(end,2), p_obj(end,3), ...
          endMarkerStyle, 'Color', endMarkerColor, 'MarkerSize', 10, 'LineWidth', 2);
end

% Loop through each trajectory
if params.show_rotation
    for idx = 1:length(trajectories)

        trajectory = trajectories{idx};
        p_obj = trajectory.pos_data;
        R = trajectory.rot_data;
        color = colors(idx,:);

        [N,M] = size(p_obj);

        % Draw orientation and cubes

        Rx = []; Ry = []; Rz = []; p = [];
        for j = round(linspace(1, N, inc))
            Rx = [Rx; R(1:3,1,j)'];
            Ry = [Ry; R(1:3,2,j)'];
            Rz = [Rz; R(1:3,3,j)'];
            p = [p; p_obj(j,:)];

            if params.show_coordinates
                T = [R(:,:,j), p_obj(j,:)'; 0 0 0 1];
                SpatialRobotModel.drawFrame(T, lineLength, lineWidth);
            end
        end

        % Draw cubes
        for j = round(linspace(1, N, inc))
            fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
            vm = repmat(p_obj(j,:), 8, 1) + [ zeros(3,1) ...
                R(1:3,1:3,j)*[lencx; 0; 0] R(1:3,1:3,j)*[lencx; lency; 0] ...
                R(1:3,1:3,j)*[0; lency; 0] R(1:3,1:3,j)*[0; 0; lencz] ...
                R(1:3,1:3,j)*[lencx; 0; lencz] R(1:3,1:3,j)*[lencx; lency; lencz] ...
                R(1:3,1:3,j)*[0; lency; lencz] ]';
            patch('Vertices', vm, 'Faces', fm, 'EdgeAlpha', 0.8, ...
                'FaceColor', color, ...
                'FaceAlpha', 0.3, 'EdgeColor', [0.15 0.15 0.90], 'LineWidth', 1.0);
        end
    end
end

% Set axis limits
margin = 0.2;
axis([minX - margin, maxX + margin, ...
    minY - margin, maxY + margin, ...
    minZ - margin, maxZ + margin]);

% Add legend if legend texts are provided
if exist('legend_texts', 'var') && ~isempty(legend_texts)
    legend(handle_plot, legend_texts, 'Location', 'best');
end

dragzoom(hFig);
end
