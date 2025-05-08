% % MATLAB code for generating random irregular shapes with central points
% close all
% clc
% clear
% % Figure setup
% % figure('Position', [100, 100, 1000, 800]);
% % hold on;
% % axis equal;
% % grid on;
% % title('Random Irregular Shapes with Central Points');
% % xlabel('X');
% % ylabel('Y');
% % xlim([-15 15]);
% % ylim([-15 15]);
% % 
% % % Set random seed for reproducibility (comment out for truly random shapes)
% % rng(123);
% % 
% % % Color palette using standard MATLAB colors
% % colors = {'b', 'r', 'g', 'm', 'c', 'y', [0.8 0.2 0.2], [0.2 0.6 0.2], [0.6 0.2 0.8], [0.8 0.8 0.2]};
% % 
% % % 1. Perlin Noise-based Shape
% % shape_center = [5, 5];
% % t = linspace(0, 2*pi, 100);
% % r = 2 + 0.8*sin(5*t) + 0.5*cos(7*t) + 0.3*sin(11*t) + 0.4*cos(13*t);
% % x = r.*cos(t) + shape_center(1);
% % y = r.*sin(t) + shape_center(2);
% % fill(x, y, colors{1}, 'FaceAlpha', 0.7);
% % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % text(shape_center(1), shape_center(2)-0.5, 'Perlin Shape', 'HorizontalAlignment', 'center');
% % 
% % % 2. Random Polygon with many vertices
% % num_vertices = randi([8, 15]);
% % shape_center = [-7, 6];
% % r = 1.5 + 1.5*rand(1, num_vertices);
% % t = linspace(0, 2*pi, num_vertices+1);
% % t = t(1:end-1);
% % x = r.*cos(t) + shape_center(1);
% % y = r.*sin(t) + shape_center(2);
% % fill(x, y, colors{2}, 'FaceAlpha', 0.7);
% % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % text(shape_center(1), shape_center(2)-0.5, 'Random Polygon', 'HorizontalAlignment', 'center');
% % 
% % % % 3. Gaussian Random Walk Boundary
% % % shape_center = [-6, -3];
% % % num_points = 150;
% % % t = linspace(0, 2*pi, num_points);
% % % % Base radius
% % % r_base = 2.5*ones(size(t));
% % % % Add random noise
% % % noise_scale = 1.2;
% % % r_noise = cumsum(noise_scale*randn(1, num_points));
% % % % Normalize the noise to have zero mean
% % % r_noise = r_noise - mean(r_noise);
% % % % Apply the noise
% % % r = r_base + r_noise;
% % % x = r.*cos(t) + shape_center(1);
% % % y = r.*sin(t) + shape_center(2);
% % % fill(x, y, colors{3}, 'FaceAlpha', 0.7);
% % % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % % text(shape_center(1), shape_center(2)-0.5, 'Random Walk', 'HorizontalAlignment', 'center');
% % 
% % % 4. Voronoi-inspired random shape
% % shape_center = [3, -4];
% % % Generate random points around center
% % num_points = 20;
% % rand_points_x = shape_center(1) + 4*(rand(1, num_points) - 0.5);
% % rand_points_y = shape_center(2) + 4*(rand(1, num_points) - 0.5);
% % % Select subset of points to form a convex hull
% % num_hull_points = randi([6, 10]);
% % hull_indices = randperm(num_points, num_hull_points);
% % hull_x = rand_points_x(hull_indices);
% % hull_y = rand_points_y(hull_indices);
% % % Sort points in counterclockwise order around their mean
% % mean_x = mean(hull_x);
% % mean_y = mean(hull_y);
% % angles = atan2(hull_y - mean_y, hull_x - mean_x);
% % [~, sorted_indices] = sort(angles);
% % hull_x = hull_x(sorted_indices);
% % hull_y = hull_y(sorted_indices);
% % % Close the shape
% % hull_x = [hull_x, hull_x(1)];
% % hull_y = [hull_y, hull_y(1)];
% % % Perturb the hull points slightly to make it more irregular
% % for i = 1:length(hull_x)-1
% %     hull_x(i) = hull_x(i) + 0.5*randn();
% %     hull_y(i) = hull_y(i) + 0.5*randn();
% % end
% % fill(hull_x, hull_y, colors{4}, 'FaceAlpha', 0.7);
% % % Plot all original random points
% % plot(rand_points_x, rand_points_y, 'k.', 'MarkerSize', 8);
% % % Mark the center
% % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % text(shape_center(1), shape_center(2)-0.5, 'Voronoi Shape', 'HorizontalAlignment', 'center');
% % 
% % % 5. Chaotic Attractor Inspired Shape
% % shape_center = [6, -7];
% % % Generate a strange attractor-like shape
% % a = 0.9; b = 0.6; c = 0.2;
% % num_points = 2000;
% % x_att = zeros(1, num_points);
% % y_att = zeros(1, num_points);
% % x_att(1) = 0.1; y_att(1) = 0.1;
% % for i = 2:num_points
% %     x_att(i) = 1 - a*x_att(i-1)^2 + b*y_att(i-1);
% %     y_att(i) = c*x_att(i-1);
% % end
% % % Scale and translate
% % x_att = x_att*1.5 + shape_center(1);
% % y_att = y_att*1.5 + shape_center(2);
% % % Create a boundary shape from the attractor
% % % Check if boundary function exists (newer MATLAB versions)
% % if exist('boundary', 'file')
% %     k = boundary(x_att', y_att', 0.5); % 0.5 is shrink factor - adjust if needed
% %     fill(x_att(k), y_att(k), colors{5}, 'FaceAlpha', 0.7);
% % else
% %     % Fallback if boundary function doesn't exist
% %     % Just plot the points
% %     scatter(x_att, y_att, 10, colors{5}, 'filled', 'MarkerFaceAlpha', 0.4);
% % end
% % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % text(shape_center(1), shape_center(2)-0.5, 'Chaotic Shape', 'HorizontalAlignment', 'center');
% % 
% % % 6. Blob with random bumps
% % shape_center = [-3, 2];
% % t = linspace(0, 2*pi, 200);
% % % Generate a series of random bumps at different frequencies
% % bump_count = randi([5, 12]);
% % bump_freqs = randi([2, 8], 1, bump_count);
% % bump_amps = 0.2 + 0.4*rand(1, bump_count);
% % bump_phases = 2*pi*rand(1, bump_count);
% % r = ones(size(t))*2;
% % for i = 1:bump_count
% %     r = r + bump_amps(i)*sin(bump_freqs(i)*t + bump_phases(i));
% % end
% % x = r.*cos(t) + shape_center(1);
% % y = r.*sin(t) + shape_center(2);
% % fill(x, y, colors{6}, 'FaceAlpha', 0.7);
% % plot(shape_center(1), shape_center(2), 'k*', 'MarkerSize', 10);
% % text(shape_center(1), shape_center(2)-0.5, 'Random Blob', 'HorizontalAlignment', 'center');
% % 
% % % 7. Amoeba-like shape with pseudopods
% % amoeba_center = [10, -2];
% % % Base circle
% % t = linspace(0, 2*pi, 200);
% % r_base = 1.5*ones(size(t));
% % % Add pseudopods
% % num_pods = randi([3, 7]);
% % pod_centers = linspace(0, 2*pi, num_pods+1);
% % pod_centers = pod_centers(1:end-1);
% % pod_widths = 0.3 + 0.5*rand(1, num_pods);
% % pod_lengths = 1 + 1.5*rand(1, num_pods);
% % r = r_base;
% % for i = 1:num_pods
% %     pod_factor = pod_lengths(i) * exp(-((t - pod_centers(i)).^2) / (2*pod_widths(i)^2));
% %     r = r + pod_factor;
% % end
% % x = r.*cos(t) + amoeba_center(1);
% % y = r.*sin(t) + amoeba_center(2);
% % fill(x, y, colors{7}, 'FaceAlpha', 0.7);
% % plot(amoeba_center(1), amoeba_center(2), 'k*', 'MarkerSize', 10);
% % text(amoeba_center(1), amoeba_center(2)-0.5, 'Amoeba Shape', 'HorizontalAlignment', 'center');
% 
% % % 8. Splatter shape with irregular protrusions
% % splatter_center = [0, 7];
% % num_splats = randi([15, 25]);
% % splat_angles = 2*pi*rand(1, num_splats);
% % splat_distances = 0.2 + 2.5*rand(1, num_splats);
% % splat_sizes = 0.2 + 0.8*rand(1, num_splats);
% % 
% % % Create a series of overlapping circles
% % for i = 1:num_splats
% %     circle_x = splatter_center(1) + splat_distances(i)*cos(splat_angles(i));
% %     circle_y = splatter_center(2) + splat_distances(i)*sin(splat_angles(i));
% %     
% %     % Draw circle
% %     t_circle = linspace(0, 2*pi, 40);
% %     x_circle = circle_x + splat_sizes(i)*cos(t_circle);
% %     y_circle = circle_y + splat_sizes(i)*sin(t_circle);
% %     fill(x_circle, y_circle, colors{8}, 'FaceAlpha', 0.3);
% % end
% % 
% % % Mark the center
% % plot(splatter_center(1), splatter_center(2), 'k*', 'MarkerSize', 10);
% % text(splatter_center(1), splatter_center(2)-0.5, 'Splatter Shape', 'HorizontalAlignment', 'center');
% 
% % Add legend
% % legend('Perlin Shape', 'Center', 'Random Polygon', '', 'Random Walk', '', 'Voronoi Shape', 'Points', '', 'Chaotic Shape', '', 'Random Blob', '', 'Amoeba Shape', '', 'Splatter Shape', '');
% 
% % 定义自由空间和障碍物
% 
% % 
% % % 定义多个障碍物
% % obstacle1 = [50 50; 120 50; 120 120; 50 120]; % 方形障碍物
% % obstacle3 = [350 50; 400 80; 380 120; 330 100]; % 梯形障碍物
% % obstacle4 = [100 250; 100 300; 150 320; 170 290; 160 250]; % 五边形(非凸)
% % obstacle5 = [300 250; 350 280; 340 330; 290 340; 270 300]; % 五边形(非凸)
% % obstacle6 = [400 300; 450 300; 450 350; 400 350]; % 小方形
% % obstacle7 = [350 170; 370 180; 380 200; 350 220; 330 200]; % 五边形(凸)
% % % 创建一个"U"形的非凸障碍物
% % obstacle8 = [180 170; 220 170; 220 180; 190 180; 190 210; 220 210; 220 220; 180 220];
% % 
% % obstacles = {obstacle1, obstacle3, obstacle4, obstacle5, obstacle6, obstacle7, obstacle8};
% % 
% % 
% 
% 
% % Generate obstacle data for CGAL Greene Decomposition
% % This script generates vertex data for several irregular shapes
% % The output can be passed to cgal_greene_decomposition function
% 
% % Clear workspace
% clear all;
% close all;
% clc;
% 
% % Set random seed for reproducibility
% rng(123);
% 
% % Create a cell array to store all obstacles
% all_obstacles = cell(5, 1);
% obstacle_names = cell(5, 1);
% 
% % 1. Perlin Noise-based Shape
% shape_center = [5, 5];
% t = linspace(0, 2*pi, 100);
% r = 2 + 0.8*sin(5*t) + 0.5*cos(7*t) + 0.3*sin(11*t) + 0.4*cos(13*t);
% x = r.*cos(t) + shape_center(1);
% y = r.*sin(t) + shape_center(2);
% % Reduce the number of points to a more manageable size
% num_points = 20;
% indices = round(linspace(1, length(x), num_points));
% perlin_shape = [x(indices)', y(indices)'];
% all_obstacles{1} = perlin_shape;
% obstacle_names{1} = 'Perlin Shape';
% 
% % 2. Random Polygon with many vertices
% num_vertices = randi([8, 15]);
% shape_center = [-7, 6];
% r = 1.5 + 1.5*rand(1, num_vertices);
% t = linspace(0, 2*pi, num_vertices+1);
% t = t(1:end-1);
% x = r.*cos(t) + shape_center(1);
% y = r.*sin(t) + shape_center(2);
% random_polygon = [x', y'];
% all_obstacles{2} = random_polygon;
% obstacle_names{2} = 'Random Polygon';
% % 
% % % 3. Voronoi-inspired random shape
% % shape_center = [3, -4];
% % num_points = 20;
% % rand_points_x = shape_center(1) + 4*(rand(1, num_points) - 0.5);
% % rand_points_y = shape_center(2) + 4*(rand(1, num_points) - 0.5);
% % num_hull_points = randi([6, 10]);
% % hull_indices = randperm(num_points, num_hull_points);
% % hull_x = rand_points_x(hull_indices);
% % hull_y = rand_points_y(hull_indices);
% % mean_x = mean(hull_x);
% % mean_y = mean(hull_y);
% % angles = atan2(hull_y - mean_y, hull_x - mean_x);
% % [~, sorted_indices] = sort(angles);
% % hull_x = hull_x(sorted_indices);
% % hull_y = hull_y(sorted_indices);
% % % Perturb the hull points slightly to make it more irregular
% % for i = 1:length(hull_x)
% %     hull_x(i) = hull_x(i) + 0.5*randn();
% %     hull_y(i) = hull_y(i) + 0.5*randn();
% % end
% % voronoi_shape = [hull_x', hull_y'];
% % all_obstacles{3} = voronoi_shape;
% % obstacle_names{3} = 'Voronoi Shape';
% 
% % 4. Blob with random bumps
% shape_center = [-3, 2];
% t = linspace(0, 2*pi, 200);
% bump_count = randi([5, 12]);
% bump_freqs = randi([2, 8], 1, bump_count);
% bump_amps = 0.2 + 0.4*rand(1, bump_count);
% bump_phases = 2*pi*rand(1, bump_count);
% r = ones(size(t))*2;
% for i = 1:bump_count
%     r = r + bump_amps(i)*sin(bump_freqs(i)*t + bump_phases(i));
% end
% x = r.*cos(t) + shape_center(1);
% y = r.*sin(t) + shape_center(2);
% % Reduce the number of points to a more manageable size
% num_points = 25;
% indices = round(linspace(1, length(x), num_points));
% blob_shape = [x(indices)', y(indices)'];
% all_obstacles{3} = blob_shape;
% obstacle_names{3} = 'Blob Shape';
% 
% % 5. Amoeba-like shape with pseudopods
% amoeba_center = [10, -2];
% t = linspace(0, 2*pi, 200);
% r_base = 1.5*ones(size(t));
% num_pods = randi([3, 7]);
% pod_centers = linspace(0, 2*pi, num_pods+1);
% pod_centers = pod_centers(1:end-1);
% pod_widths = 0.3 + 0.5*rand(1, num_pods);
% pod_lengths = 1 + 1.5*rand(1, num_pods);
% r = r_base;
% for i = 1:num_pods
%     pod_factor = pod_lengths(i) * exp(-((t - pod_centers(i)).^2) / (2*pod_widths(i)^2));
%     r = r + pod_factor;
% end
% x = r.*cos(t) + amoeba_center(1);
% y = r.*sin(t) + amoeba_center(2);
% % Reduce the number of points to a more manageable size
% num_points = 25;
% indices = round(linspace(1, length(x), num_points));
% amoeba_shape = [x(indices)', y(indices)'];
% all_obstacles{4} = amoeba_shape;
% obstacle_names{4} = 'Amoeba Shape';
% 
% % 6. Create a completely new shape: A star-crown shape
% star_center = [-8, -5];
% num_points = 9;
% t = linspace(0, 2*pi, num_points+1);
% t = t(1:end-1);
% inner_radius = 1.5;
% outer_radius = 3.5;
% % Create alternating inner and outer points
% r = zeros(1, 2*num_points);
% r(1:2:end) = outer_radius;
% r(2:2:end) = inner_radius;
% t_interleaved = zeros(1, 2*num_points);
% for i = 1:num_points
%     t_interleaved(2*i-1) = t(i);
%     t_interleaved(2*i) = t(i) + pi/num_points;
% end
% x = r.*cos(t_interleaved) + star_center(1);
% y = r.*sin(t_interleaved) + star_center(2);
% star_crown_shape = [x', y'];
% all_obstacles{5} = star_crown_shape;
% obstacle_names{5} = 'Star-Crown Shape';
% 
% % Visualize all obstacles
% figure('Position', [100, 100, 1000, 800]);
% colors = {'b', 'r', 'g', 'm', 'c', 'y'};
% for i = 1:length(all_obstacles)
%     subplot(2, 3, i);
%     points = all_obstacles{i};
%     
%     % Visualize the polygon
%     x = points(:,1);
%     y = points(:,2);
%     % Close the polygon
%     x = [x; x(1)];
%     y = [y; y(1)];
%     
%     plot(x, y, '-', 'Color', colors{i}, 'LineWidth', 2);
%     fill(x, y, colors{i}, 'FaceAlpha', 0.3);
%     
%     % Calculate and plot center point
%     center_x = mean(points(:,1));
%     center_y = mean(points(:,2));
%     hold on;
%     plot(center_x, center_y, 'k*', 'MarkerSize', 10);
%     
%     grid on;
%     title(obstacle_names{i});
%     axis equal;
% end
% 
% % Example of using the function with one of the obstacles
% fprintf('\nExample usage of cgal_greene_decomposition function:\n');
% fprintf('----------------------------------------------\n');
% fprintf('obstacle_points = all_obstacles{1}; %% Select the first obstacle\n');
% fprintf('decomposed_polygons = cgal_greene_decomposition(obstacle_points);\n');
% fprintf('\n');
% fprintf('% Then you can visualize or process the decomposed polygons\n');
% 
% % Save obstacles to a MAT file for later use
% save('obstacles_for_decomposition.mat', 'all_obstacles', 'obstacle_names');
% fprintf('\nObstacle data has been saved to "obstacles_for_decomposition.mat"\n');
% 
% % Test with the function - uncomment to use
% % Test the function with one example


% Generate expanded obstacle data for CGAL Greene Decomposition
% This script generates vertex data for several irregular shapes
% that match the scale of the provided environment

% Clear workspace
clear all;
close all;
clc;
% Generate expanded obstacle data for CGAL Greene Decomposition
% This script includes the original 5 shapes expanded to 2x their size
% alongside the 6 custom obstacles

% Clear workspace
clear all;
close all;
clc;

% Set random seed for reproducibility
rng(123);

% Create a cell array to store all obstacles
all_obstacles = cell(10, 1);
obstacle_names = cell(10, 1);

% Define the state space boundaries for reference
statespace = [-100 -100; 800 -100; 800 600; -100 600]; % The environment boundary

% Generate the original 5 obstacle shapes expanded to 2x their size
% This script creates larger versions of the original obstacles

% Clear workspace
clear all;
close all;
clc;

% Set random seed for reproducibility
rng(123);

% Create a cell array to store all obstacles
all_obstacles = cell(1,4);
obstacle_names = cell(1,4);

% Define the state space boundaries for reference
statespace = [-100 -100; 800 -100; 800 600; -100 600]; % The environment boundary

% 1. Original Perlin Noise-based Shape (expanded 2x)
shape_center = [500, 250]; % Position in environment
scale_factor = 15; % Make it 2x larger than original
t = linspace(0, 2*pi, 100);
r = (2 + 0.8*sin(5*t) + 0.5*cos(7*t) + 0.3*sin(11*t) + 0.4*cos(13*t)) * scale_factor;
x = r.*cos(t) + shape_center(1);
y = r.*sin(t) + shape_center(2);
% Reduce the number of points to a more manageable size
num_points = 20;
indices = round(linspace(1, length(x), num_points));
perlin_shape = [x(indices)', y(indices)'];
all_obstacles{1} = perlin_shape;
obstacle_names{1} = 'Perlin Shape (2x)';

% 2. Original Random Polygon with many vertices (expanded 2x)
num_vertices = randi([8, 15]);
shape_center = [700, 500]; % Position in environment
scale_factor = 15; % Make it 2x larger than original
r = (1.5 + 1.5*rand(1, num_vertices)) * scale_factor;
t = linspace(0, 2*pi, num_vertices+1);
t = t(1:end-1);
x = r.*cos(t) + shape_center(1);
y = r.*sin(t) + shape_center(2);
random_polygon = [x', y'];
all_obstacles{2} = random_polygon;
obstacle_names{2} = 'Random Polygon (2x)';

% % 3. Original Voronoi-inspired random shape (expanded 2x)
% shape_center = [390, 390]; % Position in environment
% scale_factor = 10; % Make it 2x larger than original
% num_points = 20;
% rand_points_x = shape_center(1) + 4*scale_factor*(rand(1, num_points) - 0.5);
% rand_points_y = shape_center(2) + 4*scale_factor*(rand(1, num_points) - 0.5);
% num_hull_points = randi([6, 10]);
% hull_indices = randperm(num_points, num_hull_points);
% hull_x = rand_points_x(hull_indices);
% hull_y = rand_points_y(hull_indices);
% mean_x = mean(hull_x);
% mean_y = mean(hull_y);
% angles = atan2(hull_y - mean_y, hull_x - mean_x);
% [~, sorted_indices] = sort(angles);
% hull_x = hull_x(sorted_indices);
% hull_y = hull_y(sorted_indices);
% % Perturb the hull points slightly to make it more irregular
% for i = 1:length(hull_x)
%     hull_x(i) = hull_x(i) + 0.5*scale_factor*randn();
%     hull_y(i) = hull_y(i) + 0.5*scale_factor*randn();
% end
% voronoi_shape = [hull_x', hull_y'];
% all_obstacles{3} = voronoi_shape;
% obstacle_names{3} = 'Voronoi Shape (2x)';

% 4. Original Blob with random bumps (expanded 2x)
shape_center = [400, 50]; % Position in environment
scale_factor = 15; % Make it 2x larger than original
t = linspace(0, 2*pi, 200);
bump_count = randi([5, 12]);
bump_freqs = randi([2, 8], 1, bump_count);
bump_amps = 0.2 + 0.4*rand(1, bump_count);
bump_phases = 2*pi*rand(1, bump_count);
r = ones(size(t))*2*scale_factor;
for i = 1:bump_count
    r = r + bump_amps(i)*sin(bump_freqs(i)*t + bump_phases(i))*scale_factor;
end
x = r.*cos(t) + shape_center(1);
y = r.*sin(t) + shape_center(2);
% Reduce the number of points to a more manageable size
num_points = 25;
indices = round(linspace(1, length(x), num_points));
blob_shape = [x(indices)', y(indices)'];
all_obstacles{3} = blob_shape;
obstacle_names{3} = 'Blob Shape (2x)';

% 5. Original Amoeba-like shape with pseudopods (expanded 2x)
amoeba_center = [670, 320]; % Position in environment
scale_factor = 15; % Make it 2x larger than original
t = linspace(0, 2*pi, 200);
r_base = 1.5*ones(size(t))*scale_factor;
num_pods = randi([3, 7]);
pod_centers = linspace(0, 2*pi, num_pods+1);
pod_centers = pod_centers(1:end-1);
pod_widths = 0.3 + 0.5*rand(1, num_pods);
pod_lengths = 1 + 1.5*rand(1, num_pods);
r = r_base;
for i = 1:num_pods
    pod_factor = pod_lengths(i) * exp(-((t - pod_centers(i)).^2) / (2*pod_widths(i)^2))*scale_factor;
    r = r + pod_factor;
end
x = r.*cos(t) + amoeba_center(1);
y = r.*sin(t) + amoeba_center(2);
% Reduce the number of points to a more manageable size
num_points = 25;
indices = round(linspace(1, length(x), num_points));
amoeba_shape = [x(indices)', y(indices)'];
all_obstacles{4} = amoeba_shape;
obstacle_names{4} = 'Amoeba Shape (2x)';

% Visualize all obstacles along with reference obstacles
figure('Position', [100, 100, 1200, 800]);

% First plot the state space
subplot(1, 1, 1);
hold on;
% Plot statespace boundary
plot([statespace(:,1); statespace(1,1)], [statespace(:,2); statespace(1,2)], 'k-', 'LineWidth', 2);

% Define some colors
colors = {'b', 'r', 'g', 'm', 'c'};

% Plot reference obstacles in light gray
reference_obstacles = get_reference_obstacles();
for i = 1:length(reference_obstacles)
    obs = reference_obstacles{i};
    fill([obs(:,1); obs(1,1)], [obs(:,2); obs(1,2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', [0.6 0.6 0.6]);
end

% Plot our custom obstacles
for i = 1:length(all_obstacles)
    points = all_obstacles{i};
    % Close the polygon for plotting
    x = [points(:,1); points(1,1)];
    y = [points(:,2); points(1,2)];
    
    fill(x, y, colors{i}, 'FaceAlpha', 0.5);
    plot(x, y, '-', 'Color', colors{i}, 'LineWidth', 2);
    
    % Calculate and plot center point
    center_x = mean(points(:,1));
    center_y = mean(points(:,2));
    plot(center_x, center_y, 'k*', 'MarkerSize', 10);
    text(center_x, center_y+10, obstacle_names{i}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end


axis equal;
xlim([-120 820]);
ylim([-120 620]);
xlabel('X');
ylabel('Y');

% Create separate detailed plots for each obstacle
figure('Position', [100, 100, 1000, 800]);
for i = 1:length(all_obstacles)
    subplot(2, 3, i);
    points = all_obstacles{i};
    % Close the polygon for plotting
    x = [points(:,1); points(1,1)];
    y = [points(:,2); points(1,2)];
    
    fill(x, y, colors{i}, 'FaceAlpha', 0.5);
    plot(x, y, '-', 'Color', colors{i}, 'LineWidth', 2);
    
    % Plot vertices
    plot(points(:,1), points(:,2), 'k.', 'MarkerSize', 10);
    
    % Calculate and plot center point
    center_x = mean(points(:,1));
    center_y = mean(points(:,2));
    plot(center_x, center_y, 'k*', 'MarkerSize', 10);
    
    grid on;
    title([obstacle_names{i}, ' (', num2str(size(points, 1)), ' vertices)']);
    axis equal;
end


% 
% test_all_decompositions(all_obstacles, obstacle_names)
% 
% 
% 
% 
% figure;
% for i = 1:length(all_obstacles)
% obstacle_points = all_obstacles{i};
% 
% try
%     decomposed_polygons = cgal_greene_decomposition(obstacle_points);
%     fprintf('Successfully decomposed obstacle %s\n', obstacle_names{i});
%     
%     % Visualize decomposed polygons
% 
%     hold on;
%     title(['Decomposition of ', obstacle_names{i}]);
%     
%     for i = 1:length(decomposed_polygons)
%         poly = decomposed_polygons{i};
%         fill(poly(:,1), poly(:,2), rand(1,3), 'FaceAlpha', 0.5);
%     end
%     
%     % Original shape outline
%     plot([obstacle_points(:,1); obstacle_points(1,1)], ...
%          [obstacle_points(:,2); obstacle_points(1,2)], 'k-', 'LineWidth', 2);
%     
%     axis equal;
%     grid on;
% catch e
%     fprintf('Error in decomposition: %s\n', e.message);
% end
% 
% end
% 
% all_obstacles = [all_obstacles,reference_obstacles];
% [freeSpace, V, M] = createFreeSpaceHZ(statespace, all_obstacles(1,1:14));
% plot(freeSpace)

% 
% 
% % Function to test all obstacles with the decomposition function
% function test_all_decompositions(all_obstacles, obstacle_names)
%     figure('Position', [100, 100, 1200, 800]);
%     rows = ceil(sqrt(length(all_obstacles)));
%     cols = ceil(length(all_obstacles) / rows);
%     successful = 0;
%     failed = 0;
%     
%     for i = 1:length(all_obstacles)
%         subplot(rows, cols, i);
%         hold on;
%         
%         obstacle_points = all_obstacles{i};
%         try
%             decomposed_polygons = cgal_greene_decomposition(obstacle_points);
%             fprintf('Successfully decomposed obstacle %s\n', obstacle_names{i});
%             successful = successful + 1;
%             
%             % Visualize decomposed polygons
%             colors = jet(length(decomposed_polygons));
%             for j = 1:length(decomposed_polygons)
%                 poly = decomposed_polygons{j};
%                 fill(poly(:,1), poly(:,2), colors(j,:), 'FaceAlpha', 0.5);
%             end
%             
%             % Original shape outline
%             plot([obstacle_points(:,1); obstacle_points(1,1)], ...
%                  [obstacle_points(:,2); obstacle_points(1,2)], 'k-', 'LineWidth', 2);
%             
%             title([obstacle_names{i}, ' (', num2str(length(decomposed_polygons)), ' parts)']);
%         catch e
%             fprintf('Error decomposing %s: %s\n', obstacle_names{i}, e.message);
%             failed = failed + 1;
%             
%             % Just show the original obstacle
%             fill([obstacle_points(:,1); obstacle_points(1,1)], ...
%                  [obstacle_points(:,2); obstacle_points(1,2)], 'r', 'FaceAlpha', 0.3);
%             
%             title([obstacle_names{i}, ' (FAILED)']);
%         end
%         
%         axis equal;
%         grid on;
%     end
%     
%     fprintf('\nDecomposition summary: %d successful, %d failed\n', successful, failed);
% end

% Function to get reference obstacles
function obstacles = get_reference_obstacles()
    % 1. Basic square obstacle
    obstacle1 = [50 50; 150 50; 150 150; 50 150];
    % 2. Left vertical wall of maze
    obstacle2 = [250 50; 270 50; 270 350; 250 350];
    % 3. Top horizontal wall of maze
    obstacle3 = [270 330; 580 330; 580 350; 270 350];
    % 4. Right vertical wall of maze
    obstacle4 = [560 150; 580 150; 580 330; 560 330];
    % 5. Middle horizontal wall (with gap)
    obstacle5 = [270 150; 450 150; 450 170; 270 170];
    obstacle6 = [470 150; 560 150; 560 170; 470 170];
    % 6. Circle as polygon approximation
    circle_center = [180 250];
    circle_radius = 40;
    circle_points = 7;
    obstacle7 = [];
    for i = 1:circle_points
        angle = 2*pi*(i-1)/circle_points;
        obstacle7 = [obstacle7; circle_center + circle_radius*[cos(angle), sin(angle)]];
    end
    % 7. U-shaped obstacle
    obstacle8 = [650 50; 750 50; 750 220; 700 220; 700 100; 650 100];
    % 8. Star shape
    star_center = [580 450];
    star_outer = 80;
    star_inner = 40;
    star_points = 5;
    obstacle9 = [];
    for i = 1:2*star_points
        angle = pi*(i-1)/star_points;
        if mod(i, 2) == 1
            obstacle9 = [obstacle9; star_center + star_outer*[cos(angle), sin(angle)]];
        else
            obstacle9 = [obstacle9; star_center + star_inner*[cos(angle), sin(angle)]];
        end
    end
    % 9. Narrow channels
    obstacle10 = [50 400; 200 400; 200 420; 50 420];
    obstacle11 = [50 480; 200 480; 200 500; 50 500];
    % 10. Diamond
    obstacle12 = [350 400; 380 450; 350 500; 320 450];
    % 11. Random polygon
    obstacle13 = [400 420; 440 400; 480 420; 460 470; 420 480; 390 450];

    obstacles = {obstacle1, obstacle2, obstacle3, obstacle4, obstacle5, obstacle6, ...
        obstacle7, obstacle8, obstacle9, obstacle10, obstacle11, obstacle12, obstacle13};
end


