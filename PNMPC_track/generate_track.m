clc
close all
clear all

% Example points (replace these with your actual points)
points = [

0,0;
0.1,0;
0.2,0;
0.3,0;
0.4,0;
0.5,0;
0.6,0;
0.7,0;
0.8,0;
1,0;
2,0.5;
3,1;
3,2;
2,2.5;
1,3;
0,3;
-1,2.5;
-2,2;
-2,1;
-1,0.5;
0,0;

];

% Initialize theta array
theta = zeros(size(points, 1)-1, 1);
global_theta = zeros(size(points, 1)-1, 1);

% Calculate theta for each segment
for i = 1:length(points)-1
    dx = points(i+1, 1) - points(i, 1);
    dy = points(i+1, 2) - points(i, 2);
    theta(i) = atan2d(dy, dx); % atan2d returns the angle in degrees and handles all quadrants
end

% Calculate gloabl theta
global_theta(1) = theta(1);
for i = 2:length(theta)
    delta_theta = theta(i) - theta(i-1);
    if delta_theta > 180
        delta_theta = delta_theta - 360;
    elseif delta_theta < -180
        delta_theta = delta_theta + 360;
    end
    global_theta(i) = global_theta(i-1) + delta_theta;
end

% Plotting the points and directions
figure;
hold on;
grid on;
axis equal;

% Plot points
plot(points(:,1), points(:,2), 'bo-');

% Plot directions
for i = 1:length(points)-1
    quiver(points(i, 1), points(i, 2), cosd(theta(i)), sind(theta(i)), 0.5, 'r', 'LineWidth', 1, 'MaxHeadSize', 1);
    text(points(i, 1), points(i, 2), sprintf('%.1f°', global_theta(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

Track_x = points(2:end,1);
Track_y = points(2:end,2);
Track_theta = deg2rad(global_theta);

save('track.mat','Track_x','Track_y','Track_theta');