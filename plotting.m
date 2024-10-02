clear;
%output = readtable('output/output_neg_circ_n100.csv');
output = readtable('output/output_static_n100.csv');
N=100;

X = table2array(output(:,1));
Y = table2array(output(:,2));
U = table2array(output(:,3));
V = table2array(output(:,4));
S = table2array(output(:,5));
rho = table2array(output(:,6));
Cp = table2array(output(:,7));
solid = table2array(output(:,8));



figure(1)
quiver(X, Y, U, V)

xlabel('X-axis');
ylabel('Y-axis');
title('Velocity field');
grid on
axis equal

hold on

[xq,yq] = meshgrid(...
  linspace(0,1,100),...
  linspace(0,1,100));

pq = griddata(X, Y,S,xq,yq,'cubic');
% Visualize the result

contour(xq,yq,pq,20)
c = colorbar;
c.Label.String = 'Pressure Coefficient';
center = [0.5, 0.5];
radius = 0.15;
rectangle('Position', [center(1) - radius, center(2) - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', 'black');

%rectangle('Position', [center(1) - radius, center(2) - radius, 2 * radius, 2 * radius], 'Curvature', [1, 1], 'FaceColor', 'black');
%grid on
%axis equal

figure(3)
c=0;

for i=1:N-2
    for j=1:N-2
        c=c+1;
        meshData(i,j)=solid(c,1);
    end
end

% Create a binary image where false values are white (0) and true values are black (1)
binaryImage = ~~meshData; % Use the ~ operator to invert the values

% Display the binary image with white for false and black for true
imshow(binaryImage, 'InitialMagnification', 'fit');
colormap(gca, [1 1 1; 0 0 0]); % Set the colormap to white and black
colorbar; % Optional: Display a colorbar to show the mapping

% Optional: Add title and axis labels
title('Binary Image of Mesh Data');
xlabel('X-axis');
ylabel('Y-axis');

