clear;
% Some parameters
N=100;

% Data input and preprocessing
data = readtable('output/output.csv');
X = table2array(data(:,1));
Y = table2array(data(:,2));
U = table2array(data(:,3));
V = table2array(data(:,4));
S = table2array(data(:,5));
rho = table2array(data(:,6));
Cp = table2array(data(:,7));
solid = table2array(data(:,8));

data_a = readtable('output/analytic_output.csv');
X_a = table2array(data_a(:,1));
Y_a = table2array(data_a(:,2));
S_a = table2array(data_a(:,5));
rho_a = table2array(data_a(:,6));
Cp_a = table2array(data_a(:,7));
solid_a = table2array(data_a(:,8));




%% %% QUIVER PLOT 
figure(1)
quiver(X, Y, U, V); % quiver plotting (input positions and velocities)

% Cylinder drawing
c = [0.5, 0.5];
r = 0.15;
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Velocity field','Interpreter','latex');
grid on
axis equal


%% %% CONTOUR PLOT
figure(2)
%Contour plotting
[x_grid,y_grid] = meshgrid(linspace(0,1,N),linspace(0,1,N)); 
p_grid = griddata(X, Y, S ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contour(x_grid,y_grid,p_grid,'LevelList',linspace(min(S),max(S),50))

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure Coefficient (Cp)';

% Cylinder drawing
c = [0.5, 0.5];
r = 0.15;
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Stream lines ($\psi$)','Interpreter','latex');
grid on
axis equal
colormap cool


%% %% ANALYTIC CONTOUR PLOT
figure(3)
%Contour plotting
[x_grid_a,y_grid_a] = meshgrid(linspace(0,1,N),linspace(0,1,N)); 
p_grid_a = griddata(X_a, Y_a, S_a ,x_grid_a,y_grid_a); %interpolates surface from  mesh and streamline values (cubic interpolation)
contour(x_grid_a,y_grid_a,p_grid_a,'LevelList',linspace(min(S_a),max(S_a),50))

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure Coefficient (Cp)';

% Cylinder drawing
c = [0.5, 0.5];
r = 0.15;
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Stream lines ($\psi$)','Interpreter','latex');
grid on
axis equal
colormap cool

%% %% CYLINDER PLOT
figure(5)
% Grid of solid nodes (binary)
k=0;
solid_mesh = zeros(N-2,N-2);
for i=1:N-2
    for j=1:N-2
        k=k+1;
        solid_mesh(i,j)=solid(k,1);
    end
end

% Plot mesh as a binary file
imshow(solid_mesh, 'InitialMagnification', 'fit');
colormap(gca, [1 1 1; 0 0 0]); % Set the colormap to white and black

% Plot parameters
title('Mesh nodes');
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
grid on

