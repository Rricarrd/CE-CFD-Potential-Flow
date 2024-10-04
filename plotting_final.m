clear;
% Some parameters
N = 300;
M = N;
L = 5;
H = L;
c = [L/2, H/2];
r = 0.5;


% Data input and preprocessing
data = readtable('output_final/50_output.csv');
X = table2array(data(:,1));
Y = table2array(data(:,2));
U = table2array(data(:,3));
V = table2array(data(:,4));
S = table2array(data(:,5));
rho = table2array(data(:,6));
Cp = table2array(data(:,7));
solid = table2array(data(:,8));
p = table2array(data(:,9));
T = table2array(data(:,10));


%% %% QUIVER PLOT 
figure(1)
quiver(X, Y, U, V); % quiver plotting (input positions and velocities)

% Cylinder drawing
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
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
p_grid = griddata(X, Y, S ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contour(x_grid,y_grid,p_grid,'LevelList',linspace(min(S),max(S),50))

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Streamlines [m2/s]';

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Numerical stream lines ($\psi$)','Interpreter','latex');
grid on
axis equal
colormap cool





% %% %% CYLINDER PLOT
% figure(5)
% % Grid of solid nodes (binary)
% k=0;
% solid_mesh = zeros(N-2,N-2);
% for i=1:N-2
%     for j=1:N-2
%         k=k+1;
%         solid_mesh(i,j)=solid(k,1);
%     end
% end
% 
% % Plot mesh as a binary file
% imshow(solid_mesh, 'InitialMagnification', 'fit');
% colormap(gca, [1 1 1; 0 0 0]); % Set the colormap to white and black
% 
% % Plot parameters
% title('Mesh nodes');
% xlabel('X-axis [m]');
% ylabel('Y-axis [m]');
% grid on

%% %% CP PLOT 
figure(6)
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
Cp_grid = griddata(X, Y, Cp ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contourf(x_grid,y_grid, Cp_grid); % quiver plotting (input positions and velocities)


% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure Coefficient (Cp)';

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Cp');
grid on
axis equal

%% %% T PLOT 
figure(7)
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
T_grid = griddata(X, Y, T ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contourf(x_grid,y_grid, T_grid); % quiver plotting (input positions and velocities)

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Temperature [K]';

% Cylinder drawing
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

%% %% P PLOT 
figure(8)
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
p_grid = griddata(X, Y, p ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contourf(x_grid,y_grid, p_grid); % quiver plotting (input positions and velocities)

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure [Pa]';

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Pressure field','Interpreter','latex');
grid on
axis equal