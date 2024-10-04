clear;
% Some parameters
N = 200;
M = N;
L = 5;
H = L;
c = [L/2, H/2];
r = 0.5;
contours = 5;


% Data input and preprocessing
data = readtable('output_final/200_output.csv');
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
contourf(x_grid,y_grid, Cp_grid,contours*2); % quiver plotting (input positions and velocities)


% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure Coefficient (Cp)';
colormap spring

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
contourf(x_grid,y_grid, T_grid,contours); % quiver plotting (input positions and velocities)

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Temperature [K]';
colormap spring

% Cylinder drawing
phi = linspace(0, 2*pi);
x_r = r*cos(phi) + c(1);
y_r = r*sin(phi) + c(2);
patch(x_r,y_r,'black');

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Temperature field','Interpreter','latex');
grid on
axis equal

%% %% P PLOT 
figure(8)
[x_grid,y_grid] = meshgrid(linspace(0,L,M),linspace(0,H,N)); 
p_grid = griddata(X, Y, p ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
contourf(x_grid,y_grid, p_grid,contours); % quiver plotting (input positions and velocities)

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Pressure [Pa]';
colormap spring

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

%% RESULTS CONVERGENCE
L = 0;
H = 0;
R = 0;
N = 0;
M = 0;
Cl = 0;
Cd = 0;
Circ = 0;

for n=[15,25,50,75,100,150,200]
   name = 'output_final\';
   elems = sprintf('%i_results.csv',n);
   data = table2array(readtable(append(name,elems)));
   L = [L, data(1)];
   H = [H, data(2)];
   R = [R, data(3)];
   N = [N, data(4)];
   M = [M, data(5)];
   Cl = [Cl, data(6)];
   Cd = [Cd, data(7)];
   Circ = [Circ, data(8)];
    
end
figure(9)
plot(N, Cl)
xlabel('Node number');
ylabel('Lift coefficient (Cl)');
title('Lift coefficient vs mesh nodes','Interpreter','latex');

figure(10)
plot(N, Cd)
xlabel('Node number');
ylabel('Drag coefficient (Cd)');
title('Drag coefficient vs mesh nodes','Interpreter','latex');

figure(11)
plot(N, Circ)
xlabel('Node number');
ylabel('Circulation (m/s)');
title('Circulation vs mesh nodes','Interpreter','latex');

