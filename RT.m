clear
close all
clc

%%  Code for plano-convex lens

n           = 1.5168; %Index of refraction of lens
radius      = 20; %Radius of spherical surface
thickness   = 5; %Central thickness of lens
dz          = 0.01; %Step size for computation purposes
aperture    = 8;
number_rays = 30;
dy          = (2*aperture + 1)/number_rays;
y           = -aperture:dy:aperture; %Field of view

%% Ray matrix
[raymatrix,x_front,x_optaxis,zmax] = ...
    plano_convex(n,radius,thickness,dz,y);


%% Figure1

front_lens = sqrt(radius^2 - (x_front-radius).^2); 
figure(1)
plot(40-x_optaxis,raymatrix','r') %Rays
hold on

%Plano-convex lens
line([40-thickness 40-thickness], [max(front_lens) -max(front_lens)],'color','black')
plot(40-x_front,front_lens,'black',40-x_front,-front_lens,'black') %Lens front surface 
plot(x_optaxis,zeros(1,length(x_optaxis)),'k--') %Optical axis
hold off


%Semi-transparent mirror
y_line = [-15 15];
m = sin(pi/4)/cos(pi/4);
deltax = 150;
x_line = (y_line/m)+deltax;
%Reflector_Enlarger = 5;
% line([(150-(thickness+Reflector_Enlarger)) (150+(thickness+Reflector_Enlarger))], sin(pi/4)*[-10-Reflector_Enlarger 10+Reflector_Enlarger],'color','black')
line(x_line, y_line,'color','black') 
hold on


%Rays in front of semi-transparent lens
y_0 = raymatrix(:,1); %beginning
x_0 = 40*ones(1,29);
y_1 = raymatrix(:,1); %end
x_1 = (y_1/m)+deltax;

x_ray = [x_0; x_1'];
y_ray = [raymatrix(:,1)'; raymatrix(:,1)'];
plot(x_ray, y_ray, 'r')

%Right mirror
x_2 = [200 200];
y_2 = [-15 15];
plot(x_2, y_2, 'black')

%Upper mirror
x_3 = [135 165];
y_3 = [50 50];
plot(x_3, y_3, 'black');

%Transmitted rays
x_2 = 200*ones(1,29);
x_trans = [x_1'; x_2];
y_trans = y_ray;
plot (x_trans, y_trans, 'b');

%Reflected rays
y_3 = 50*ones(1,29); 
x_reflect = [x_1'; x_1'];
y_reflect = [raymatrix(:,1)'; y_3];
plot (x_reflect, y_reflect, 'b');

%Interference rays
y_4 = -100*ones(1,29);
x_interference = x_reflect;
y_interference = [raymatrix(:,1)'; y_4];
plot (x_interference, y_interference, 'g');


hold on
grid on
axis equal

%% Intensitätsverteilung Theorie (Wave parameters)


E0 = 1; %amplitude of the electric field

lambda = 10^(-6); %wavelength in [m]
k = (2*pi)/lambda; %wave number

x = 0; %point 
f = 1; 
omega = 2*pi*f; %angular frequency


maxTime = 2;
maxCount = 1000;
count = 1:maxCount;
t = count/maxCount*maxTime; %time

delta_phi = 0.5*pi; %phase shift 
q_1 = 1:29;
q_2 = (delta_phi*ones(29,1));
phi = q_1'.*q_2; %phase


E = E0*cos(omega*t+phi); %electric field
V = sum(E);

I = (V).^2; %intensity


figure (2)

plot(t,I, 'b'); 
hold on;
plot(t,E, 'r');
hold off
axis ([0 maxTime -1 4]);
grid on
xlabel ('Time (s)');
ylabel ('Intensity');


%%

function ray_air = plane_refract_ray(y, slope, thickness, n, z)
theta1  = atan(slope);
theta2  = asin(n*sin(theta1));
slope2  = tan(theta2);
ray_air = -(thickness-z)*slope2 + y;
end
 
%  Refraction spherical interface
function [ray, slope, z] = sphere_refract_ray(y, radius, thickness, n, dz)
sag     = radius - sqrt(radius^2-y^2);
z       = sag:dz:thickness;
sinphi1 = y/radius;
sinphi2 = sinphi1/n;
phi1    = asin(sinphi1);
phi2    = asin(sinphi2);
theta   = phi2-phi1;
slope   = tan(theta);
ray     = -slope*(sag-z) + y;
end

function [raymatrix, z_front, z_optaxis, zmax] = plano_convex(...
    n, radius, thickness, dz, y)

power = (n-1)/radius;
f     = 1/power;

zmax = floor(f+0.1*f);
z_front = 0:dz:thickness-dz;
z_back = thickness:dz:zmax-dz;
z_optaxis = [z_front, z_back];
y( y == 0 ) = 10^(-10);
raymatrix   = zeros(length(y), length(z_optaxis));

%  Ray Tracing
for i = 1:length(y)
    
    % Refraction spherical surface
    [ray_lens, slope, x_lens] = sphere_refract_ray(y(i), radius, ...
        thickness, n, dz);
    
    % Refraction plane surface
    ray_air = plane_refract_ray(ray_lens(end), slope, thickness, ...
        n, z_back);
    
    % Incoming ray
    x_front_air = 0:dz:x_lens(1)-dz;
    ray_front_air = y(i)*ones(1, length(x_front_air));
    
    % Create matrix of rays
    if length(ray_lens)+length(ray_air)+length(x_front_air) <= length(z_optaxis)
        raymatrix(i,:) = ...
            [ray_front_air, ray_lens, ray_air];
    else
        raymatrix(i,:) = ...
            [ray_front_air, ray_lens(1:length(ray_lens)-1), ray_air];
    end
end
end


