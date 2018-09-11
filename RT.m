clear
close all
clc

%%  Code for plano-convex lens

n           = 1.46; %Index of refraction of lens
radius      = 20;   %Radius of spherical surface
thickness   = 5; %Central thickness of lens
dz          = 0.01; %Step size for computation purposes
aperture    = 8;
number_rays = 30;
dy          = (2*aperture + 1)/number_rays;
y           = -aperture:dy:aperture; %Field of view

%% Ray matrix
[raymatrix,x_front,x_optaxis,zmax] = ...
    plano_convex(n,radius,thickness,dz,y);


%% Figure

front_lens = sqrt(radius^2 - (x_front-radius).^2); 
figure(1)
plot(40-x_optaxis,raymatrix','r') %Rays
hold on

%Plano-convex lens
line([40-thickness 40-thickness], [max(front_lens) -max(front_lens)],'color','black')
plot(40-x_front,front_lens,'black',40-x_front,-front_lens,'black') %Lens front surface 
plot(x_optaxis*5,zeros(1,length(x_optaxis)),'k--') %Optical axis
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
x_rm = [200 200];
y_rm = [-15 15];
plot(x_rm, y_rm, 'black')

%Upper mirror
x_um = [135 165];
y_um = [50 50];
plot(x_um, y_um, 'black');

%Transmitted rays
x_2 = x_rm(1,1).*ones(1,29);
x_trans = [x_1'; x_2];
y_trans = y_ray;
plot (x_trans, y_trans, 'b');

%Reflected rays
y_3 = y_um(1,1).*ones(1,29); 
x_reflect = [x_1'; x_1'];
y_reflect = [raymatrix(:,1)'; y_3];
plot (x_reflect, y_reflect, 'b');

%Interference rays
detektor_posx = [140 160];
detektor_posy = [-100 -100];
plot(detektor_posx, detektor_posy, 'black')

y_4 = detektor_posy(1)*ones(1,29);
x_interference = x_reflect;
y_interference = [raymatrix(:,1)'; y_4];
plot (x_interference, y_interference, 'g');


hold on
grid on
axis equal

%% Intensitätsverteilung Theorie (Wave parameters)

f0 = 5; % frequency
A = 2; % amplitude
lambda = 468; %nm
k = (2*pi)/lambda; 


maxWaves = 10; % number of waves %erstmal nur 10 aber eigentlich 29!!
 
maxTime = 5;
maxCount = 1000;
x = 10^(-4):10^(-4):0.1;

% Generate waves
t = linspace(-f0/50,f0/50,maxCount); %nm
df=0;
for waves = 1: maxWaves
       df = df+10;
        w = 2*pi*(f0+df);
I(waves,:)= (A*cos(w.*t+k.*x)).^2;
end

I_sum = sum(I)-20;
% Plot
figure(2);
%{
for count = 1: maxWaves
    plot(t,I(count,:));
    hold on;
end
%}
plot(t,I_sum);
hold off;
grid on
%set(gca, FontSize, 9, FontWeight, Bold, LineWidth, 1);
xlabel('Time (s)');
ylabel('Intensity');


%% Intensitätsverteilung


for i = 0:200        
    if i <= x_line 
        int = I;     % intensity of ray 
    else 
        int = 0.5*I;
    end
end


%%  Schnittpunkt: Linse 

% Eingabe radius, y - Höhe, x_pos

% Höhe ermitteln:
x_pos   = 35; % Position der Linse
sag     = radius - sqrt(radius.^2-y.^2); % Verlauf konvexe Seite
x_linse = (x_pos + thickness) - sag; % Höhenwerte der Schnittpunkte


sinphi1 = y/radius;   % Winkel in der Linse
sinphi2 = sinphi1/n;
phi1    = asin(sinphi1);
phi2    = asin(sinphi2);
theta   = phi1 - phi2;

G       = thickness - sag; % Länge der Ak
y_diff  = G.*tan(theta); % Diff zwischen Höhe am Ende der Linse und Anfang
y_Schni = y - y_diff;

%Schnittpunkt bei 
%figure(1)
%plot (x_pos*ones(1,29),y_Schni, 'og')
%hold on

% Winkel ermitteln:
x_Anfang = -5; % Strahlenanfang

beta1 = atan(abs(y_Schni)/(x_pos-x_Anfang)); %rad
beta2 = beta1*180/pi; %deg


%% Schnittpunkt: Strahlenteiler

x_Schnitt_ST = deltax*ones(1,29) + y;
y_Schnitt_ST = y;
%plot(x_Schnitt_ST, y_Schnitt_ST, 'og')

%% Schnittpunkt: rechter Spiegel

x_Schnitt_rm = (x_rm(1)).*ones(1,29);
y_Schnitt_rm = y;
%plot(x_Schnitt_rm, y_Schnitt_rm, 'og')

%% Schnittpunkt: oberer Spiegel

x_Schnitt_um = x_Schnitt_ST;
y_Schnitt_um = (y_um(1)).*ones(1,29);
%plot(x_Schnitt_um, y_Schnitt_um, 'og')

%% Schnittpunkt: Detektor

x_Schnitt_D = x_Schnitt_ST;
y_Schnitt_D = y_4;
%plot(x_Schnitt_D, y_Schnitt_D, 'og')

%% Rays

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

zmax = floor(f+0.035*f);
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
