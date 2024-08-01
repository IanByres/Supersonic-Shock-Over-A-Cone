clc; clear; close all; tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           ARO3111 Computer Assignment
%        Supersonic Flow Past Un-Yawed Cone
%                 by Ian Byres
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
The following program takes the free-stream Mach number and the conical
shock angle to estimate flowfield properties after the shock, as well as
the cone angle that corresponds to the given shock and free-stream
conditions. This program uses the inverse numerical approach described in
Chapter 10 of John D. Anderson's Modern Compressible Flow, 4th Edition.
%}

%% Inputs (EDIT HERE)

M1 = 6; % Upstream Mach number
k = 1.4; % Ratio of specific heats (assume constant throughout)
thetaS = 12; % Shock Angle in [deg]
delTheta = -0.1; % Iteration increment size [deg] (must be negative)

%% ===== Initial setup (DO NOT EDIT PAST HERE) =====

%{ 
The following code computes the Mach number immediately after the shock
wave away from the cone surface using the locally planar approximation
where the shock components can be evaluated using the 2D oblique shock
relations on the ray immediately after the shock.
%}

thetaS = deg2rad(thetaS);
delTheta = deg2rad(delTheta);

Mn1 = M1*sin(thetaS); % Normal upstream Mach number to shock wave
Mn2 = sqrt((1+((k-1)/2)*(Mn1^2))/(k*Mn1^2 - (k-1)/2)); % Normal Mach number after shock away from cone
deltaShock = atan( (2*(M1^2 * sin(thetaS)^2 - 1))/(tan(thetaS)*(2+(M1^2)*(k+cos(2*thetaS)))) ); % Flow turn angle of shock
M2 = Mn2/sin(thetaS - deltaShock); % Mach number after shock away from cone surface

rP2P1 = ( 2*k*Mn1^2 - (k-1) )/(k+1); % Ratio of P2 to freestream pressure
rRho2Rho1 = ( (k+1)*Mn1^2 )/( (k-1)*Mn1^2 + 2 ); % Ratio of Rho2 to freestream density
rT2T1 = ( (2*k*Mn1^2 - (k-1))*((k-1)*Mn1^2 + 2) )/( (k+1)^2 * Mn1^2 ); % Ratio of T2 to freestream temperature

rP2Ptot2 = (1 + ((k-1)/2)*M2^2)^(-k/(k-1)); % Ratio of P2 to stagnation pressure
rRho2Rhotot2 = (1 + ((k-1)/2)*M2^2)^(-1/(k-1)); % Ratio of Rho2 to stagnation density
rT2Ttot2 = (1 + ((k-1)/2)*M2^2)^(-1); % Ratio of T2 to stagnation temperature

Vs = sqrt((((k-1)/2)*M2^2)/(1+((k-1)/2)*M2^2)); % Nondimensionalized velocity after shock
Vrs = Vs*cos(thetaS-deltaShock); % Nondimensionalized radial component of velocity after shock
Vts = -Vs*sin(thetaS-deltaShock); % Nondimensionalized angular component of velocity after shock

%% ===== Integrator Called =====

%{ 
Creates integrating span of angles and initial condition matrix for the
integrator. Calls ode45 as Runge-Kutta integrator with event at boundary of
cone to stop integration when Vt = 0 (we have reached the cone surface).
%}

itSpan = thetaS:delTheta:0; % Set of angles to integrate over
initialConditions = [Vrs, Vts]; % Initial Conditions for integrator
options = odeset('Events',@boundaryCheck); % Event function call for integrator

[theta, X] = ode45(@(theta,x) TaylorMaccoll(theta,x,k),itSpan,initialConditions,options); % Calls runge-kutta integrator

%% ===== Results Manipulation =====

%{ 
In this section the integrated results are manipulated to find flowfield
properties and other results.
%}

V = sqrt(X(:,1).^2 + X(:,2).^2); % Nondimensionalized velocities on each ray
M = (sqrt(2)*V)./sqrt( k - 1 - k*V.^2 + V.^2 ); % Nondimensionalized Mach numbers on each ray
rPP1 = ((1 + ((k-1)/2).*M.^2).^(-k/(k-1))) .* (rP2Ptot2)^-1 .* rP2P1; % Pressure ratio to free-stream on each ray
rRhoRho1 = ((1 + ((k-1)/2).*M.^2).^(-1/(k-1))) .* (rRho2Rhotot2)^-1 .* rRho2Rho1; % Density ratio to free-stream on each ray
rTT1 = ((1 + ((k-1)/2).*M.^2).^(-1)) .* (rT2Ttot2)^-1 .* rT2T1; % Temperature ratio to free-stream on each ray

freeStrAngle = zeros(2,1); % Creating additional data points for plots
freeStrCond = ones(2,1); % Creating additional data points for plots
freeStrMach = [M1;M1];
for i = 1:2 % Creating additional data points for plots
    freeStrAngle(i) = thetaS + deg2rad(2-i);
end
V = [freeStrCond;V];
M = [freeStrMach;M];
rPP1 = [freeStrCond;rPP1];
rRhoRho1 = [freeStrCond;rRhoRho1];
rTT1 = [freeStrCond;rTT1];
theta = [freeStrAngle;theta];

thetaCone = theta(end); % Computed cone angle corresponding to shock conditions
Mc = M(end); % Computed Mach number on cone surface
Cpc = (rPP1(end) - 1)/((k/2)*M1^2); % Pressure coefficient on cone surface


%% ===== Displaying Results =====

%{
Here I display the results of the program as plots of Mach number and
property ratios along each ray of the flowfield past the shock.
%}

toc;
fprintf("\n");
fprintf("Freestream Mach: %1g \n",M1);
fprintf("Ratio of Specific Heats: %1g \n",k);
fprintf("Shock Angle [deg]: %1g \n \n",rad2deg(thetaS));
fprintf("Cone Angle [deg]: <strong>%1g</strong> \n",rad2deg(thetaCone));
fprintf("Mach Number Along Cone Surface: <strong>%1g</strong> \n",Mc);
fprintf("Pressure Coefficient on Cone Surface: <strong>%1f</strong> \n",Cpc);


f = figure(1); hold on;
f.WindowState = 'maximized';

min = 0.95;
mx = 1.05;

subplot(2,2,1);
plot(rad2deg(theta),M);
set (gca,'XDir','reverse');
title("Mach Number vs Angle");
ylabel("Mach Number");
xlabel("Angle Theta [deg]");
ylim([min (mx*max(M))]);

subplot(2,2,2);
plot(rad2deg(theta),rPP1);
set (gca,'XDir','reverse');
title("Pressure Ratio P/P1 vs Angle");
ylabel("Pressure Ratio P/P1");
xlabel("Angle Theta [deg]");
ylim([min (mx*max(rPP1))]);

subplot(2,2,3);
plot(rad2deg(theta),rRhoRho1);
set (gca,'XDir','reverse');
title("Density Ratio Rho/Rho1 vs Angle");
ylabel("Density Ratio Rho/Rho1");
xlabel("Angle Theta [deg]");
ylim([min (mx*max(rRhoRho1))]);

subplot(2,2,4);
plot(rad2deg(theta),rTT1);
set (gca,'XDir','reverse');
title("Temperature Ratio T/T1 vs Angle");
ylabel("Temperature Ratio T/T1");
xlabel("Angle Theta [deg]");
ylim([min (mx*max(rTT1))]);


function X = TaylorMaccoll(theta,x,k)
    %{
    The following function is the state vector and its derivative from the
    Taylor-Maccoll equation.
    %}
    Vr = x(1);
    Vt = x(2);

    eq = ((Vt^2)*Vr - ((k-1)/2)*(1-(Vr^2)-(Vt^2))*(2*Vr + Vt*cot(theta)))/(((k-1)/2)*(1-(Vr^2)-(Vt^2)) - (Vt^2));

    X = [Vt;eq];
end


function [check, isterminal, direction] = boundaryCheck(~,x)
    %{
    The following function triggers a termination of the integrator when
    the Vtheta component reaches 0 (corresponding to having reached the
    surface of the cone).
    %}
    Vt = x(2);
    check = (Vt < 0); % event occurs when Vt drops below 0
    isterminal = 1; % terminate integration when event occurs
    direction = []; % direction does not matter (default)
    

end












