clear; format short e;
% See Shaughnessy's fluids textbook p. 841 for example.
% All conservative assumptions are made to minimize head loss, which
% maximizes flow rate, which maximizes the expected pressure spike.

% Point 1 is the surface of the water in the tank.
% Point 2 is the water flowing out the end of the tubing.


%% Calculate velocity

% Parameters
rho = 1000; %kg/m^3
g = 9.8; %m/s^2
mu = 1.002e-3; %Ns/m^2

% Known states
p1 = (100 + 14.7) * 6894.76; %psi atmospheric to Pa
%p2 = 0 * 6894.76; %psi atmospheric to Pa
p2=14.7*6894.76; %psi atmospheric to Pa (ambient; not vacuum)
z1 = 0.25; %m %use the level of the tubing as datum, assume water level in tank is 0.25 m above tubing
z2 = 0; %m
V1 = 0; %m/s %assume velocity of water surface is negligible

% Head loss sources
%D = 1/4 * .0254; %in to m ***Check actual ID of 1/4" tubing***
D=0.18*0.0254; %http://www.mcmaster.com/#standard-metal-tubing/=10zskyw
A = pi/4 * D^2; %m^2
L = 60 * .0254; %in to m %total length of pipe
%e = 0; %for smooth tubes ***check what it is for the tubing we're using***
e=0.045/1000; % (mm to m) Commercial steel (pg.815 table 13.3)
K_in = .04; % Fig 13.20, case d (conservative case)

n_handvalve = 1;
n_fastvalve = 1;
n_elbow = 0;
n_tee = 1;
n_cross = 1;

LeD_valve = 3; % Table 13.6, ball valve, fully open
LeD_elbow = 30; % Table 13.6, elbow, 90 degree standard
LeD_tee = 20; % Table 13.6, standard tee, flow through run
LeD_cross = 20; % ***Temporarily assuming cross is the same as tee***

Le_tot = D * (LeD_valve*n_handvalve + LeD_valve*n_fastvalve + ...
    LeD_elbow*n_elbow + LeD_tee*n_tee + LeD_cross*n_cross);

% Estimate frictionless exit velocity
H = z1-z2; %m %height difference between two points
dP = p1-p2; %Pa %pressure difference between two points
V2a = sqrt(2* (g*H + dP/rho)) %m/s %frictionless exit velocity

% Iterations
Re = rho * V2a * D / mu;

if Re < 2300 %laminar flow
    fprintf('Laminar')
    if L/D > 0.06 * Re; % from p.794
        fprintf('Fully developed');
    else
        fprintf('Not fully developed');
    end
    alpha = 2; %See page 823 (a=1 for turbulent; a=2 for laminar)
    f = 64/Re; %Eqn 13.19 
    V2b = sqrt(2* (g*H + dP/rho) / (alpha + f/D * (L+Le_tot) + K_in))
    Q=V2b*A
end

if Re > 2300 %transition or turbulent flow
    fprintf('Turbulent')
    if L/D > 4.4 * Re ^(1/6); % from p.794
        fprintf('Fully developed');
    else
        fprintf('Not fully developed');
    end
    alpha = 1;
    f = (-2 * log10(e/D / 3.7065 - 5.0452/Re * log((e/D)^1.1098 / 2.8257 + ...
        5.8506 / Re^0.8981))) ^(-2); %Chen Eqn 3.19b
    V2b = sqrt(2* (g*H + dP/rho) / (alpha + f/D * (L+Le_tot) + K_in))
    Q=V2b*A %m^3/s
end


%% Calculate pressure spike without device
c = 1482; %m/s
p_wh = rho*V2b*c * 1.45038e-4 %Pa to psi


