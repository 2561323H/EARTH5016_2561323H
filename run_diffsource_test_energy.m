%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

%NN = [100,200,400];


%for nn = 1:3 


% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target no. of columns in x-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image
test = 'no';


% units = value of each pixel (colour)
% D = original depth
% Nz = target no. of rows in z-direction

[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit taken from provided excel file and
% DOI: 10.1023/B:NARR.0000032647.41046.e7.


matprop = [
        % unit  conductivity  density  heat capacity  heat production      % porosity
          1	    3.6788	        2697.6	    1172 	    4.172e-6              1.03            % Granite phase 1, 4.172alternative Cp 1172
          2	    2.465           2700	    979	        2.9e-6                1.2          % Basement gneiss, alternative Cp 979
          3	    3.2197	        2703.5	    1172	    5.575e-6              1.03         % Granite phase 2
          4	    0.77	        1942.3	    740	        0.75e-6               25.3             % Sand
          5	    0.77	        2648	    740	        0.95e-6               32.5          % Gravel
          6	    0.924	        2081.7	    860	        1.43e-6               0.55           % Clay, mudstone
          7	    1.67	        1916	    910	        1.3e-6                17         % Silt
          8	    0.919	        1909.78	    740	        1.2e-6                21.2        % Mud, silt, sand
          9	    1e-6            1   	    1000	    0                     100];           % air/water
        

switch test 
    
    case 'no' 

        
        % get coefficient fields based on spatial distribution of rock units from image
        
        rho    = reshape(matprop(units,3),Nz,Nx); % density
        Cp     = reshape(matprop(units,4),Nz,Nx); % specific heat capacity
        kT     = reshape(matprop(units,2),Nz,Nx); % conductivity
        Hr     = reshape(matprop(units,5),Nz,Nx); % heat rate
        
        

        % calculate heat diffusivity [m2/s]
        k0 = kT ./ rho ./ Cp;

    case 'yes'

        rho    = 2400*ones(Nz,Nx); % density
        Cp     = 1000*ones(Nz,Nx); % specific heat capacity
        kT     = ones(Nz,Nx); % conductivity
        Hr     = ones(Nz,Nx); % heat rate

        % calculate heat diffusivity [m2/s]
        k0 = kT ./ rho ./ Cp;
        

end

% set model parameters
dTdz = [0, 35/1000];  % set boundary condition
T0  = 10;              % surface temperature degree C
Tair = 10;             % air temperature degree C
nop   = 5000;          % output figure produced every 'nop' steps
wT   = 10;         % initial temperature peak width [m]
yr    = 3600*24*365;  % seconds per year [s]
tend  = 1e6 *yr;       % stopping time [s]
CFL   = 1/5;         % Time step limiter

%*****  RUN MODEL
run('./diffsource_test_energy.m');

