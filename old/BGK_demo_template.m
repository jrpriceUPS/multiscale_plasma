%%%--------------------------------------------------------------------%%%
%
% Template for creating demos/testcases for the BGK_Func function.
% Specifying domain paremeters allows for flexible testcase creation.
% Default values are included in the template (as comments)
%
% The BGK_Func function takes the mass, charge, initial distribution, and
% all necessary spatial and time domain parameters.
%
% Species data:
%   species.Nsp     -   number of species
%   species.mass    -   mass of each species                    [kg]
%   species.charge  -   charge number of each species
%
%   f0              -   initial distribution function (Nsp x Nx x Nv x Nv)
%
% Spatial and time domain:
%   params.Nx       -   number of sptatial points in discretization
%   params.Lx       -   width of the spatial domain             [m]
%   params.Nv       -   number of velocity points in discretization
%   params.Lv       -   maximum velocity to integrate           [m/s]
%   params.dt       -   timestep to use                         [s]
%   params.Ndt      -   number of timesteps to take
%   params.order    -   first-order (1) or second-order(2)
%   params.CFL      -   CFL number to choose time step
%
% Plotting and timing parameters:
%   params.plotting -   whether of not to plot the results (boolean)
%   params.Njumps   -   if plotting, how often to plot data
%   params.timing   -   whether to output detailed timing info
%   params.legends  -   cell array of string names of species for plotting
%   params.axissize -   matrix of y-axis limits for density (row 1),
%                       temperature (row 2), and velocity (row 3)
%
%%%--------------------------------------------------------------------%%%

% clear variables and set up paths
clear all, close all
addpath('../BGK/')
addpath('../MD/')
addpath('../dependencies/')

% constans
EV_TO_J =   1.60217653e-19;         % EV to J (for temperature conversion)
Mp      =   1.67262158e-27;         % Mass of a proton


% species properties (mass and charge should be column vectors (Nsp x 1))
species.Nsp     =   1;                          % default: 1
species.mass    =   [1]'*Mp;                    % default: 1*Mp
species.charge  =   [1]';                       % default: 1

% user specified spatial and time parameters (Lv is derived)
params.Nx   =   32;                             % default: 32
params.Lx   =   800e-9;                         % default: 800e-9
params.Nv   =   60;                             % default: 60
params.dt   =   1e-15;                          % default: 1e-15
params.Ndt  =   1000;                           % default: 1000

% plotting and timing parameters
params.plotting =   0;                          % default: 0
params.Njumps   =   10;                         % default: 10
params.timing   =   0;                          % default: 0
params.legends  =   {};                         % default: {} (unspecified)
params.axissize =   [];                         % default: [] (unspecified)


%-------------------------------------------------------------------------%
% stuff to put in maxwellian for initial distribution function, f0
% template has boring, uniform dummy values: CHANGE THIS PART
moments.n   =   ones(species.Nsp,params.Nx);                % density
moments.u   =   zeros(species.Nsp,params.Nx,2);             % velocity
moments.T   =   100*ones(species.Nsp,params.Nx)*EV_TO_J;    % temperature
%-------------------------------------------------------------------------%

% we can now get Lv and fill the maxwellian
params.Lv   =   10 * ...
    sqrt(max(max(moments.T ./ repmat(species.mass,[1,params.Nx]))));
f0 = maxwellian([],species,moments,params);

% set up axis and legend if not specified
if isempty(params.legends)
    legends = cell(species.Nsp,1);
    for i = 1:species.Nsp
        params.legends{i} = sprintf('Species %d', i);
    end
end
if isempty(params.axissize)
    params.axissize = NaN;
end

% call BGK_func
BGK_Func(species, f0, params);