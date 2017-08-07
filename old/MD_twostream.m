%%%---------------------------------------------------------------------%%%
%
% TWO STREAM TESTCASE
% 
% MD test case for a two-stream instability problem. 
% Treat as two identical species, with one half of the particles traveling
% to the right at the thermal speed and the other half traveling to the
% left at the thermal speed.
%
% Everything is in dimensionless units based on the ion-circle radius,
% mass, charge, and plasma frequency of arbitrary reference parameters.
%
% The MD_Func function takes the species properties, distribution, and all
% necessary parameters to set up the 2D-2V MD domain.
%
% Species data:
%   species.Nsp     -   number of species
%   species.Np      -   total number of particles
%   species.mass    -   particle mass of each species species
%   species.charge  -   charge number for each species (Z*e = Q)
%
% Distribution data (for each species at each x-location):
%   distribution.f  -   distribution function with which to initialize
%   distribution.n  -   number density          (from 0th moment)
%   distribution.u  -   bulk velocity           (from 1st moment)
%   distribution.T  -   temperature (as energy) (from 2nd moment)
%
% Spatial and time domain:
%   params.Nx       -   number of x-direction gridpoints
%   params.Lx       -   length of domain in the x-direction
%   params.Nv       -   number of velocity space gridpoints
%   params.Lv       -   max velocity in velocity domain
%   params.dt       -   non-dimensional timestep to take (based on max(wp))
%   params.Ndt      -   number of timesteps to take
%   params.k        -   screening parameter, k = a/lambda
%
% Nearest neighbor:
%   params.cutoff   -   cutoff radius in terms of k
%   params.Ndiv     -   number of divisions per cutoff radius sized bin
%
% Equilibration:
%   params.gamma    -   equilibration parameter
%   params.Ndteq    -   number of timesteps to take for equilibration
%
% Data output/dumping:
%   params.fnames   -   filenames for the data dump
%   params.Njumps   -   number of timesteps between data dumps
%   params.Nsplits  -   number of files to split output into, requires:
%                       mod(Ndt,Nsplits) == 0, mod(Ndt/Nsplits,Njumps) == 0
%   params.timing   -   whether to output detailed timing data:
%                       0 - only output time when stages end
%                       1 - output bulk data at end
%                       2 - timing at every timestep
%
%%%--------------------------------------------------------------------%%%

% clear variables and set up paths
clear all%, close all
addpath('../BGK/')
addpath('../MD/')
addpath('../dependencies/')

% species properties (N, mass, charge are Nsp x 1 column vectors)
species.Nsp     =   2;                      % default: 1
species.Np      =   48000;                  % default: 5000
species.mass    =   [1 1];                      % default: 1
species.charge  =   [1 1];                      % default: 1
species.mass    =   reshape(species.mass,[species.Nsp],1);
species.charge  =   reshape(species.charge,[species.Nsp],1);

% user defined spatial and time parameters
params.Nx   =   32;                         % default: 10
params.Lx   =   300*pi;                     % default: 150
params.Nv   =   60;                         % default: 60
params.dt   =   0.00025;                       % default: 0.01
params.Ndt  =   120000;                       % default: 1000
params.k    =   1;                          % default: 1

% nearest neighbor list parameters
params.cutoff   =   5;                      % default: 5
params.Ndiv     =   1;                      % default: 1

% equlibration parameters
params.gamma    =   3;                      % default: 1
params.Ndteq    =   1000;                   % default: 1000

% data output parameters
params.fnames   =   {};                     % default: [] (use default)
params.Njumps   =   100;                    % default: 10
params.Nsplits  =   8;                      % default: 1
params.timing   =   2;                      % default: 2

% stuff to put in maxwellian for initial distribution function, based on
% the selected Gamma, and n = n0. (CHANGE TO SOMETHING INTERESTING)
Gamma   = 0.1;
navg    = 1/(pi);
T0      = 1/(2*Gamma);
u0      = sqrt(T0);
distribution.n  =   navg/2 * ones(species.Nsp,params.Nx);
distribution.T  =   T0 * ones(species.Nsp,params.Nx);
distribution.u  =   zeros(species.Nsp,params.Nx,2);
distribution.u(1,:) = u0;
distribution.u(2,:) = -u0;

% we can now get Lv and fill the maxwellian
params.Lv   =   10 * ...
    sqrt(max(max(distribution.T ./ repmat(species.mass,[1,params.Nx]))));
distribution.f = maxwellian([],species,distribution,params);

% set up filenames if unspecified
if isempty(params.fnames)
    params.fnames{1} = 'energy.dat';
    params.fnames{2} = 'movieMD.dat';
end

% call MD_Func
output = MD_Func(species,distribution,params);