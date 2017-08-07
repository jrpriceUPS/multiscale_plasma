%%%---------------------------------------------------------------------%%%
%
% MD demo case for a two-species plasma where each species is initialized
% with its own mass, charge, and initial conditions. 
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
%   params.dt       -   max non-dimensional timestep to take for measuring
%   params.accCoeff -   dt = min(dt, accCoeff/max(|acc|))
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
% Velocity Resample
%   params.refinement   -   refinement level when interpolating f (2^N)
%   params.Nmoments     -   number of moments to require good sample
%   params.abserr       -   absolute error allowed in resample (1xNmoments)
%   params.relerr       -   percent error allowed in resample (1xNmoments)
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
species.Nsp     =   2;                      % default: 2
species.Np      =   10000;                  % default: 10000
species.mass    =   [1 2];                  % default: [1 2]
species.charge  =   [1 2];                  % default: [1 2]

% make sure mass and charge vectors are correctly formatted
species.mass    =   reshape(species.mass,[species.Nsp],1);
species.charge  =   reshape(species.charge,[species.Nsp],1);

% user defined spatial and time parameters
params.Nx       =   10;                     % default: 10
params.Lx       =   100*pi;                 % default: 100*pi
params.Nv       =   60;                     % default: 60
params.dt       =   0.01;                   % default: 0.01
params.accCoeff =   0.5;                    % default: 0.5
params.Ndt      =   1000;                   % default: 1000
params.k        =   1;                      % default: 1

% nearest neighbor list parameters
params.cutoff   =   5;                      % default: 5
params.Ndiv     =   1;                      % default: 1

% equlibration parameters
params.gamma    =   1;                      % default: 1
params.Ndteq    =   500;                    % default: 500

% velocity resample parameters
params.refinement   =   4;                              % default: 4
params.Nmoments     =   4;                              % default: 4
params.abserr       =   [0.05 0.05 0.05 0.05];
params.relerr       =   [0.05 0.05 0.05 0.05];

% data output parameters
params.fnames   =   {};                     % default: {} (use default)
params.Njumps   =   10;                     % default: 10
params.Nsplits  =   1;                      % default: 1
params.timing   =   2;                      % default: 2

% stuff to put in maxwellian for initial distribution function, based on
% the criteria that k = 1
Gamma = [0.1; 0.1];
distribution.n  =   1/(pi*species.Nsp) * ones(species.Nsp,params.Nx);
distribution.T  =   repmat(species.charge.^2 ./ Gamma,[1,params.Nx]) .* ...
    sqrt(pi*distribution.n) / 2;
distribution.u  =   zeros(species.Nsp,params.Nx,2);

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