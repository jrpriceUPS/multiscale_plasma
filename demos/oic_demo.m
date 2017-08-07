%%%---------------------------------------------------------------------%%%
%
% Template for creating demos/testcases for the MD_Func function.
% Specifying the parameters allows for flexible testcase creationg.
% Default values are included in the template (as comments).
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
%   params.timing   -   whether to output detailed timing data:
%                       0 - only output time when stages end
%                       1 - output bulk data at end
%                       2 - timing at every timestep
%
%%%--------------------------------------------------------------------%%%
%clear all;close all;
k=12; %number of initial distributions

% clear variables and set up paths
%clear all%, close all
addpath('../BGK/')
addpath('../MD/')
addpath('../dependencies/')

plotparams.colors=distinguishable_colors(k);

% species properties
species.Np      =   4096;                   % default: 4096 (if grid, must be perfect square)
species.mass    =   1;                      % default: 1
species.charge  =   1;                      % default: 1

% user defined spatial and time parameters (Lv is derived)
params.Lx     =   (4*pi*species.Np/3)^(1/3);  % default: for well-defined a
params.Ly     =   params.Lx;                  % default: same as Lx (must be same for grid)
params.Lz     =   params.Lx;                  % default: same as Lx (must be same for grid)
params.dt     =   0.01;                       % default: 0.001
params.Ndt    =   500;                         % default: 1000
params.k      =   0.1;                          % default: 1
params.Nr     =   100;      % how many bins to use in radial distribution
params.record =   200;      % when to start recording   

% temperature
Gamma           =   1;
species.T       =   species.charge^2 / (3*Gamma);


%halton distribution, 2, 3, 5
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   3;
oic.h3        =   5;
outputhalton23 = order_induced_cooling(species,params,oic);

output(1).dt=outputhalton23.dt;
output(1).TList=outputhalton23.TList;
plotparams.legends{1}='Halton sequence, {2,3,5}';

%halton distribution, 2, 3, 7
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   3;
oic.h3        =   7;
outputhalton25 = order_induced_cooling(species,params,oic);

output(2).dt=outputhalton25.dt;
output(2).TList=outputhalton25.TList;
plotparams.legends{2}='Halton sequence, {2,3,7}';

%halton distribution, 2, 5, 7
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   7;
outputhalton35 = order_induced_cooling(species,params,oic);

output(3).dt=outputhalton35.dt;
output(3).TList=outputhalton35.TList;
plotparams.legends{3}='Halton sequence, {2,5,7}';

%halton distribution, 3, 5, 7
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   7;
outputhalton57 = order_induced_cooling(species,params,oic);

output(4).dt=outputhalton57.dt;
output(4).TList=outputhalton57.TList;
plotparams.legends{4}='Halton sequence, {3,5,7}';

%uniform distribution
oic.method    =   'uniform';
outputuniform = order_induced_cooling(species,params,oic);

%output.dt(5)=outputuniform.dt;
output(11).dt=outputuniform.dt;
output(11).TList=outputuniform.TList;
plotparams.legends{11}='Uniform';

%grid
oic.method    =   'grid';
outputgrid = order_induced_cooling(species,params,oic);

output(12).dt=outputgrid.dt;
output(12).TList=outputgrid.TList;
plotparams.legends{12}='Grid';

%halton distribution, 3, 5, 11
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   11;
outputhalton3511 = order_induced_cooling(species,params,oic);

output(5).dt=outputhalton3511.dt;
output(5).TList=outputhalton3511.TList;
plotparams.legends{5}='Halton sequence, {3,5,11}';

%halton distribution, 3, 7, 11
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   7;
oic.h3        =   11;
outputhalton3711 = order_induced_cooling(species,params,oic);

output(6).dt=outputhalton3711.dt;
output(6).TList=outputhalton3711.TList;
plotparams.legends{6}='Halton sequence, {3,7,11}';

%halton distribution, 2, 5, 11
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   11;
outputhalton2511 = order_induced_cooling(species,params,oic);

output(7).dt=outputhalton2511.dt;
output(7).TList=outputhalton2511.TList;
plotparams.legends{7}='Halton sequence, {2,5,11}';

%halton distribution, 3, 5, 13
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   13;
outputhalton3513 = order_induced_cooling(species,params,oic);

output(8).dt=outputhalton3513.dt;
output(8).TList=outputhalton3513.TList;
plotparams.legends{8}='Halton sequence, {3,5,13}';

%halton distribution, 5, 7, 13
oic.method    =   'halton';
oic.h1        =   5;
oic.h2        =   7;
oic.h3        =   13;
outputhalton5713 = order_induced_cooling(species,params,oic);

output(9).dt=outputhalton5713.dt;
output(9).TList=outputhalton5713.TList;
plotparams.legends{9}='Halton sequence, {5,7,13}';

%halton distribution, 2, 5, 13
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   13;
outputhalton2513 = order_induced_cooling(species,params,oic);

output(10).dt=outputhalton2513.dt;
output(10).TList=outputhalton2513.TList;
plotparams.legends{10}='Halton sequence, {2,5,13}';

