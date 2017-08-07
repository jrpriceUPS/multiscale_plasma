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
params.Ndt    =   0;                         % default: 1000
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
output1 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {2,3,5}';

figure
plot(output1.rRange,output1.g);
title(legend);


%halton distribution, 2, 3, 7
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   3;
oic.h3        =   7;
output2 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {2,3,7}';

figure
plot(output2.rRange,output2.g);
title(legend);

%halton distribution, 2, 5, 7
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   7;
output3 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {2,5,7}';

figure
plot(output3.rRange,output3.g);
title(legend);

%halton distribution, 3, 5, 7
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   7;
output4 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {3,5,7}';

figure
plot(output4.rRange,output4.g);
title(legend);

%uniform distribution
oic.method    =   'uniform';
output5 = order_induced_cooling(species,params,oic);

legend='Uniform';

figure
plot(output5.rRange,output5.g);
title(legend);

%grid
oic.method    =   'grid';
output6 = order_induced_cooling(species,params,oic);

legend='Grid';

figure
plot(output6.rRange,output6.g);
title(legend);

%halton distribution, 3, 5, 11
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   11;
output7 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {3,5,11}';

figure
plot(output7.rRange,output7.g);
title(legend);

%halton distribution, 3, 7, 11
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   7;
oic.h3        =   11;
output8 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {3,7,11}';

figure
plot(output8.rRange,output8.g);
title(legend);

%halton distribution, 2, 5, 11
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   11;
output9 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {2,5,11}';

figure
plot(output9.rRange,output9.g);
title(legend);

%halton distribution, 3, 5, 13
oic.method    =   'halton';
oic.h1        =   3;
oic.h2        =   5;
oic.h3        =   13;
output10 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {3,5,13}';

figure
plot(output10.rRange,output10.g);
title(legend);

%halton distribution, 5, 7, 13
oic.method    =   'halton';
oic.h1        =   5;
oic.h2        =   7;
oic.h3        =   13;
output11 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {5,7,13}';

figure
plot(output11.rRange,output11.g);
title(legend);

%halton distribution, 2, 5, 13
oic.method    =   'halton';
oic.h1        =   2;
oic.h2        =   5;
oic.h3        =   13;
output12 = order_induced_cooling(species,params,oic);

legend='Halton sequence, {2,5,13}';

figure
plot(output12.rRange,output12.g);
title(legend);

