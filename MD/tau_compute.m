function tau=tau_compute(species,distribution,timeseries,params)
%%%---------------------------------------------------------------------%%%
%
% Given time series data from MD and distribution data from BGK, compute
% the tau such that the average rate of change of the entropy at each
% spatial position for each species matches between models.
%
%
% INPUTS:
%
% Species data:
%   species.Nsp       -   number of species
%   species.Np        -   total number of particles
%   species.mass      -   particle mass of each species
%   species.charge    -   charge number for each species (Z*e = Q)
%
% Distribution data (for each species at each x-location):
%   distribution.f    -   distribution function from BGK
%   distribution.n    -   number density          (from 0th moment)
%   distribution.u    -   bulk velocity           (from 1st moment)
%   distribution.T    -   temperature (as energy) (from 2nd moment)
%
% Time series data (for each species at each x-location):
%   timeseries.tList  -   list of time of measurement
%   timeseries.HList  -   list of entropies at time
%
% OUTPUT:
%   tau               -   the tau parameters needed to match the two models
%
%%% -----------------------------------------------------------------------

% particles and species
Nsp     =   species.Nsp;

%distribution
f       =   distribution.f;

%time series
tList   =   timeseries.tList;
HList   =   timeseries.HList;

% spatial and velocity parameters
Nx      =   params.Nx;
Nv      =   params.Nv;
Lv      =   params.Lv;

% system
vlist   = linspace(-Lv,Lv,Nv);
dv      = vlist(2)-vlist(1);


% set integration weights
wtN         =   dv*ones(Nv,1);
wtN(1)      =   0.5*wtN(1);
wtN(end)    =   0.5*wtN(end);
W           =   permute(repmat(wtN*wtN.',[1,1,Nsp,Nx]),[3,4,1,2]);  % wi*wj

% set the temperature to the arithmetic average of the temperature of each
% species in each grid cell

distribution.T = repmat(sum(distribution.T.*distribution.n)./sum(distribution.n),Nsp,1);


% set equilibrium distribution
f_eq    =  maxwellian([],species,distribution,params);


% set locations with zero probability to 1 so log(f_k) does not blow up
f(f==0) =  1;


% compute entropy production rate in BGK sense
BGK_H =  sum(sum(W.*((f_eq-f).*log(f)),4),3);


% for each species and position, fit the rate of change of entropy in MD
% using linear regression
MD_H  =  zeros(Nsp,Nx);
for sp=1:Nsp
    for ix=1:Nx
        P = polyfit(tList.',squeeze(HList(sp,ix,:)),1);
        MD_H(sp,ix) = P(1);
    end
end

% compute tau
tau = BGK_H./MD_H;

% if there are no particles of species k at point x, set tau to zero
tau(distribution.n==0) = 0;