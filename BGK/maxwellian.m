function fm = maxwellian(repeatedparams, species, moments, params)
%%%---------------------------------------------------------------------%%%
%
% Creates a maxwellian distribution array of size Nsp x Nx x Nv x Nv
% using data about the species, known moments, and domain. This is designed
% for consistency and compatibility with the variables used in BGK_func.
%
% Everything is in dimensionless units based on the ion-circle radius,
% mass, charge, and plasma frequency of arbitrary reference parameters.
% 
% Since the equation doesn't actually change, this also will work with
% dimensional units.
%
% INPUTS:
%
% Species data:
%   species.Nsp     -   number of species
%   species.mass    -   mass of each species
%   species.charge  -   charge number of each species
%
% Moments:
%   moments.n       -   density                 (Nsp x Nx)
%   moments.u       -   bulk velocity           (Nsp x Nx x 2)
%   moments.T       -   temperature (energy)    (Nsp x Nx)
%
% Parameters:
%   params.Nx       -   number of sptatial points in discretization
%   params.Lx       -   width of the spatial domain
%   params.Nv       -   number of velocity points in discretization
%   params.Lv       -   maximum velocity to integrate
%
% Repeated Parameters (to save time)
%   v1R             -   velocity repeated to be Nsp x Nx x Nv x Nv
%   v2R             -   "transpose" of v1R
%   massR           -   mass repeated to be Nsp x Nx x Nv x Nv
%   nR              -   density repeated to be Nsp x Nx x Nv x Nv
%   uR              -   velocity repeated to be Nsp x Nx x Nv x Nv x 2
%   TR              -   temperature repeated to be Nsp x Nx x Nv x Nv
%
% OUTPUT:
% fm    -   the maxwellian with the specified moments
%%%---------------------------------------------------------------------%%%
if nargin > 1
    Nsp = species.Nsp;
    Nx  = params.Nx;
    Nv  = params.Nv;
    Lv  = params.Lv;

    % set up repeated arrays for vectorized operation
    v       =   linspace(-Lv,Lv,Nv);
    v1      =   squeeze(repmat(v,[1,1,Nv]));
    v2      =   v1';
    v1R     =   permute(repmat(v1,[1,1,Nsp,Nx]),[3,4,1,2]);
    v2R     =   permute(repmat(v2,[1,1,Nsp,Nx]),[3,4,1,2]);
    massR   =   repmat(species.mass,[1,Nx,Nv,Nv]);
    nR      =   repmat(moments.n,[1,1,Nv,Nv]);
    uR      =   permute(repmat(moments.u,[1,1,1,Nv,Nv]),[1,2,4,5,3]);
    TR      =   repmat(moments.T,[1,1,Nv,Nv]);
else
    v1R     =   repeatedparams.v1R;
    v2R     =   repeatedparams.v2R;
    massR   =   repeatedparams.massR;
    nR      =   repeatedparams.nR;
    uR      =   repeatedparams.uR;
    TR      =   repeatedparams.TR;
end

% construct the maxwellian
fm  =   massR .* nR ./ (2*pi*TR) .* exp(-0.5 * massR .* ...
    ((v1R - uR(:,:,:,:,1)).^2 + (v2R - uR(:,:,:,:,2)).^2) ./ TR);

% for regions with zero density, set the distribution to zero
fm(nR==0) = 0;