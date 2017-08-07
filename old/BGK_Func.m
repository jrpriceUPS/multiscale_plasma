function [f, n, u, T] = BGK_Func(species, f0, params)
%%%---------------------------------------------------------------------%%%
%
% Solves the 1D-2V BGK equation for n particle species based on an initial
% distribution function on a one-dimensional periodic domain with
% two-dimensional velocities.
%
% Returns the final distribution function plus its 0th through 2nd moments.
%
% INPUTS:
% 
% Species data (mass and charge should be Nsp x 1 column vectors):
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
%
% Plotting and timing parameters:
%   params.plotting -   whether of not to plot the results (boolean)
%   params.Njumps   -   if plotting, how often to plot data
%   params.timing   -   whether to output detailed timing info
%   params.legends  -   cell array of string names of species for plotting
%   params.axissize -   matrix of y-axis limits for density (row 1),
%                       temperature (row 2), and velocity (row 3)
%
% OUTPUTS:
%
%   f   -   the final distribution function after Ndt timesteps
%   n   -   final density at each spatial point                 [m^-3]
%   u   -   final bulk velocity at each spatial point           [m/s]
%   T   -   final temperature at each spatial point             [J]
%
%%%--------------------------------------------------------------------%%%
if nargin ~= 3
    error('Wrong number inputs. Must have species, f0, and parameters.')
end

A=[];
B=[];
relTol=1e-12;
absTol=1e-12;
T_e=[];

% begin timing
tic; start = toc;

% useful constants
EV_TO_J =   1.60217653e-19;     % convert eV to joules
EPS_0   =   8.854187817e-12;    % background permittivity       [F/m]
E_0     =   1.602176565e-19;    % elementary charge             [C]

% species data
Nsp     =   species.Nsp;
mass    =   species.mass;
charge  =   species.charge;

% spatial domain data
Nx  =   params.Nx;
Lx  =   params.Lx;
dx  =   Lx/Nx;
x   =   0.5*dx:dx:(Lx-0.5*dx);

% velocity space data
Nv  =   params.Nv;
Lv  =   params.Lv;
v   =   linspace(-Lv,Lv,Nv);
dv  =   v(2) - v(1);

% time domain data
dt      =   params.dt;
Ndt     =   params.Ndt;
Njumps  =   params.Njumps;

% other setup for plotting and timing
plotting    =   params.plotting;
timing      =   params.timing;
legends     =   params.legends;
axissize    =   params.axissize;
colors      =   distinguishable_colors(length(mass));
if plotting
    nfig    =   figure(1);  set(nfig,'OuterPosition',[50,550,500,500]);
    ufig    =   figure(2);  set(ufig,'OuterPosition',[550,550,500,500]);
    Tfig    =   figure(3);  set(Tfig,'OuterPosition',[550,50,500,500]);
end


% set integration weights
wtN         =   dv*ones(Nv,1);
wtN(1)      =   0.5*wtN(1);
wtN(end)    =   0.5*wtN(end);

% create repeated arrays from data for vectorized calculations and put in
% cell for maxwellians
W   =   permute(repmat(wtN*wtN.',[1,1,Nsp,Nx]),[3,4,1,2]);       % wi*wj
WV1 =   permute(repmat(wtN*(wtN.*v')',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*vi*wj
WV2 =   permute(repmat((wtN.*v')*wtN',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*wj*vj

v1      =   squeeze(repmat(v,[1,1,Nv]));
v2      =   v1';
v1R     =   permute(repmat(v1,[1,1,Nsp,Nx]),[3,4,1,2]);
v2R     =   permute(repmat(v2,[1,1,Nsp,Nx]),[3,4,1,2]);
v1adv   =   permute(repmat(v1,[1,1,Nx]),[3,1,2]);
massR   =   repmat(mass,[1,Nx]);
chargeR =   repmat(charge, [1,Nx]);

% for single species collision maxwellian
repeatedparams.v1R = v1R;
repeatedparams.v2R = v2R;
repeatedparams.massR = repmat(massR,[1,1,Nv,Nv]);
repeatedparams.nR = 0;
repeatedparams.uR = 0;
repeatedparams.TR = 0;
% for cross-species collision maxwellians
repeatedparamsX.v1R = v1R(1,:,:,:);
repeatedparamsX.v2R = v2R(1,:,:,:);
repeatedparamsX.massR = 0;
repeatedparamsX.nR = 0;
repeatedparamsX.uR = 0;
repeatedparamsX.TR = 0;

% initialize variables
f   =   f0;                         % distribution function
n   =   zeros(Nsp,Nx);              % density
u   =   zeros(Nsp,Nx,2);            % bulk velocity
T   =   zeros(Nsp,Nx);              % temperature

ne0 =   0;                          % domain electron density
F   =   zeros(Nsp,Nx);              % poisson forces
phi =   zeros(Nx,1);                % poisson potential
phinonlin = -0.55*ones(Nx,1);

fadv    =   zeros(size(f));         % advective component of df/dt
fcoll   =   zeros(size(f));         % collisional component of df/dt

% intialize timing variables
if timing
    avgmoments  =   0;
    avgadvect   =   0;
    avgpoisson  =   0;
    avgcollide  =   0;
end


%-------------------------------------------------------------------------%
% MAIN LOOP --------------------------------------------------------------%
%-------------------------------------------------------------------------%
for tstep = 1:Ndt
    
    % compute the moments
    getMoments;
    
    % if plotting, plot every Njumps
    if plotting && mod(tstep-1,Njumps) == 0
        plotMoments;
    end
    
    % get the field and forces using the poisson equation
    nonlinearPoisson;
    
    % compute the advection step using the force and velocities
    advectOne;
    
    % check for NaN and abort if detected
    if sum(isnan(fadv)) ~= 0
        error('NaN detected in advective component.')
    end
    
    % compute the collisional step
    BGKCollide;
    
    % update f
    f = f + fadv + fcoll;
    
end % end main loop

% get final moments and plot final distribution
tstep = tstep + 1;
getMoments;

if plotting
    plotMoments;
end

% timing output
fprintf('Total time elapsed: %f sec\n', toc - start)
if timing
    fprintf('Average moment calculation time: %f sec\n', avgmoments)
    fprintf('Average advection time: %f sec\n', avgadvect)
    fprintf('Average Poisson time: %f sec\n', avgpoisson)
    fprintf('Average collision time: %f sec\n', avgcollide)    
end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% NESTED FUNCTIONS -------------------------------------------------------%
%-------------------------------------------------------------------------%
    
%-------------------------------------------------------------------------%
% Moment Calculation %
    function getMoments
        if timing, start0 = toc; end
        
        % compute the density
        n = sum(sum(W.*f,4),3);

        % compute the velocity
        u(:,:,1) = sum(sum(WV1 .* f,4),3);
        u(:,:,2) = sum(sum(WV2 .* f,4),3);
        u = u ./ repmat(n,[1,1,2]);

        % compute the temperature
        u1R = repmat(u(:,:,1),[1,1,Nv,Nv]);
        u2R = repmat(u(:,:,2),[1,1,Nv,Nv]);
        vmu1 = v1R - u1R;
        vmu2 = v2R - u2R;
        T = sum(sum(W .* (vmu1.^2 + vmu2.^2) .*f,4),3) .* massR ./ (2*n);
        T(n == 0) = 0;
        
        if timing, avgmoments = avgmoments + (toc-start0) / (Ndt+1); end
    end % getMoments
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Linear Poisson Solver %
%
% Solves linearized version nonlinear 1D Poisson equation
%
% - phi_xx + g(phi) = s
%
% with periodic BCs
% where g is the nonlinear function
%
% n_{e0} exp(E_0 phi / T)
%
% and s is a source from the ion terms
%
% The discretization is
% A_ij phi_j + g(phi_i) = s_i
% where A is the usual three point stencil for Laplacian
%
    function nonlinearPoisson
        if timing, start0 = toc; end
        
        % set source
        source = sum(chargeR .* n, 1);
        source = source';
        
        % find ne0 and A if this is the first time through
        if tstep == 1
            ne0 = sum(source);
            
            % construct the second derivative discretization matrix
            A = (-diag(ones(Nx-1,1),-1) - diag(ones(Nx-1,1),1) + ...
                2*diag(ones(Nx,1))) / (dx*dx);
            
            % set periodic boundary conditions
            A(1,end) = A(1,2);
            A(end,1) = A(end,end-1);
            
            % set electron temperature (presently arbitrary)
            T_e = 200 * EV_TO_J;
            
            B=A-E_0^2*ne0/(T_e*EPS_0)*eye(Nx);
        end
        
        % compute phi
        phi=B\(ne0-source);
        
        %the force is the x derivative of phi, use a second order
        %discretization to find this
        F = (phi([2:Nx 1]) - phi([Nx 1:Nx-1])) / (2*dx);
        F = (4*pi * E_0^2 * chargeR ./ massR) .* repmat(F',[Nsp,1]);

        if timing, avgpoisson = avgpoisson + (toc - start0) / Ndt; end
    end % nonlinPoisson
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% First order in space advection scheme %
    function advectOne
        if timing, start0 = toc; end
        
        %initialize advection array
        fadv = zeros(size(f));
        
        %update advection using upwinding and downwinding
        for i=1:Nsp
            fadv(i,:,1:Nv/2,:) = -dt.*v1adv(:,1:Nv/2,:)/dx .*...
                squeeze(f(i,[2:Nx 1],1:Nv/2,:) - f(i,:,1:Nv/2,:));
            
            fadv(i,:,Nv/2+1:Nv,:) = -dt.*v1adv(:,Nv/2+1:Nv,:)/dx .*...
                squeeze(f(i,:,Nv/2+1:Nv,:) - f(i,[Nx 1:Nx-1],Nv/2+1:Nv,:));
        end
        
        %using the forces computed in the poisson subroutine, advance the
        %advection term using upstreaming and downstreaming
        Fpos = F > 0;
        Fneg = F < 0;
        Frep = repmat(F,[1,1,Nv,Nv]);
        for i=1:Nsp
            % F > 0
            fadv(i,Fpos(i,:),2:Nv,:) = fadv(i,Fpos(i,:),2:Nv,:) + ...
                dt/dv*Frep(i,Fpos(i,:),2:Nv,:) .* ...
                (f(i,Fpos(i,:),2:Nv,:) - f(i,Fpos(i,:),1:Nv-1,:));
            fadv(i,Fpos(i,:),1,:) = fadv(i,Fpos(i,:),1,:) + ...
                dt/dv*Frep(i,Fpos(i,:),1,:) .* f(i,Fpos(i,:),1,:);
            % F < 0
            fadv(i,Fneg(i,:),1:Nv-1,:) = fadv(i,Fneg(i,:),1:Nv-1,:) + ...
                dt/dv*Frep(i,Fneg(i,:),1:Nv-1,:) .* ...
                (f(i,Fneg(i,:),2:Nv,:) - f(i,Fneg(i,:),1:Nv-1,:));
            fadv(i,Fneg(i,:),Nv,:) = fadv(i,Fneg(i,:),Nv,:)...
                -dt/dv*Frep(i,Fneg(i,:),Nv,:) .* f(i,Fneg(i,:),Nv,:);
        end
        
        if timing, avgadvect = avgadvect + (toc - start0) / Ndt; end
    end % advectOne
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% BGK collisional subroutine (will be replaced by MD) %
    function BGKCollide
        if timing, start0 = toc; end
        
        % initialize colllision array
        fcoll = zeros(size(f));
        
        % compute bmax^2
        bmax2 = 1 ./ sum(chargeR.^2.*E_0.*E_0.*n ./ (EPS_0*T));
        
        % compute self-collision parameters
        bmin = chargeR.^2.*E_0.*E_0 ./ T;
        Lambda = repmat(bmax2,[Nsp,1]) ./ (bmin.*bmin);
        tau = 3*sqrt(2*massR) .* (2*pi*T).^1.5 .* EPS_0.^2 ./ ...
            ((E_0*chargeR).^4 .* 0.5.*log(1 + Lambda.^2));
        lambda_ii = n ./ tau;
        
        % compute self-collisions using relaxation to maxwellian
        repeatedparams.nR = repmat(n,[1,1,Nv,Nv]);
        repeatedparams.uR = permute(repmat(u,[1,1,1,Nv,Nv]),[1,2,4,5,3]);
        repeatedparams.TR = repmat(T,[1,1,Nv,Nv]);
        fcoll = fcoll + dt*repmat(lambda_ii,[1,1,Nv,Nv]) .* ...
            (maxwellian(repeatedparams) - f);
        
        
        % loop over all cross-species pairs for interspecies collisions
        for i = 1:Nsp
            for j = i+1:Nsp
                % reduced mass
                mu = mass(i) * mass(j) / (mass(i) + mass(j));
                
                % mixture temperature
                mixT = (n(i,:).*T(i,:) + n(j,:).*T(j,:)) ./ ...
                    (n(i,:) + n(j,:));
                
                % compute cross-collision parameters
                bmin = charge(i)*charge(j)*E_0*E_0 ./ mixT;
                Lambda = bmax2 ./ (bmin .* bmin);
                tau = 6*sqrt(2) * (pi*mixT./mu).^1.5 .* EPS_0*EPS_0 .* ...
                    mu.*mass(i) ./ ((charge(i)*charge(j)*E_0*E_0)^2 .* ...
                    0.5.*log(1 + Lambda.*Lambda));
                lambda_ij = n(j,:) ./ tau;
                lambda_ji = lambda_ij * mass(j)/mass(i) .* n(i,:)./n(j,:);
                
                % compute mixT_max and mixU_max
                nR = repmat(n,[1,1,2]);
                lambda_ijRu = repmat(lambda_ij,[1,1,2]);
                lambda_jiRu = repmat(lambda_ji,[1,1,2]);
                mixU_max = ...
                    (mass(i) * nR(i,:,:) .* lambda_ijRu .* u(i,:,:) + ...
                    mass(j) * nR(j,:,:) .* lambda_jiRu .* u(j,:,:)) ./ ...
                    (mass(i) * nR(i,:,:) .* lambda_ijRu + ...
                    mass(j) * nR(j,:,:) .* lambda_jiRu);
                mixT_max = (lambda_ij .* n(i,:) .* T(i,:) + ...
                    lambda_ji .* n(j,:) .* T(j,:) - ...
                    (mass(i) * n(i,:) .* lambda_ij .* ...
                    sum(mixU_max.^2 - u(i,:,:).^2,3) + ...
                    mass(j) * n(j,:) .* lambda_ji .* ...
                    sum(mixU_max.^2 - u(j,:,:).^2,3)) / 2) ./ ...
                    (n(i,:) .* lambda_ij + n(j,:).*lambda_ji);
                
                % get collisional contributions
                repeatedparamsX.uR = ...
                    permute(repmat(mixU_max,[1,1,1,Nv,Nv]),[1,2,4,5,3]);
                repeatedparamsX.TR = repmat(mixT_max,[1,1,Nv,Nv]);
                lambda_ijRv = repmat(lambda_ij,[1,1,Nv,Nv]);
                lambda_jiRv = repmat(lambda_ji,[1,1,Nv,Nv]);
                
                % contribution to species i
                repeatedparamsX.nR = repmat(n(i,:),[1,1,Nv,Nv]);
                repeatedparamsX.massR = mass(i);
                fcoll(i,:,:,:) = fcoll(i,:,:,:) + dt*lambda_ijRv .* ...
                    (maxwellian(repeatedparamsX) - f(i,:,:,:));
                
                % contribution to species j
                repeatedparamsX.nR = repmat(n(j,:),[1,1,Nv,Nv]);
                repeatedparamsX.mass = mass(j);
                fcoll(j,:,:,:) = fcoll(j,:,:,:) + dt*lambda_jiRv .* ...
                    (maxwellian(repeatedparamsX) - f(j,:,:,:));
            end
        end
        
        if timing, avgcollide = avgcollide + (toc-start0) / Ndt; end
    end % BGKCollide
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Make Plots %
    function plotMoments
        
        % density plot
        figure(1), clf, hold on
        for i = 1:Nsp
            plot(x,n(i,:),'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('n [m^{-2}]')
        if ~isnan(axissize)
            axis([0 Lx axissize(1,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Timestep %d', tstep-1))
        
        % x-velocity plot
        figure(2), clf, hold on
        for i = 1:Nsp
            plot(x,u(i,:,1),'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('u [m/s]')
        if isnan(axissize)==0
            axis([0 Lx axissize(3,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Timestep %d', tstep-1))
        
        % temperature plot
        figure(3), clf, hold on
        for i = 1:Nsp
            plot(x,T(i,:)/EV_TO_J,'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('T [eV]')
        if isnan(axissize)==0
            axis([0 Lx axissize(2,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Timestep %d', tstep-1))
        
        %pause(0.1)
        
    end % plotMoments
%-------------------------------------------------------------------------%

end % BGK_Func
