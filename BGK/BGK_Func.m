function distribution = BGK_Func(species, f0, params)
%%%---------------------------------------------------------------------%%%
%
% Solves the 1D-2V BGK equation for n particle species based on an initial
% distribution function on a one-dimensional periodic domain with
% two-dimensional velocities.
%
% Everything is in dimensionless units based on the ion-circle radius,
% mass, charge, and plasma frequency.
%
% Returns the final distribution function plus its 0th through 2nd moments.
%
% Uses strang operator splitting and RK2 for the time evolution.
%
% INPUTS:
% 
% Species data (mass and charge should be Nsp x 1 column vectors):
%   species.Nsp     -   number of species
%   species.mass    -   mass of each species
%   species.charge  -   charge number of each species
%
%   f0              -   initial distribution function (Nsp x Nx x Nv x Nv)
%
% Spatial and time domain:
%   params.Nx           -   number of sptatial points in discretization
%   params.Lx           -   length of the spatial domain
%   params.Nv           -   number of velocity points in discretization
%   params.Lv           -   maximum velocity to integrate
%   params.Ndt          -   number of timesteps to take
%   params.order        -   first-order (1) or second-order (2)
%   params.minmodtheta  -   parameter for minmod (second order only)
%   params.CFL          -   CFL number to choose time step
%   params.k            -   screening parameter, k = a/lambda
%
% Plotting and timing parameters:
%   params.plotting -   whether of not to plot the results (boolean)
%   params.Njumps   -   if plotting, how often to plot data
%   params.saveplot -   whether to save the plots to image files
%   params.timing   -   whether to output detailed timing info
%   params.legends  -   cell array of string names of species for plotting
%   params.axissize -   matrix of y-axis limits for density (row 1),
%                       temperature (row 2), and velocity (row 3)
%
% OUTPUTS:
%
%   distribution.f  -   the final distribution function after Ndt timesteps
%   distribution.n  -   final density at each spatial point
%   distribution.u  -   final bulk velocity at each spatial point
%   distribution.T  -   final temperature at each spatial point
%
%%%--------------------------------------------------------------------%%%
if nargin ~= 3
    error('Wrong number inputs. Must have species, f0, and parameters.')
end

% begin timing
tic; start = toc;

% species data
Nsp     =   species.Nsp;
mass    =   species.mass;
charge  =   species.charge;

% spatial domain data
Nx  =   params.Nx;
Lx  =   params.Lx;
dx  =   Lx/Nx;
x   =   0.5*dx:dx:(Lx-0.5*dx);
k   =   params.k;

% velocity space data
Nv  =   params.Nv;
Lv  =   params.Lv;
v   =   linspace(-Lv,Lv,Nv);
dv  =   v(2) - v(1);

% time domain data
CFL     =   params.CFL;
Ndt     =   params.Ndt;
Njumps  =   params.Njumps;

% order in time and space
order       = params.order;
minmodtheta = params.minmodtheta;
if order ~= 1 && order ~= 2
    error('Must specify either first or second order evolution.')
end

% other setup for plotting and timing
plotting    =   params.plotting;
saveplot    =   params.saveplot;
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

% create repeated arrays of integration weights for computing double
% integrals
W   =   permute(repmat(wtN*wtN.',[1,1,Nsp,Nx]),[3,4,1,2]);       % wi*wj
WV1 =   permute(repmat(wtN*(wtN.*v')',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*vi*wj
WV2 =   permute(repmat((wtN.*v')*wtN',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*wj*vj

% create repeated matrices for computing Maxwellians and componentwise
% multiplication
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

if Nsp>1
    % for cross-species collision maxwellians
    repeatedparamsX.v1R = v1R(1,:,:,:);
    repeatedparamsX.v2R = v2R(1,:,:,:);
    repeatedparamsX.massR = 0;
    repeatedparamsX.nR = 0;
    repeatedparamsX.uR = 0;
    repeatedparamsX.TR = 0;
end

% initialize variables
t   =   0;                          % time

f   =   f0;                         % distribution function
n   =   zeros(Nsp,Nx);              % density
u   =   zeros(Nsp,Nx,2);            % bulk velocity
T   =   zeros(Nsp,Nx);              % temperature

F   =   zeros(Nsp,Nx);              % poisson forces

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
% if plotting, plot initial distribution
tstep = 0;
if plotting
        getMoments;
        plotMoments;
end
    
% get first timestep using CFL criteria
getMoments;
ElectricField;
wp = sqrt(chargeR.^2 .* (pi*n).^(3/2) ./ massR);
dt = CFL * min([dx/Lv dv/max(abs(F(:))) 1/max(wp(:))]);

% if second order, we need to do an initial velocity-advection half-step
if order == 2
    velocityAdvectTwo(dt/4);
    fold = f;
    f = f + fadv;
    velocityAdvectTwo(dt/2);
    f = fold + fadv;
end

for tstep = 1:Ndt
    % if plotting, plot every Njumps
    if plotting && mod(tstep,Njumps) == 0
        getMoments;
        plotMoments;
    end
    
    % get timestep using CFL criteria
    if tstep > 1 % already covered for first step
        getMoments;
        ElectricField;
        wp = sqrt(chargeR.^2 .* (pi*n).^(3/2) ./ massR);
        dt = CFL * min([dx/Lv dv/max(abs(F(:))) 1/max(wp(:))]);
    end
    
    if order == 1
        % compute the velocity-space advection step using the forces
        velocityAdvectOne(dt);
        f = f + fadv;
        
        % compute the space advection step using the velocities
        spaceAdvectOne(dt);
        f = f + fadv;
        
        % compute the collisional step
        BGKCollide(dt);
        f = f + fcoll;
        
    elseif order == 2
        % compute half-step in BGK term
        BGKCollide(dt/4);
        fold = f;
        f = f + fcoll;
        BGKCollide(dt/2);
        f = fold + fcoll;
        
        % compute full step in space advection
        spaceAdvectTwo(dt/2);
        fold = f;
        f = f + fadv;
        spaceAdvectTwo(dt);
        f = fold + fadv;
        
        % compute half-step in BGK term
        BGKCollide(dt/4);
        fold = f;
        f = f + fcoll;
        BGKCollide(dt/2);
        f = fold + fcoll;
        
        % compute full step in velocity advection (except last step)
        if tstep < Ndt
            velocityAdvectTwo(dt/2);
            fold = f;
            f = f + fadv;
            velocityAdvectTwo(dt);
            f = fold + fadv;
        end
    end
       
    t = t + dt;
end % end main loop

% last half-step in velocity advection (second order)
if order == 2
    velocityAdvectTwo(dt/4);
    fold = f;
    f = f + fadv;
    velocityAdvectTwo(dt/2);
    f = fold + fadv;
end

% get final moments and plot final distribution
getMoments;
if plotting
    plotMoments;
end

distribution.f = f;
distribution.n = n;
distribution.u = u;
distribution.T = T;

% timing output
fprintf('Total time elapsed: %f sec\n', toc - start)
if timing
    fprintf('Average moment calculation time: %f sec\n', avgmoments/order)
    fprintf('Average advection time: %f sec\n', avgadvect/order)
    fprintf('Average Poisson time: %f sec\n', avgpoisson/order)
    fprintf('Average collision time: %f sec\n', avgcollide/order)    
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
% Electric Field Solver %
%
% Finds the electric field given the charge distribution by solving the
% screened Poisson equation
%
%d^2(phi)/dx^2 - 1/lambda^2 phi = -1/2 rho
%
    function ElectricField
        if timing, start0 = toc; end
        
        % Charge density of ions
        source  =  sum(chargeR.*n,1);
        
        % Standard discretization of Poisson operator
        A = 1/(dx^2)*(-diag(ones(Nx-1,1),-1) - diag(ones(Nx-1,1),1) ...
            + 2*diag(ones(Nx,1))) ;
        
        % Periodic BCs
        A(1,end) = A(1,2);
        A(end,1) = A(end,end-1);
        
        B = A - k^2*eye(Nx);
        
        
        % the force is the x derivative of phi, use a second order
        % discretization to find this
        phi =  B/(source/(2/pi));
        
        E   =  (phi([2:Nx 1])-phi([Nx 1:Nx-1]))/(2*dx);
        F   =  -chargeR ./ massR .* repmat(E',[Nsp,1]);
        
        if timing, avgpoisson = avgpoisson + (toc - start0) / Ndt; end
    end % ElectricField
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% First order in space advection scheme % 
% solves dt/dt + vx*df/dx = 0
    function spaceAdvectOne(dt)
        if timing, start0 = toc; end
        
        % initialize advection array
        fadv = zeros(size(f));
        
        % update advection using upwinding and downwinding
        for i = 1:Nsp
            fadv(i,:,1:Nv/2,:) = -dt.*v1adv(:,1:Nv/2,:)/dx .*...
                squeeze(f(i,[2:Nx 1],1:Nv/2,:) - f(i,:,1:Nv/2,:));
            
            fadv(i,:,Nv/2+1:Nv,:) = -dt.*v1adv(:,Nv/2+1:Nv,:)/dx .*...
                squeeze(f(i,:,Nv/2+1:Nv,:) - f(i,[Nx 1:Nx-1],Nv/2+1:Nv,:));
        end
        
        errorCheck;
        
        if timing, avgadvect = avgadvect + (toc - start0) / Ndt; end
    end % spaceAdvectOne
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% First order in velocity space advection scheme %
% solves dt/dt + F*df/dvx = 0
    function velocityAdvectOne(dt)
        if timing, start0 = toc; end
        
        % get the moments and the Poisson force terms
        getMoments;
        ElectricField;
        
        % initialize advection array
        fadv = zeros(size(f));
        
        % using the forces computed in the poisson subroutine, advance the
        % advection term using upstreaming and downstreaming
        Fpos = F > 0;
        Fneg = F < 0;
        Frep = repmat(F,[1,1,Nv,Nv]);
        
        for i = 1:Nsp
            % F > 0
            fadv(i,Fpos(i,:),2:Nv,:) = ...
                dt/dv*Frep(i,Fpos(i,:),2:Nv,:) .* ...
                (f(i,Fpos(i,:),2:Nv,:) - f(i,Fpos(i,:),1:Nv-1,:));
            fadv(i,Fpos(i,:),1,:) = dt/dv*Frep(i,Fpos(i,:),1,:) .* ...
                fadv(i,Fpos(i,:),1,:);
            % F < 0
            fadv(i,Fneg(i,:),1:Nv-1,:) = ...
                dt/dv*Frep(i,Fneg(i,:),1:Nv-1,:) .* ...
                (f(i,Fneg(i,:),2:Nv,:) - f(i,Fneg(i,:),1:Nv-1,:));
            fadv(i,Fneg(i,:),Nv,:) = -dt/dv*Frep(i,Fneg(i,:),Nv,:) .*... 
                fadv(i,Fneg(i,:),Nv,:);
        end
        
        errorCheck;
        
        if timing, avgadvect = avgadvect + (toc - start0) / Ndt; end
    end % velocityAdvectOne
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Second order in space advection scheme %
% solves df/dt + vx*df/dx = 0
    function spaceAdvectTwo(dt)
        if timing, start0 = toc; end
        
        % initialize advection array
        fadv = zeros(size(f));
        
        % calculate all x-slopes using minmod (periodic)
        xslopes = zeros(size(f));
        minmod  = zeros([size(f) 3]);
        minmod(:,:,:,:,1) = f - f(:,[Nx 1:Nx-1],:,:);
        minmod(:,:,:,:,2) = f(:,[2:Nx 1],:,:) - f;
        minmod(:,:,:,:,3) = (f(:,[2:Nx 1],:,:) - f(:,[Nx 1:Nx-1],:,:)) ...
            / (2*minmodtheta);
        mmodmin = min(minmod,[],5);
        mmodmax = max(minmod,[],5);
        indicator = sum(sign(minmod),5);
        xslopes(indicator == 3)  = mmodmin(indicator == 3) * minmodtheta;
        xslopes(indicator == -3) = mmodmax(indicator == -3) * minmodtheta;
        
        % update advection using upwinding and downwinding
        for i = 1:Nsp
            % v < 0
            fadv(i,:,1:Nv/2,:) = -dt.*v1adv(:,1:Nv/2,:)/dx .* ...
                squeeze(f(i,[2:Nx 1],1:Nv/2,:) - ...
                0.5*xslopes(i,[2:Nx 1],1:Nv/2,:) - ...
                f(i,:,1:Nv/2,:) + ...
                0.5*xslopes(i,:,1:Nv/2,:));
            % v > 0
            fadv(i,:,Nv/2+1:Nv,:) = -dt.*v1adv(:,Nv/2+1:Nv,:)/dx .* ...
                squeeze(f(i,:,Nv/2+1:Nv,:) + ...
                0.5*xslopes(i,:,Nv/2+1:Nv,:) - ...
                f(i,[Nx 1:Nx-1],Nv/2+1:Nv,:) - ...
                0.5*xslopes(i,[Nx 1:Nx-1],Nv/2+1:Nv,:));
        end
        
        errorCheck;
        
        if timing, avgadvect = avgadvect + (toc - start0) / Ndt; end
    end % spaceAdvectTwo
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Second order in velocity-space advection scheme %
% solves df/dt + F*df/dvx = 0
    function velocityAdvectTwo(dt)
        if timing, start0 = toc; end
        
        % get moments and calculate poisson term
        getMoments;
        ElectricField;
        
        % initialize advection array
        fadv = zeros(size(f));
        
        % calculate all vx-slopes using minmod (f = 0 if j < 0 or j > Nv)
        vslopes = zeros(size(f));
        minmod  = zeros([size(f) 3]);
        minmod(:,:,2:Nv,:,1)    = f(:,:,2:Nv,:) - f(:,:,1:Nv-1,:);
        minmod(:,:,1,:,1)       = f(:,:,1,:);
        minmod(:,:,1:Nv-1,:,2)  = f(:,:,2:Nv,:) - f(:,:,1:Nv-1,:);
        minmod(:,:,Nv,:,2)      = -f(:,:,Nv,:);
        minmod(:,:,2:Nv-1,:,3)  = (f(:,:,3:Nv,:) - f(:,:,1:Nv-2,:)) / 2;
        minmod(:,:,1,:,3)       = f(:,:,2,:) / 2;
        minmod(:,:,Nv,:,3)      = -f(:,:,Nv-1,:) / 2;
        minmodmin = min(minmod,[],5);
        minmodmax = max(minmod,[],5);
        indicator = sum(sign(minmod),5);
        vslopes(indicator == 3)  = minmodmin(indicator == 3);
        vslopes(indicator == -3) = minmodmax(indicator == -3);
        
        % using the forces computed in the poisson subroutine, advance the
        % advection term using upstreaming and downstreaming
        Fpos = F > 0;
        Fneg = F < 0;
        Frep = repmat(F,[1,1,Nv,Nv]);
        
        for i = 1:Nsp
            % F > 0
            fadv(i,Fpos(i,:),2:Nv,:) = ...
                dt/dv*Frep(i,Fpos(i,:),2:Nv,:) .* ...
                (f(i,Fpos(i,:),2:Nv,:) + ...
                0.5*vslopes(i,Fpos(i,:),2:Nv,:) - ...
                f(i,Fpos(i,:),1:Nv-1,:) - ...
                0.5*vslopes(i,Fpos(i,:),2:Nv,:));
            fadv(i,Fpos(i,:),1,:) = ...
                dt/dv*Frep(i,Fpos(i,:),1,:) .* ...
                (f(i,Fpos(i,:),1,:) + 0.5*vslopes(i,Fpos(i,:),1,:));
            % F < 0
            fadv(i,Fneg(i,:),1:Nv-1,:) = ...
                dt/dv*Frep(i,Fneg(i,:),1:Nv-1,:) .* ...
                (f(i,Fneg(i,:),2:Nv,:) - ...
                0.5*vslopes(i,Fneg(i,:),2:Nv,:) - ...
                f(i,Fneg(i,:),1:Nv-1,:) + ...
                0.5*vslopes(i,Fneg(i,:),1:Nv-1,:));
            fadv(i,Fneg(i,:),Nv,:) = ...
                dt/dv*Frep(i,Fneg(i,:),Nv,:) .* ...
                (f(i,Fneg(i,:),Nv,:) - 0.5*vslopes(i,Fneg(i,:),Nv,:));
        end
        
        errorCheck;
        
        if timing, avgadvect = avgadvect + (toc - start0) / Ndt; end
    end % velocityAdvectTwo
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% BGK collisional subroutine (will be replaced by MD) %
% Dummy version, just for testing the dimensionless version, use 1/wp for
% the tau value.
    function BGKCollide(dt)
        if timing, start0 = toc; end
        
        % get moments
        getMoments;
        
        % get measure of tau (use 1/wp for now)
        wp = sqrt(chargeR.^2 .* (pi*n).^(3/2) ./ massR);
        
        % initialize colllision array
        fcoll = zeros(size(f));
        
        % if tau = 1/wp for now
        % compute self-collisions based on this hand-wavy method
        repeatedparams.nR = repmat(n,[1,1,Nv,Nv]);
        repeatedparams.uR = permute(repmat(u,[1,1,1,Nv,Nv]),[1,2,4,5,3]);
        repeatedparams.TR = repmat(T,[1,1,Nv,Nv]);
        fcoll = fcoll + dt*repmat(wp,[1,1,Nv,Nv]) .* ...
            (maxwellian(repeatedparams) - f);
        
        % cross-species collisions
        for i = 1:Nsp
            for j = i+1:Nsp
                % no idea what to put here...YET!
            end
        end

        errorCheck;
        
        if timing, avgcollide = avgcollide + (toc-start0) / Ndt; end
    end % BGKCollide
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Check for Invalid Data %
    function errorCheck
        err = '';
        if sum(isnan(fadv(:))) ~= 0
            err = sprintf('%s%s\n',err,...
            'NaN detected in advective component.');
        end
        if sum(imag(fadv(:)) ~= 0) ~= 0
            err = sprintf('%s%s\n',err,...
                'Imaginary value detected in advective component.');
        end
        if sum(imag(fcoll(:)) ~= 0) ~= 0
            err = sprintf('%s%s\n',err,...
                'Imaginary value detected in collisional component.');
        end
        if isempty(err) ~= 1
            error(err(1:end-1));
        end            
    end % errorCheck
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Make Plots of zeroth through second moments (n, u, T) %
    function plotMoments
        
        % density plot
        nfig = figure(1); clf, hold on
        for i = 1:Nsp
            plot(x,n(i,:),'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('n')
        if ~isnan(axissize)
            axis([0 Lx axissize(1,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Time %.4f', t))
        if saveplot
            saveas(nfig,sprintf('../output/images/nBGK%d.png',tstep),'png')
        end
        
        % x-velocity plot
        ufig = figure(2); clf, hold on
        for i = 1:Nsp
            plot(x,u(i,:,1),'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('u')
        if isnan(axissize)==0
            axis([0 Lx axissize(3,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Time %.4f', t))
        if saveplot
            saveas(ufig,sprintf('../output/images/uBGK%d.png',tstep),'png')
        end
        
        % temperature plot
        Tfig = figure(3); clf, hold on
        for i = 1:Nsp
            plot(x,T(i,:),'Color',colors(i,:))
        end
        xlabel('x')
        ylabel('T')
        if isnan(axissize)==0
            axis([0 Lx axissize(2,:)]);
        end
        legend(legends{1:Nsp})
        title(sprintf('Time %.4f', t))
        if saveplot
            saveas(Tfig,sprintf('../output/images/TBGK%d.png',tstep),'png')
        end
        
        pause(0.05)
    end % plotMoments
%-------------------------------------------------------------------------%

end % BGK_Func
