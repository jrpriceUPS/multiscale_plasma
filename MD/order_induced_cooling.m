function output = order_induced_cooling(species, params, oic)

output=0;

% particles
Np       =   species.Np;

% spatial and velocity parameters
k       =   params.k;
Lx      =   params.Lx;
Ly      =   params.Ly;
Lz      =   params.Lz;
T       =   species.T;

% time parameters
wp2     =   species.charge^2 .* ...
    (pi*Np/(Lx*Ly*Lz)).^(3/2) ./ species.mass;
dt      =   params.dt / sqrt(max(wp2(:)));
Ndt     =   params.Ndt;
tGlobal =   0;
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% INITIALIZE VARIABLES ---------------------------------------------------%
%-------------------------------------------------------------------------%

% system
boxSize = [Lx Ly Lz];
halfBox = 0.5 * boxSize;
dim = 3;

% vectors for particles
mass    =   species.mass*ones(Np,1);        % mass              m/m0
charge  =   species.charge*ones(Np,1);        % charge            Z/Z0

repmass =   repmat(mass,[1,dim]);
pos     =   zeros(Np,dim);      % position          x/a0
vel     =   zeros(Np,dim);      % velocity          v/(a0*w0)
acc     =   zeros(Np,dim);      % acceleration      a/(a0*w0^2)

% observables
potE    =   zeros(Ndt+1,1);        % potential energy  PE/(a0^2*w0^2*m0)
kinE    =   zeros(Ndt+1,1);        % kinetic energy    KE/(a0^2*w0^2*m0)
totE    =   zeros(Ndt+1,1);        % total energy      TE/(a0^2*w0^2*m0)
TList   =   zeros(Ndt+1,1);        % temperature       T/(a0^2*w0^2*m0)

% distribution function
maxr    =   sqrt(Lx^2+Ly^2+Lz^2);
rRange  =   linspace(0,maxr,params.Nr);
deltar      =   rRange(2)-rRange(1);
rCount  =   zeros(size(rRange));
numbersamples   =   params.Ndt - params.record;

% multi-scale (cells)
cellSize    =   boxSize(1);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% INITIALIZE THE SYSTEM --------------------------------------------------%
%-------------------------------------------------------------------------%

% prepare the initial state based on the input moments


if strcmp('halton',oic.method)
    
    h1 = oic.h1;
    h2 = oic.h2;
    h3 = oic.h3;
    SetInitialStateHalton;

elseif strcmp('uniform',oic.method)
    
    SetInitialStateUniform;
    
elseif strcmp('grid',oic.method)
    
    SetInitialStateGrid;
    
end


% initial measurement and movie frame
Forces;
Measurements;
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% MAIN SIMULATION --------------------------------------------------------%
%-------------------------------------------------------------------------%
% main loop
for t = 1:Ndt
    
    % track global timestep
    tGlobal = tGlobal + 1;
    
    % advance particle position and velocity
    Evolution;
    
    % record measureables
    Measurements;
    
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% OUTPUT -----------------------------------------------------------------%
%-------------------------------------------------------------------------%

output.TList    =   TList;
output.kinE     =   kinE;
output.potE     =   potE;
output.totE     =   totE;
output.dt       =   dt;
output.g        =   3*rCount/(Np^2)./((rRange+deltar).^3-rRange.^3);
output.rRange   =   rRange;
output.pos      =   pos;
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Set Initial Particle Positions and Velocities Halton %
    function SetInitialStateHalton
        
        % prevent particles from being too close too eachother on bin edges
        buf     =   0;%1/6;
        scale   =   [cellSize - 2*buf, boxSize(2) - 2*buf, boxSize(3) - 2*buf];
        
        % put particles into cells
        rng(1)
        vel     = randn(size(vel));             % normal distribution
        
        % generate positions using Halton, linearly transform
        for i=1:Np
            pos(i,1) = Halton(i,h1);
            pos(i,2) = Halton(i,h2);
            pos(i,3) = Halton(i,h3);
        end
        
        pos(:,1) = buf + scale(1)*pos(:,1);
        pos(:,2) = buf + scale(2)*pos(:,2);
        pos(:,3) = buf + scale(3)*pos(:,3);
        
        vel(:,1) = vel(:,1) * sqrt(T/species.mass);
        vel(:,2) = vel(:,2) * sqrt(T/species.mass);
        vel(:,3) = vel(:,3) * sqrt(T/species.mass);
    end % SetInitialState
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Set Initial Particle Positions and Velocities Halton %
    function SetInitialStateUniform
        
        % prevent particles from being too close too eachother on bin edges
        buf     =   0;%1/6;
        scale   =   [cellSize - 2*buf, boxSize(2) - 2*buf, boxSize(3) - 2*buf];
        
        % put particles into cells
        rng(1)
        vel     = randn(size(vel));             % normal distribution
        
        % generate positions uniformly
        pos(:,1)=scale(1)*rand(Np,1);
        pos(:,2)=scale(2)*rand(Np,1);
        pos(:,3)=scale(3)*rand(Np,1);
        
        vel(:,1) = vel(:,1) * sqrt(T/species.mass);
        vel(:,2) = vel(:,2) * sqrt(T/species.mass);
        vel(:,3) = vel(:,3) * sqrt(T/species.mass);
    end % SetInitialState
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Set Initial Particle Positions and Velocities Grid %
    function SetInitialStateGrid
        
        if Lx==Ly
            
            ParticlesPerSide = round(Np^(1/3));
            posx=linspace(1,Lx,ParticlesPerSide);
            posx=posx+Lx/ParticlesPerSide/2;
            posx(end)=Lx/ParticlesPerSide/2;
            
            posy=linspace(1,Ly,ParticlesPerSide);  
            posy=posy+Ly/ParticlesPerSide/2;
            posy(end)=Ly/ParticlesPerSide/2;
            
            posz=linspace(1,Lz,ParticlesPerSide);  
            posz=posz+Lz/ParticlesPerSide/2;
            posz(end)=Lz/ParticlesPerSide/2;
            
            
            [x,y,z]=meshgrid(posx,posy,posz);
            
            pos(:,1)=x(:);
            pos(:,2)=y(:);
            pos(:,3)=z(:);
            
        end
        
        % put particles into cells
        rng(1)
        vel     = randn(size(vel));             % normal distribution
        
        vel(:,1) = vel(:,1) * sqrt(T/species.mass);
        vel(:,2) = vel(:,2) * sqrt(T/species.mass);
        vel(:,3) = vel(:,3) * sqrt(T/species.mass);
    end % SetInitialState
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Simulation Phase Evolution %
    function Evolution
        % Velocity Verlet time integration
        vel = vel + 1/2*dt*acc;
        pos = pos + vel*dt;
        
        % Periodic boundary conditions
        for d = 1:dim
            Filter = pos(:,d) < 0 | pos(:,d) > boxSize(d);
            pos(Filter,d) = pos(Filter,d) - sign(pos(Filter,d))*boxSize(d);
        end
        
        % rebin particles and compute the forces
        Forces;
        
        % integrate velocity
        vel = vel + 1/2*dt*acc;
        
    end % Evolution
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Calculate Forces %
    function Forces
        
        % reset acceleration
        acc(:,:) = 0;
        
        for i=1:Np-1
            
            dr = zeros(Np-i,dim);
            
            for d = 1:dim
                dr_d = -pos(i+1:Np,d) + pos(i,d);
                Mask = abs(dr_d) > halfBox(d);
                dr_d(Mask) = dr_d(Mask) - sign(dr_d(Mask)) * boxSize(d);
                dr(:,d) = dr_d;
            end
            
            dd = sqrt(sum(dr.^2,2));
            
            Vij = 1/3 * repmat(charge(i),[Np-i,1]) .* ...
                charge(i+1:Np) ./ dd .* exp(-k*dd);
            
            potE(tGlobal+1) = potE(tGlobal+1) + sum(Vij);
            
            Fd = Vij .* (1./dd + k);
            
            Fr = zeros(Np-i,dim);
            
            for d = 1:dim
                Fr(:,d) = Fd .* (dr(:,d) ./ dd);
            end
            
            acc(i,:) = acc(i,:) + sum(Fr,1)/mass(i);
            acc(i+1:Np,:) = acc(i+1:Np,:) - Fr./repmass(i+1:Np,:);
            
            if tGlobal>params.record
                rCount = rCount + hist(dd,rRange)/numbersamples;
            end
            
        end
        
    end % Forces
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% Store Measureable Quantities
    function Measurements
        tIdx = tGlobal + 1;
        
        TList(tIdx)   = sum(sum(repmass.*vel.^2))/ (3*Np);
        
        % save kinetic, potential, and total energy during main simulation
        kinE(tIdx) = 0.5*sum(sum(repmass.*vel.^2));
        totE(tIdx) = kinE(tIdx) + potE(tIdx);
    end % Measurements
%-------------------------------------------------------------------------%

end