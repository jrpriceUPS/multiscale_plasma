function output = MD_Func(species, distribution, params)
%%%---------------------------------------------------------------------%%%
%
% 2D Molecular dynamics simulation for the multiscale model. Supports
% n-species in two dimensions.
%
% Uses units in terms of the atomic circle radius, plasma freqency, and
% mass and charge at some reference conditions.
%
% The domain size is determined by the length of the domain, Lx, the number
% of particles, Np, and the average 2D density.
%
% Nearest neighbor lists and a cutoff radius is used to calculate the
% interatomic potential and the forces in ~O(N) time.
%
% Saves data on the set intervals based on "dt," but can take variable time
% steps to prevent unreasonable accelerations when particles get too close.
%
% INPUTS:
%
% Species data:
%   species.Nsp     -   number of species
%   species.Np      -   total number of particles
%   species.mass    -   particle mass of each species
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
%   params.accCoeff -   dt = min(dtmax, accCoeff/max(|acc|))
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
% OUTPUT:
%   output.nList    -   density measured at each cell at each step
%   output.uList    -   velocity measured at each cell at each step
%   output.TList    -   temperature measured at each cell at each step
%   output.HList    -   entropy measured at each cell at each step
%
%   output.kinE     -   total kinetic energy measured at each step
%   output.potE     -   total potential energy measured at each step
%   output.totE     -   total energy measured at each step
%
%   output.dt       -   actual timestep taken
%   output.tActual  -   actual duration of simulation time
%
%%% -----------------------------------------------------------------------

% initialize timing counters
tic
start0 = toc;
timing = params.timing;
if timing
    avgForcesTime   =   0;
    avgBinTime      =   0;
    avgMeasureTime  =   0;
    avgMovieTime    =   0;
end
    

%-------------------------------------------------------------------------%
% SET UP DOMAIN PARAMETERS -----------------------------------------------%
%-------------------------------------------------------------------------%

% particles and species
Nsp     =   species.Nsp;
Np      =   species.Np;
N       =   round(Np * sum(distribution.n,2) / sum(distribution.n(:)));
spIdx   =   [0; cumsum(N)]; % start and end indices for each species

% spatial and velocity parameters
Nx      =   params.Nx;
navg    =   sum(distribution.n(:)) / Nx;
k       =   params.k;
Lx      =   params.Lx;
Ly      =   Np / (Lx*navg); % navg = Np / (Lx*Ly)
Nv      =   params.Nv;
Lv      =   params.Lv;

% time parameters
dt          =   params.dt;
accCoeff    =   params.accCoeff;
Ndt         =   params.Ndt;
Ndteq       =   params.Ndteq;
Ndt_        =   Ndt + Ndteq;
tGlobal     =   0;          % overall timestep (equilibration + simulation)
ttotal      =   0;          % total time elapsed
nextWrite   =   dt;         % time at next write

% desired moments
n0      =   distribution.n;
u0      =   distribution.u;
T0      =   distribution.T;

% nearest neighbor
rcut    =   params.cutoff * k;
Ndiv    =   params.Ndiv;

% equilibration
gamma   =   params.gamma;

% velocity resample
refinement  =   params.refinement;
Nmoments    =   params.Nmoments;
abserr      =   params.abserr;
relerr      =   params.relerr;

% data output parameters
Njumps  =   params.Njumps;
Nsplits =   params.Nsplits;
timing  =   params.timing;
energyFileID    =   fopen(sprintf('../output/%s', params.fnames{1}), 'w');
movieFileID     =   '';
movieNumber     =   0;
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% INITIALIZE VARIABLES ---------------------------------------------------%
%-------------------------------------------------------------------------%

% system
boxSize = [Lx Ly];
vlist   = linspace(-Lv,Lv,Nv);
dv      = vlist(2)-vlist(1);
vgrid   = [-inf,vlist+dv/2];vgrid(end)=inf;
halfBox = 0.5 * boxSize;
dim = 2;

% vectors for particles
mass    =   zeros(Np,1);        % mass              m/m0
charge  =   zeros(Np,1);        % charge            Z/Z0
pSp     =   zeros(Np,1);        % species number (for measurements)
for sp = 1:Nsp
    index           =   spIdx(sp)+1:spIdx(sp+1);
    mass(index)     =   species.mass(sp);
    charge(index)   =   species.charge(sp);
    pSp(index)      =   sp;
end
repmass =   repmat(mass,[1,2]);
pos     =   zeros(Np,dim);      % position          x/a0
vel     =   zeros(Np,dim);      % velocity          v/(a0*w0)
acc     =   zeros(Np,dim);      % acceleration      a/(a0*w0^2)

% observables
potE    =   zeros(Ndt_+1,1);        % potential energy            PE/(a0^2*w0^2*m0)
kinE    =   zeros(Ndt_+1,1);        % kinetic energy              KE/(a0^2*w0^2*m0)
totE    =   zeros(Ndt_+1,1);        % total energy                TE/(a0^2*w0^2*m0)
f_MD    =   zeros(Nv,Nv);           % temporary MD version of f   same as distribution.f
nList   =   zeros(Nsp,Nx,Ndt_+1);   % density                     n*a0^2
uList   =   zeros(Nsp,Nx,2,Ndt_+1); % velocity                    u/(a0*w0)
TList   =   zeros(Nsp,Nx,Ndt_+1);   % temperature                 T/(a0^2*w0^2*m0)
HList   =   zeros(Nsp,Nx,Ndt+1);     % Entropy                     unknown units at present
tList   =   zeros(1,Ndt+1);          % List of actual times

% set integration weights
wtN         =   dv*ones(Nv,1);
wtN(1)      =   0.5*wtN(1);
wtN(end)    =   0.5*wtN(end);
W           =   wtN*wtN.';       % wi*wj

%TEMPORARY
W2   =   permute(repmat(wtN*wtN.',[1,1,Nsp,Nx]),[3,4,1,2]);
WV1 =   permute(repmat(wtN*(wtN.*vlist')',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*vi*wj
WV2 =   permute(repmat((wtN.*vlist')*wtN',[1,1,Nsp,Nx]),[3,4,2,1]); % wi*wj*vj
%TEMPORARY

% multi-scale (cells)
cellSize    =   boxSize(1) / Nx;
nbParticles =   zeros(Nsp,Nx);

% nearest neighbor lists
nbBins  =   floor(Ndiv*boxSize/rcut);   % number bins in each direction
nbBins(nbBins < 3) = 1;                 % must have 3x3 or more to work
Nbins   =   prod(nbBins);               % total number of bins
binSize =   boxSize./nbBins;            % length of bin sides

% construct connectivity matrix
connectivity    =   NNConnectivity;
nInBin          =   zeros(1,Nbins);
% we assume that far fewer particles than this will end up in any bin
binToParticles  =   zeros(ceil(Np/Nbins*10), Nbins);

% equilibration
whichCell   =   zeros(Np,1);    % confine particle to cell to equilibrate
particleU   =   zeros(Np,dim);  % bulk velocity of each particle's cell
particleT   =   zeros(Np,1);    % temperature of each particle's cell
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% INITIALIZE THE SYSTEM --------------------------------------------------%                                                               
%-------------------------------------------------------------------------%

% prepare the initial state based on the input moments
SetInitialState;
fprintf('Setup: %f sec\n', toc - start0)
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% EQUILIBRATION ----------------------------------------------------------%
%-------------------------------------------------------------------------%
start0 = toc;
equilibration_complete = 0;

% determine which cell each particle is in and record the bulk velocity and
% temperature corresponding to that cell
whichCell = ceil(pos(:,1) / cellSize);
for sp = 1:Nsp
    index = spIdx(sp)+1:spIdx(sp+1);
    particleU(index,:) = squeeze(u0(sp,whichCell(index),:));
    particleT(index) = squeeze(T0(sp,whichCell(index)));
end
repT = repmat(particleT,[1,2]);

% initially bin particles and calculate forces
BinParticles;
Forces;
dtuse = min(dt, 0.5/max(sqrt(sum(acc.^2,2))));

% initial measurement and movie frame
Measurements(0);
MovieWrite;

% equilibration loop
while ttotal < Ndteq * dt
    
    if timing == 2, start2 = toc; end
    writeFlag = 0;
    
    % reduce friction halfway through equilibration
    if abs(ttotal - Ndteq * dt / 2) < 1e-10, gamma = gamma / 2; end
    
    % advance particle position and velocity using Langevin equations
    LangevinEvolution;
    ttotal = ttotal + dtuse;
    
    % record measureables
    if abs(ttotal - nextWrite) < 1e-10
        writeFlag = 1;
        nextWrite = nextWrite + dt;
        % track global timetep
        tGlobal = tGlobal + 1;
        Measurements(0);
        fprintf('writing measurement\n')
    end
    
    % record particle data for movie
    if writeFlag && mod(tGlobal, Njumps) == 0
        MovieWrite;
        fprintf('writing movie\n')
    end
    
    % plot u, T averages every 100 timesteps
    if tGlobal >= 200 && mod(tGlobal,50) == 0
        figure(1)
        plot(squeeze(mean(uList(:,:,1,tGlobal-199:tGlobal),4))')
        figure(2)
        plot(squeeze(mean(TList(:,:,tGlobal-199:tGlobal),3))')
        pause(0.05);
    end
    
    if timing == 2
        fprintf('Equilibration time %d: %f sec\n', tGlobal, toc - start2)
        fprintf('dtuse = %f\n', dtuse)
    end
end

% velocity resample with rejection to match distribution
for sp = 1:Nsp
    for ix = 1:Nx
        Filter = (whichCell == ix) & (pSp == sp);
        if sum(Filter)~=0
        vel(Filter,:) = ...
            VelocitySample(squeeze(distribution.f(sp,ix,:,:)),...
            sum(Filter),Lv,refinement,Nmoments,abserr,relerr);
        end
    end
end
ttotal = 0;
equilibration_complete=1;
Measurements(1);
tGlobal = Ndteq + 1;
MovieWrite;
tGlobal = Ndteq;

fprintf('Equilibration: %f sec\n', toc - start0)
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% MAIN SIMULATION --------------------------------------------------------%
%-------------------------------------------------------------------------%
start0 = toc;
nextWrite = dt;
dtuse = min(dt, accCoeff/max(sqrt(sum(acc.^2,2))));
% main loop
while ttotal < Ndt * dt
    if timing == 2, start2 = toc; end
    writeFlag = 0;
    
    % advance particle position and velocity
    Evolution;
    ttotal = ttotal + dtuse;
    
    % record measureables
    if abs(ttotal - nextWrite) < 1e-10
        writeFlag = 1;
        nextWrite = nextWrite + dt;
        % track global timetep
        tGlobal = tGlobal + 1;
        Measurements(1);
        fprintf('writing measurement\n')
    end
    
    % record particles data for movie
    if writeFlag && mod(nextWrite/dt - 1, Njumps) == 0
        MovieWrite;
        fprintf('writing movie\n')
    end
    
    % plot u, T averages every 100 timesteps
    if mod(tGlobal,50) == 0 && tGlobal >= 200
        figure(1)
        plot(squeeze(mean(uList(:,:,1,tGlobal-199:tGlobal),4))')
        figure(2)
        plot(squeeze(mean(TList(:,:,tGlobal-199:tGlobal),3))')
        pause(0.05);
    end
    
    if timing == 2
        fprintf('Simulation time %f: %f sec\n', ttotal, toc - start2)
        fprintf('dtuse = %f\n', dtuse)
    end
end
fprintf('Main Simulation: %f sec\n', toc - start0)
fclose(movieFileID);
if timing
    fprintf('Average time to calculate forces: %f sec\n', avgForcesTime)
    fprintf('Average time to bin particles: %f sec\n', avgBinTime)
    fprintf('Average time to measure moments: %f sec\n', avgMeasureTime)
    fprintf('Average time to write movie file: %f sec\n', avgMovieTime)
end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% OUTPUT -----------------------------------------------------------------%
%-------------------------------------------------------------------------%
output.nList    =   nList;
output.uList    =   uList;
output.TList    =   TList;
output.HList    =   HList;
output.tList    =   tList;
output.kinE     =   kinE;
output.potE     =   potE;
output.totE     =   totE;
output.dt       =   dt;
output.tActual  =   ttotal;
species.N       =   N;
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% NESTED FUNCTIONS -------------------------------------------------------%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Set Initial Particle Positions and Velocities %
    function SetInitialState
        if timing, start1 = toc; end
        
        % determine the number of particles to put in each cell
        nbParticles = n0 ./ repmat(sum(n0,2),[1,Nx]) .* repmat(N,[1,Nx]);
        nbParticles = round(nbParticles);
        
        % if not exactly N particles were placed, add or subtract from the
        % cells with the most particles
        for sp = 1:Nsp
            nbParticlesTmp = nbParticles(sp,:);
            while sum(nbParticles(sp,:)) ~= N(sp)
                [~,idx] = max(nbParticlesTmp);
                nbParticlesTmp(idx) = 0;
                if sum(nbParticles(sp,:)) < N(sp)
                    nbParticles(sp,idx) = nbParticles(sp,idx) + 1;
                elseif sum(nbParticles(sp,:)) > N(sp)
                    nbParticles(sp,idx) = nbParticles(sp,idx) - 1;
                end
            end
            clear nbParticlesTmp;
        end
        
        % put particles into cells
        posIdx  = ones(Nsp,1) + spIdx(1:Nsp);   % track indices in pos
        vel     = randn(size(vel));             % normal distribution    
        nInCell = sum(nbParticles,1);           % total number in each cell
        
        for cell = 1:Nx
            % generate positions using Halton, linearly transform
            posList = zeros(nInCell(cell),2);
            for p = 1:nInCell(cell)
                posList(p,1) = Halton(p,2);
                posList(p,2) = Halton(p,3);
            end
            posList(:,1) = cellSize*posList(:,1) + cellSize*(cell-1);
            posList(:,2) = boxSize(2)*posList(:,2);
            
            % assign positions to particles and scale velocity
            posListIdx = 1;
            for sp = 1:Nsp
                % set position
                index1 = posIdx(sp):(posIdx(sp)+nbParticles(sp,cell)-1);
                index2 = posListIdx:(posListIdx+nbParticles(sp,cell)-1);                
                pos(index1,:) = posList(index2,:);
                
                % scale velocity
                vel(index1,1) = vel(index1,1) .* ...
                    sqrt(T0(sp,cell) ./ mass(index1));
                vel(index1,2) = vel(index1,2) .* ...
                    sqrt(T0(sp,cell) ./ mass(index1));
                
                % update tracking indices
                posIdx(sp) = posIdx(sp) + nbParticles(sp,cell);
                posListIdx = posListIdx + nbParticles(sp,cell);
            end
        end
        
        if timing, fprintf('Set initial state: %f sec\n', toc-start1), end
    end % SetInitialState
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Calculate Forces %
    function Forces
        if timing, start1 = toc; end
        
        % reset acceleration
        acc(:,:) = 0;
        
        % loop over the bins (skipping empty bins)
        for bin = find(nInBin ~= 0)
            % get list of particles to compare
            p1 = binToParticles(1:nInBin(bin),bin);
            p2 = nonzeros(binToParticles(:,connectivity(:,bin)));
            
            % initialize distance and force arrays
            dr = zeros(length(p1),length(p2),dim);
            Fr = zeros(length(p1),length(p2),dim);
            
            % components of relative distance in box
            for d = 1:dim
                dr_d = repmat(pos(p1,d),[1,length(p2)]) - ...
                    repmat(pos(p2,d)',[length(p1),1]);
                Mask = abs(dr_d) > halfBox(d);
                dr_d(Mask) = dr_d(Mask) - sign(dr_d(Mask)) * boxSize(d);
                dr(:,:,d) = dr_d;
            end
            
            % norm of relative distances
            dd = sqrt(sum(dr.^2,3));
            dd(dd > rcut) = 0;
            ddInv = 1 ./ dd;
            ddInv(dd == 0) = 0;
            
            % energy of particle pairs
            Vij = 1/2 * repmat(charge(p1),[1,length(p2)]) .* ...
                repmat(charge(p2)',[length(p1),1]) .* ddInv .* exp(-k*dd);
            
            % save potential energy of system (Vij is double counting)
            potE(tGlobal+1) = potE(tGlobal+1) + 1/2*sum(sum(Vij));
            
            % calculate forces
            Fd = Vij .* (ddInv + k);
            for d = 1:dim
                Fr(:,:,d) = Fd .* (dr(:,:,d) .* ddInv);
            end
            
            % compute acceleration of particles
            acc(p1,:) = reshape(sum(Fr,2),[length(p1),2]) ./ repmass(p1,:);
        end
        
        if timing
            avgForcesTime = avgForcesTime + (toc - start1) / (Ndt_ + 1);
        end
    end % Forces
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Simulation Phase Evolution %
    function Evolution
        % Velocity Verlet time integration
        vel = vel + 1/2*dtuse*acc;
        pos = pos + vel*dtuse;
        
        % Periodic boundary conditions
        for d = 1:dim
            Filter = pos(:,d) < 0 | pos(:,d) > boxSize(d);
            pos(Filter,d) = pos(Filter,d) - sign(pos(Filter,d))*boxSize(d);
        end
        
        % rebin particles and compute the forces
        BinParticles;
        Forces;
        dtuse = min(dt, 0.5/max(sqrt(sum(acc.^2,2))));
        if ttotal + dtuse > nextWrite + 1e-10
            dtuse = nextWrite - ttotal;
        end
        
        % integrate velocity
        vel = vel + 1/2*dtuse*acc;
        
    end % Evolution
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Equilibration Phase Evolution %
    function LangevinEvolution
        % Velocity Verlet time integration for Langevin equation
        vel = vel + 1/2*dtuse*acc - 1/2*dtuse*gamma*vel + ...
            sqrt(dtuse*gamma*repT ./ (repmass)) .* randn(Np,2);
        if sum(abs(dtuse*vel(:) > cellSize)) > 0
            error('Particles moved entire cell length. Reduce timestep.');
        end
        pos = pos + dtuse*vel;
        
        % periodic boundary condition in the y-direction
        Filter = pos(:,2) < 0 | pos(:,2) > boxSize(2);
        pos(Filter,2) = pos(Filter,2) - sign(pos(Filter,2))*boxSize(2);
        
        % record which cell each particle is in after advancing time and
        % determine relections
        whichCell2 = ceil(pos(:,1) / cellSize);
        reflect = whichCell2 - whichCell;
        
        % reflect particles that leave their cell
        FilterR = reflect == 1;     % particle exited right side of cell
        FilterL = reflect == -1;    % particle exited left side of cell
        pos(FilterR,1) = 2*whichCell(FilterR)*cellSize - pos(FilterR,1);
        pos(FilterL,1) = 2*(whichCell(FilterL) - 1)*cellSize - ...
            pos(FilterL,1);
        vel(FilterR,1) = -vel(FilterR,1);
        vel(FilterL,1) = -vel(FilterL,1);
        
        % rebin particles and compute the forces
        BinParticles;
        Forces;
        dtuse = min(dt, 0.5/max(sqrt(sum(acc.^2,2))));
        if ttotal + dtuse > nextWrite - 1e-10
            dtuse = nextWrite - ttotal;
        end
            
        
        % update velocity with Velocity Verlet using Langevin thermostat
        vel = 1 / (1 + 1/2*gamma*dtuse) * (vel + 1/2*dtuse*acc + ...
            sqrt(dtuse*gamma*repT ./ (repmass)) .* randn(Np,2));
        
    end % LangevinEvolution
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Set Up Nearest Neighbor Connectivity %
    function connectivity = NNConnectivity
        % prevent unnecessary evaluations if binSize is factor of rcut
        fix = 1e-10 * (mod(rcut./binSize,1) == 0);
        
        % max y-offset for each x-offset
        xLen = floor(rcut/binSize(1) + 1 - fix(1));
        yLen = floor(sqrt(rcut^2 - ((0:xLen-1)*binSize(1)).^2) / ...
            binSize(2) + 1 - fix(2));
        
        % use symmetry over 4 quadrents and add zero-offset pieces
        Nneighbors = 4*sum(yLen) + 2*xLen + 2*yLen(1) + 1;
        offsets = zeros(Nneighbors,2);
        ctr = 1;
        for i = -xLen:xLen
            iFix = (i == 0);  % "ylen(0)" is same as yLen(1)
            for j = -yLen(abs(i) + iFix):yLen(abs(i) + iFix)
                offsets(ctr,:) = [j i];
                ctr = ctr + 1;
            end
        end
        
        % build connectivity matrix
        connectivity = zeros(Nneighbors,Nbins);
        for bin = 1:Nbins
            i = mod(rem(bin-1,nbBins(1)) + offsets(:,1), nbBins(1)) + 1;
            j = mod(ceil(bin/nbBins(1)) + offsets(:,2) - 1, nbBins(2)) + 1;
            connectivity(:,bin) = (j-1)*nbBins(1) + i;
        end
    end % NNConnectivity
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Put Particles in Bins for Nearest Neighbor %
    function BinParticles
        if timing, start1 = toc; end
        
        % refresh bin lists
        nInBin(:) = 0;
        binToParticles(:,:) = 0;
        
        % place particles in bins
        for p = 1:Np
            bin = ceil(nbBins(1) * pos(p,1) / boxSize(1)) + ...
                floor(nbBins(2) * pos(p,2) / boxSize(2)) * nbBins(1);
            nInBin(bin) = nInBin(bin) + 1;
            binToParticles(nInBin(bin),bin) = p;
        end
        
        if timing
            avgBinTime = avgBinTime + (toc - start1) / (Ndt_ + 1);
        end
    end % BinParticles
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Store Measureable Quantities
    function Measurements(writeflag)
        if timing, start1 = toc; end
        
        tIdx = tGlobal + 1;
        % track moments throughout
        cell = ceil(pos(:,1) / cellSize);
        cellArea = boxSize(2)*cellSize;
        
        for ix = 1:Nx
            for sp = 1:Nsp
                Filter = (cell == ix) & (pSp == sp);
                
                if sum(Filter)~=0
                    tmpN = sum(Filter);
                    nList(sp,ix,tIdx)   = tmpN / cellArea;
                    uList(sp,ix,1,tIdx) = sum(vel(Filter,1)) / tmpN;
                    uList(sp,ix,2,tIdx) = sum(vel(Filter,2)) / tmpN;
                    vmu1 = vel(Filter,1) - squeeze(uList(sp,ix,1,tIdx));
                    vmu2 = vel(Filter,2) - squeeze(uList(sp,ix,2,tIdx));
                    TList(sp,ix,tIdx)   = sum(mass(Filter) .*...
                        (vmu1.^2 + vmu2.^2)) / (2*tmpN - 2*(tIdx < Ndteq));
                    
                    %best way to do this next step would be with hist3, but I
                    %don't have the statistical toolbox. Instead I'll use code
                    %developed externally which might be unreliable.
                    %
                    %I think I might have the units on f_MD correct now,
                    %but I am not certain
                    
                    f_MD                = histcn(vel(Filter,:),vgrid,vgrid)*Nx*navg/Np;
                    f_MD(f_MD==0)       = 1;
                    
                    if equilibration_complete==1
                        HList(sp,ix,tIdx-Ndteq)   = sum(sum(W.*f_MD.*log(f_MD)));
                        tList(tIdx-Ndteq)         = ttotal;
                    end
                
                % naively change recordings in the case of no particles
                % of a particular species in the box. I am not sure if
                % this is correct - Jake
                else
                    nList(sp,ix,tIdx) = 0;
                    uList(sp,ix,1,tIdx) = 0;
                    uList(sp,ix,2,tIdx) = 0;
                    TList(sp,ix,tIdx)   = 0;
                    
                    if equilibration_complete==1
                        HList(sp,ix,tIdx-Ndteq)   = 0;
                        tList(tIdx-Ndteq)         = ttotal;
                    end
                    
                end
                    
            end
        end
        
        % save kinetic, potential, and total energy during main simulation
        kinE(tIdx) = 0.5*sum(sum(repmass.*vel.^2));
        totE(tIdx) = kinE(tIdx) + potE(tIdx);
        if writeflag
            fprintf(energyFileID, '%13.6E %13.6E %13.6E %13.6E\n', ...
                tGlobal-Ndteq, totE(tIdx), kinE(tIdx), potE(tIdx));
        end
        
        if timing
            avgMeasureTime = avgMeasureTime + (toc - start1) / (Ndt_ + 1);
        end
    end % Measurements
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% Output Movie File for Ovito Playback
    function MovieWrite
        % stuff to split movie into equilibration phase plus number of
        % splits (if needed to avoid huge output files in long simulation)
        if timing, start1 = toc; end
        if tGlobal == 0
            movieFileID = fopen(sprintf('../output/MD_Movie%d.dat', ...
                movieNumber), 'w');
            movieNumber = movieNumber + 1;
        elseif tGlobal == Ndteq + 1;
            fclose(movieFileID);
            movieFileID = fopen(sprintf('../output/MD_Movie%d.dat', ...
                movieNumber), 'w');
            movieNumber = movieNumber + 1;
        elseif tGlobal > Ndteq+1 && tGlobal < Ndt_ && ...
                mod(tGlobal-Ndteq,Ndt/Nsplits) == 0
            fclose(movieFileID);
            movieFileID = fopen(sprintf('../output/MD_Movie%d.dat', ...
                movieNumber), 'w');
            movieNumber = movieNumber + 1;
        end
        
        
        fprintf(movieFileID, '%u\n', Np);
        fprintf(movieFileID, '%1s\n', ' ');
        % species, 2D position, 2D velocity, kinetic energy
        for p = 1:Np
            fprintf(movieFileID, ...
                '%4u %13.6E %13.6E %13.6E %13.6E %13.6E\n', ...
                pSp(p), pos(p,1), pos(p,2), vel(p,1), vel(p,2), ...
                1/2*mass(p)*sum(vel(p,:).^2));
        end
        
        if timing
            avgMovieTime = avgMovieTime + (toc - start1) / (Ndt_ / Njumps);
        end
    end % MovieWrite
%-------------------------------------------------------------------------%

end