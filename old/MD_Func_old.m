function MD_Func(species, distribution, domain, params, equil)
%
%MD_Func(species, distribution, domain, params)
%
%%% Molecular dynamics part of the multi-scale model ----------------------
%   Uses atomic units.
%   n-species allowed with different masses and charges
%   Uses nearest neighbor list to calculate potential in ~O(N) time
%   Takes the following input parameters (passed as structs):
%
%   species:
%       N         -   number of particles in each species
%       mass      -   mass of each particle species             [kg]
%       charge    -   charge of each particle species
%
%   distribution:
%       f       -   distribution function to initialize with
%       n       -   density     (from 0th moment)               [m^-3]
%       u       -   velocity    (from 1st moment)               [m/s]
%       T       -   temperature (from 2nd moment) (actually kT) [J]
%
%   domain:
%       Lx      -   length of domain (x-direction)              [m]
%       AR      -   aspect ratio of domain (Lx/Ly)
%       lambda  -   electron screening debye length             [m]
%       dt      -   time step to take                           [s]
%       Ndt     -   number of timesteps to take
%       Nx      -   number of space points
%       Nv      -   number of velocity points
%       dv      -   space between velocities
%
%   params:
%       cutoff  -   cutoff radius number of debye lengths
%       spacing -   minimum particle spacing (number of ion-circle radii)
%       Ndiv    -   number of sub-divisions per rcut-length bin
%       Njumps  -   number of jumps between data dumping
%       timing  -   how much timing data to display:
%                   0 - quiet operation (no output)
%                   1 - bulk timing info (setup, avg time on forces, etc)
%                   2 - timing info at every time step
%
%   equil:
%       dteq  -   maximum timestep in equilibration time    [s]
%       gamma -   parameter in Langevin equation
%       Ndteq -   number of time steps in equilibration phase
%
%%% -----------------------------------------------------------------------

%initialize time counters
tic
start = toc;

%open data files
energyFileID = fopen('energy2.dat','w');
movieFileID = fopen('movieMD.dat','w');


%--------------------------------------------------------------------------
%||||| Conversion to Atomic Units |||||------------------------------------
%--------------------------------------------------------------------------

AtomicLengthPerMeter = 5.2917721092e-11;            %length
AtomicTimePerSecond = 2.418884326505e-17;           %time
AtomicMassPerKg = 9.10938291e-31;                   %mass

%--------------------------------------------------------------------------
%||||| Parameters |||||----------------------------------------------------
%--------------------------------------------------------------------------
% Numerical

boxSize = [ domain.Lx/AtomicLengthPerMeter ...
    , domain.Lx/domain.AR/AtomicLengthPerMeter ];  % Dimension and Size of the box
dt = domain.dt/AtomicTimePerSecond;                % Time-step
Ndt = domain.Ndt;                                  % Number of integration steps
Njumps = params.Njumps;                            % Jumps between measurements

dteq = equil.dteq/AtomicTimePerSecond;             % Equilibration time step
Ndteq = equil.Ndteq;                               % Number of equilibration time steps
gamma = equil.gamma;                               % Equilibration parameter gamma


% Particles

nSp = length(species.mass);             % Number of species types
N = species.N;                          % Total number of species
Np = sum(N);                            % Total number of particles

% Mass and charge of each species
masslist =  species.mass/AtomicMassPerKg;
chargelist = species.charge;

lambda = domain.lambda/AtomicLengthPerMeter;     % Debye length

T = distribution.T;  % kT (J)

% Multi-Scale

Nx = domain.Nx;
Nv = domain.Nv;
%Lv = domain.Lv;
dv = domain.dv;
%dvfine = domain.dvfine;
dvfine = dv/100;
Lv=3.1093e+06;

% Nearest neighbor
cutoff = params.cutoff/AtomicLengthPerMeter; % cutoff radius for comparisons
Ndiv = params.Ndiv;                          % Number bin subdivisions
safetyfactor = params.safetyfactor;          % particles to allow for in bin


%--------------------------------------------------------------------------
%||||| Initialize |||||----------------------------------------------------
%--------------------------------------------------------------------------
% System
halfBox = 0.5*boxSize;
Dim = length(boxSize);                  % Dimension of the system
%
% Vectors for particles
mass   = zeros(Np,1);                   % mass
charge = zeros(Np,1);                   % charge
pos = zeros(Np,Dim);                    % position r(x,y,z)
vel = zeros(Np,Dim);                    % velocity v(x,y,z)
acc = zeros(Np,Dim);                    % acceleration
%
%generate vectors the same size as the number of particles containing the
%mass and charges of the corresponding particle
currentparticle=0;
for i=1:nSp
    mass(currentparticle+1:currentparticle+N(i))=masslist(i);
    charge(currentparticle+1:currentparticle+N(i))=chargelist(i);
    currentparticle=currentparticle+N(i);
end

%
% Observables
potentialE = 0.0;
kinE = 0.0;
potE = 0.0;
totE = 0.0;
density = zeros(nSp,Nx);
velocity = zeros(nSp,Nx,2);
temperature = zeros(nSp,Nx);

% Multi-Scale
cellSize = boxSize(1)./Nx;
cellCenters = (0:Nx-1).*cellSize + 0.5*cellSize - ...
    0.5*boxSize(1);
nbOfPart=[];


% Nearest neighbors
nbBins = floor(Ndiv*boxSize/cutoff);    % Number bins in each direction
nbBins(nbBins < 3) = 1;                 % must have 3x3 or more to work
Nbins = prod(nbBins);                   % total bins
lbin = boxSize./nbBins;                 % length of bin sides

%construct connectivity matrix
[~, connectivity] = NNconnectivity;
nInBin = zeros(1,Nbins);
% assumption: this is >> number of particles than will end up in bin
binToParticles = zeros(ceil(Np/Nbins*safetyfactor), Nbins);

%--------------------------------------------------------------------------
%----- End initialize -----------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%||||| MAIN |||||----------------------------------------------------------
%--------------------------------------------------------------------------
% Prepare the new initial state according to input moments
SetNewInitialState

%note which cell each particle resides in initially (during equilibration,
%particles are confined to their starting cell)
whichcell=zeros(Np,1);

%for each particle, record the bulk velocity and temperature they should
%have based upon their location initially
bulkvelocity=zeros(Np,2);
CellT=zeros(Np,1);
currentspecies=0;

if nSp==1
    for i=1:Nx
        index=cellCenters(i)-cellSize/2<pos(:,1)&pos(:,1)<cellCenters(i)+cellSize/2;
        whichcell(index)=i;
        bulkvelocity(index,1)=distribution.u(i,1);
        bulkvelocity(index,2)=distribution.u(i,2);
        CellT(index)=T(i);
    end
else
    for j=1:nSp
        for i=1:Nx
            index2=cellCenters(i)-cellSize/2<pos(currentspecies+1:currentspecies+N(j),1)&pos(currentspecies+1:currentspecies+N(j),1)<cellCenters(i)+cellSize/2;
            index=zeros(Np,1);
            index(currentspecies+1:currentspecies+N(j))=index2;
            index=logical(index);
            whichcell(index)=i;
            bulkvelocity(index,1)=distribution.u(j,i,1);
            bulkvelocity(index,2)=distribution.u(j,i,2);
            CellT(index)=T(j,i);
        end
        currentspecies=currentspecies+N(j);
    end
end


currentspecies=0;

% Initial call to forces
forcesTime = 0;
BinParticles
Forces

if nSp==1
    temperaturelist=zeros(Ndteq,Nx);
else
    temperaturelist=zeros(Ndteq,nSp,Nx);
end

%run equilibration phase to place all particles such that their
%distribution matches the desired density, bulk velocity, and temperature
MovieWrite
Equilibrate

%the following serves only to visualize the results of equilibration to
%make sure it functioned correctly, the final code will have this removed

if nSp==1
    for i=1:Nx
        index=cellCenters(i)-cellSize/2<pos(:,1)&pos(:,1)<cellCenters(i)+cellSize/2;
        density(i)=sum(index);
        velocity(i)=mean(vel(index,1));
        velocity2=mean(vel(index,2));
        temperature(i)=masslist/(2*density(i)-2)...
            *(sum((vel(index,1)-velocity(i)).^2)+sum((vel(index,2)-velocity2).^2));
    end
else
    currentspecies=0;
    for j=1:nSp
        for i=1:Nx
            index2=cellCenters(i)-cellSize/2<pos(currentspecies+1:currentspecies+N(j),1)&pos(currentspecies+1:currentspecies+N(j),1)<cellCenters(i)+cellSize/2;
            index=zeros(Np,1);
            index(currentspecies+1:currentspecies+N(j))=index2;
            index=logical(index);
            density(j,i)=sum(index);
            velocity(j,i,1)=mean(vel(index,1));
            velocity(j,i,2)=mean(vel(index,2));
            temperature(j,i)=masslist(j)/(2*density(j,i)-2)...
                *(sum((vel(index,1)-velocity(j,i,1)).^2)+sum((vel(index,2)-velocity(j,i,2)).^2));
        end
        currentspecies=currentspecies+N(j);
    end
end


%VelocitySample;

if nSp==1
    for i=1:Nx
        index=cellCenters(i)-cellSize/2<pos(:,1)&pos(:,1)<cellCenters(i)+cellSize/2;
        density(i)=sum(index);
        velocity(i)=mean(vel(index,1));
        velocity2=mean(vel(index,2));
        temperature(i)=masslist/(2*density(i)-2)...
            *(sum((vel(index,1)-velocity(i)).^2)+sum((vel(index,2)-velocity2).^2));
    end
else
    currentspecies=0;
    for j=1:nSp
        for i=1:Nx
            index2=cellCenters(i)-cellSize/2<pos(currentspecies+1:currentspecies+N(j),1)&pos(currentspecies+1:currentspecies+N(j),1)<cellCenters(i)+cellSize/2;
            index=zeros(Np,1);
            index(currentspecies+1:currentspecies+N(j))=index2;
            index=logical(index);
            density(j,i)=sum(index);
            velocity(j,i,1)=mean(vel(index,1));
            velocity(j,i,2)=mean(vel(index,2));
            temperature(j,i)=masslist(j)/(2*density(j,i)-2)...
                *(sum((vel(index,1)-velocity(j,i,1)).^2)+sum((vel(index,2)-velocity(j,i,2)).^2));
        end
        currentspecies=currentspecies+N(j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Equilibration complete%
%%%%%%%%%%%%%%%%%%%%%%%%

%initial call to forces
Forces
time = 0;
fprintf('End of equilibration phase: %f sec\n', toc - start)
forcesTime = 0;

%loop through the time steps
for i=1:Ndt
    start = toc;
    
    %advance time
    time = time +1;
    
    %advance particle position and velocity
    Evolution;
    
    %record measurable quantities
    Measurements;
    
    %record particles for movie
    if(mod(time,Njumps) == 0)
        MovieWrite;
    end
    
    %print time for step
    fprintf('End of step %d: %f sec\n', time, toc - start);
end
toc
fprintf('time spent on binning and forces: %f sec\n', forcesTime)
%
% Compute final distribution function
GetFinalDistribution
%--------------------------------------------------------------------------
%----- End MAIN -----------------------------------------------------------
%--------------------------------------------------------------------------


%----------------------------------------------------------------------
%||||| NNconnectivity |||||--------------------------------------------
%----------------------------------------------------------------------
    function [Nneighbors, connectivity] = NNconnectivity
        %%% build nearest neighbor connectivity %%%
        % prevent unnecessary evaluations if lbin is factor of rcut
        fix = 1e-10 * (mod(cutoff./lbin,1) == 0);
        
        % get max y-offset for each x-offset
        xLen = floor(cutoff/lbin(1) + 1 - fix(1));
        yLen = floor(sqrt(cutoff^2 - ((0:xLen-1)*lbin(1)).^2) / lbin(2) + ...
            1 - fix(2));
        
        % use symmetry over 4 quadrants and add zero-offset components
        Nneighbors = 4*sum(yLen) + 2*xLen + 2*yLen(1) + 1;
        offsets = zeros(Nneighbors,2);
        n = 1;
        for ix = -xLen:xLen
            ixFix = (ix == 0); % "yLen(0)" is same as yLen(1)
            for jx = -yLen(abs(ix)+ixFix):yLen(abs(ix)+ixFix)
                offsets(n,:) = [jx ix];
                n = n+1;
            end
        end
        
        % build connectivity matrix (modular arithmetic to represent bins
        % as 1-indexed 1D vector)
        connectivity = zeros(Nneighbors,Nbins);
        for bin = 1:Nbins
            ix = mod(rem(bin-1,nbBins(1)) + offsets(:,1), nbBins(1)) + 1;
            jx = mod(ceil(bin/nbBins(1)) + offsets(:,2) - 1, ...
                nbBins(2)) + 1;
            connectivity(:,bin) = (jx-1)*nbBins(1) + ix;
        end
    end
%----------------------------------------------------------------------
%----- End NNconnectivity ---------------------------------------------
%----------------------------------------------------------------------



%----------------------------------------------------------------------
%||||| SetNewInitialState |||||----------------------------------------
%----------------------------------------------------------------------
    function SetNewInitialState
        
        % Determine the number of particles to put in each cell
        if nSp==1
            ratioOfPart = distribution.n/sum(distribution.n)*Np;
        else
            ratioOfPart = distribution.n./(squeeze(repmat(sum(distribution.n,2),[1,1,Nx]))).*squeeze(repmat(N,[1,1,Nx]));
        end
        nbOfPart = round(ratioOfPart);
        
        %if not exactly N particles were placed, add or subtract from the
        %cells with the most particles
        
        for i=1:nSp
            nbOfPartTemp = nbOfPart(i,:);
            while(sum(nbOfPart(i,:)) ~= N(i))
                [~,idx] = max(nbOfPartTemp);
                if(sum(nbOfPart(i,:)) < N(i))
                    nbOfPartTemp(idx) = 0;
                    nbOfPart(i,idx) = nbOfPart(i,idx) +1;
                elseif(sum(nbOfPart(i,:)) > N(i))
                    nbOfPartTemp(idx) = 0;
                    nbOfPart(i,idx) = nbOfPart(i,idx) -1;
                end
            end
            clear nbOfPartTemp;
        end
        
        totalperbin=sum(nbOfPart,1);
        
        buffer=1/100*cellSize;
        b=cellSize/2-buffer;
        
        buffer=1/100*boxSize(2);
        c=boxSize(2)/2-buffer;
        
        
        currentparticles=zeros(1,nSp);
        
        if nSp~=1
            for i=1:nSp
                currentparticles(i)=sum(N(1:i-1));
            end
        end
        
        vel = randn(size(vel));                 % Normal velocity distribution
        
        for x=1:Nx
            poslist=zeros(totalperbin(x),2);
            for j=1:totalperbin(x)
                poslist(j,1)=Halton(j,2);
                poslist(j,2)=Halton(j,3);
            end
            poslist(:,1)=2*b*poslist(:,1)+cellCenters(x)-b;
            poslist(:,2)=2*c*poslist(:,2)-c;
            
            poslistposition=0;
            for i=1:nSp
                index=currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x);
                index2=poslistposition+1:poslistposition+nbOfPart(i,x);
                
                pos(index,:) = poslist(index2,:);
                
                if nSp==1
                    vel(index,1)=vel(index,1).*sqrt(distribution.T(x)/masslist(i))+distribution.u(x,1);
                    vel(index,2)=vel(index,2).*sqrt(distribution.T(x)/masslist(i))+distribution.u(x,2);
                else
                    
                    vel(index,1)=vel(index,1).*sqrt(distribution.T(i,x)/masslist(i))+distribution.u(i,x,1);
                    vel(index,2)=vel(index,2).*sqrt(distribution.T(i,x)/masslist(i))+distribution.u(i,x,2);
                end
                poslistposition=poslistposition+nbOfPart(i,x);
                currentparticles(i)=currentparticles(i)+nbOfPart(i,x);
            end
        end
    end
%----------------------------------------------------------------------
%----- End SetNewInitialState -----------------------------------------
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%||||| BinParticles |||||----------------------------------------------
%----------------------------------------------------------------------
    function BinParticles
        
        %refresh bin lists
        nInBin(:) = 0;
        binToParticles(:,:) = 0;
        
        %place particles in appropriate bins
        for pn = 1:Np
            bin = ceil(nbBins(1)*(pos(pn,1)+halfBox(1))/boxSize(1)) + ...
                floor(nbBins(2)*(pos(pn,2)+halfBox(2))/boxSize(2)) * ...
                nbBins(1);
            nInBin(bin) = nInBin(bin) + 1;
            binToParticles(nInBin(bin),bin) = pn;
        end
    end
%----------------------------------------------------------------------
%----- End BinParticles -----------------------------------------------
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%||||| Forces |||||----------------------------------------------------
%----------------------------------------------------------------------
    function Forces
        start0 = toc;
        
        % Initialize variables
        potentialE = 0.0;
        acc = zeros(Np,Dim);
        
        % Loop over the bins (skip empty bins)
        for bin = find(nInBin ~= 0)
            % get particles to be compared
            p1 = binToParticles(1:nInBin(bin),bin);
            p2 = nonzeros(binToParticles(:,connectivity(:,bin)));
            
            % initialize matrices
            dr = zeros(length(p1),length(p2),Dim);
            Fr = zeros(length(p1),length(p2),Dim);
            
            % components of relative distance in box
            for d = 1:Dim
                dr(:,:,d) = repmat(pos(p1,d),1,length(p2)) - ...
                    repmat(pos(p2,d)',length(p1),1);
                for ix = 1:length(p1)
                    Mask = abs(dr(ix,:,d)) > halfBox(d);
                    dr(ix,Mask,d) = dr(ix,Mask,d) - ...
                        sign(pos(p1(ix),d))*boxSize(d);
                end
            end
            
            % norm of relative distancces
            dd = sqrt(sum(dr.^2,3));
            dd(dd > cutoff) = 0;
            ddInv = 1 ./ dd;
            ddInv(dd == 0) = 0;
            
            % energy of particle pairs
            Vij = repmat(charge(p1),1,length(p2)) .* ...
                repmat(charge(p2)',length(p1),1) .* ...
                exp(-dd/lambda) .* ddInv;
            
            % save total energy of system (Vij double counts it)
            potentialE = potentialE + sum(sum(Vij))/2;
            
            % amplitudes of forces
            Fd = Vij .* (ddInv + 1/lambda);
            for d = 1:Dim
                Fr(:,:,d) = Fd .* (dr(:,:,d) .* ddInv);
            end
            
            %compute acceleration of particles
            acc(p1,:) = reshape(sum(Fr,2),size(acc(p1,:))) ...
                ./ repmat(mass(p1),1,2);
        end
    end
%----------------------------------------------------------------------
%----- End force ------------------------------------------------------
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%||||| Equilibrate |||||-----------------------------------------------
%----------------------------------------------------------------------
    function Equilibrate
        %equilibration stage to bring particles into correct distribution
        for i=1:Ndteq
            %reset reflect so we definitely go into while loop
            reflect=2;
            i
            %reset dtuse to its default for equilibration phase
            dtuse=10*dteq;
            
            %complete a step that does not have a particle moving enough to
            %pass through multiple cells
            while sum(reflect>1 | reflect<-1)>0
                %if particles have moved through multiple cells, decrease
                %the time step by a factor of ten and try again
                dtuse=dtuse/10;
                
                %use Verlet velocity algorithm to compute new half velocity
                %and new position
                %
                %velocity updated with Langevin thermostat such that
                %particles in a particular cell have the correct mean
                %velocity and move toward correct temperature
                velnew = vel+dtuse/2*(acc-gamma*(vel-bulkvelocity))+sqrt(...
                    dtuse*gamma*squeeze(repmat(CellT,1,1,2))./squeeze(repmat(mass,1,1,2))...
                    ).*randn(Np,2);
                posnew = pos +velnew*dtuse;
                
                
                %periodic boundary conditions in y direction
                Filter = abs(posnew(:,2)) > halfBox(2);
                posnew(Filter,2) = posnew(Filter,2) -sign(posnew(Filter,2))*boxSize(2);
                
                %record which cell each particle is in after the move
                whichcell2=zeros(Np,1);
                for j=1:Nx
                    index=cellCenters(j)-cellSize/2<posnew(:,1)&posnew(:,1)<cellCenters(j)+cellSize/2;
                    whichcell2(index)=j;
                end
                
                %periodic boundary conditions in x direction (the left
                %boundary is handled by the zero initialization of whichcell2)
                whichcell2(posnew(:,1)>cellCenters(Nx)+cellSize/2)=Nx+1;
                
                %mark which particles have moved into a new cell (they
                %should actually be reflected
                reflect=whichcell2-whichcell;
                
            end
            
            %once we successfully move particles, mark new velocities and
            %positions
            vel=velnew;
            pos=posnew;
            
            %compute which need to be reflected backwards and forwards
            Filter1=reflect==1;
            Filter2=reflect==-1;
            
            %compute left endpoint for use in reflections
            leftend=cellCenters(1)-cellSize/2;
            
            %reflect backwards
            if sum(Filter1~=0)
                pos(Filter1,1)=2*(leftend+whichcell(Filter1)*cellSize)-pos(Filter1,1);
            end
            
            %reflect forwards
            if sum(Filter2~=0)
                pos(Filter2,1)=2*(leftend+(whichcell(Filter2)-1)*cellSize)-pos(Filter2,1);
            end
            
            %place particles in bins
            BinParticles;
            
            %compute forces at new positions
            Forces;
            
            %update velocity with Verlet velocity algorithm using Langevin
            %thermostat again
            vel =1./(1+gamma*dtuse./2).*...
                (vel+dtuse/2*acc+dtuse/2*gamma*bulkvelocity...
                +sqrt(dtuse*gamma*squeeze(repmat(CellT,1,1,2))./squeeze(repmat(mass,1,1,2))).*randn(Np,2));
            
            % write to the movie file
            MovieWrite;
            
            %compute the moments in each cell for comparison against
            %desired values
            %
            %this section might be removed in the final code once we know
            %it works right, or we might compare these computed values
            %against the desired ones to decide when to stop the
            %equilibration phase - density and velocity seem to not change
            %at all at present, which is what we hope for
            if nSp==1
                for k=1:Nx
                    index=cellCenters(k)-cellSize/2<pos(:,1)&pos(:,1)<cellCenters(k)+cellSize/2;
                    density(k)=sum(index);
                    velocity(k)=mean(vel(index,1));
                    velocity2=mean(vel(index,2));
                    temperature(k)=masslist/(2*density(k)-2)...
                        *(sum((vel(index,1)-velocity(k)).^2)+sum((vel(index,2)-velocity2).^2));
                end
            else
                currentspecies=0;
                for j=1:nSp
                    for k=1:Nx
                        index2=cellCenters(k)-cellSize/2<pos(currentspecies+1:currentspecies+N(j),1)&pos(currentspecies+1:currentspecies+N(j),1)<cellCenters(k)+cellSize/2;
                        index=zeros(Np,1);
                        index(currentspecies+1:currentspecies+N(j))=index2;
                        index=logical(index);
                        density(j,k)=sum(index);
                        velocity(j,k,1)=mean(vel(index,1));
                        velocity(j,k,2)=mean(vel(index,2));
                        temperature(j,k)=masslist(j)/(2*density(j,k)-2)...
                            *(sum((vel(index,1)-velocity(j,k,1)).^2)+sum((vel(index,2)-velocity(j,k,2)).^2));
                    end
                    currentspecies=currentspecies+N(j);
                end
            end
            

            %figure(3)
            %plot(velocity(2,:,1));
            %figure(4)
            %plot(velocity(2,:,2));
            %figure(5)
            %plot(squeeze(temperature(1,:)))
            %figure(6)
            %plot(squeeze(temperature(2,:)))
            
            
            if mod(i,50)==0
                figure(1)
                plot(squeeze(temperaturelist(i-50+1:i,1,:)))
                figure(2)
                plot(squeeze(temperaturelist(i-50+1:i,2,:)))
            end
            
            %record temperature for plotting
            if nSp==1
                temperaturelist(i,:)=temperature;
            else
                temperaturelist(i,:,:)=temperature;
            end
        end
    end
%----------------------------------------------------------------------
%----- End Equilibrate ------------------------------------------------
%----------------------------------------------------------------------

%----------------------------------------------------------------------
%||||| VelocitySample |||||--------------------------------------------
%----------------------------------------------------------------------
    function VelocitySample
        
        xlist=-Lv:dv:Lv;
        xfine=-Lv:dvfine:Lv;
        
        vel=zeros(size(vel));
        
        currentparticles=zeros(1,nSp);
        
        if nSp~=1
            for i=1:nSp
                currentparticles(i)=sum(N(1:i-1));
            end
        end
        
        for x=1:Nx
            for i=1:nSp
                if nSp==1
                    
                    
                    
                    vel(currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x),1)=vel(currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x),1).*sqrt(distribution.T(x)/masslist(i))+distribution.u(x,1);
                    vel(currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x),2)=vel(currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x),2).*sqrt(distribution.T(x)/masslist(i))+distribution.u(x,2);
                else
                    
                    velmat=squeeze(distribution.f(i,x,:,:));
                    velmatxraw=sum(velmat,1);
                    velmatyraw=sum(velmat,2);
                    
                    velmatx=interp1(xlist,velmatxraw,xfine);
                    velmaty=interp1(xlist,velmatyraw,xfine);
                    
                    
                    
                    index=currentparticles(i)+1:currentparticles(i)+nbOfPart(i,x);
                    
                    vel(index,1)=xfine(randp(velmatx,length(index),1));
                    vel(index,2)=xfine(randp(velmaty,length(index),1));
                end
                currentparticles(i)=currentparticles(i)+nbOfPart(i,x);
            end
        end
        
    end
%----------------------------------------------------------------------
%----- End VelocitySample ---------------------------------------------
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%||||| Measurements |||||----------------------------------------------
%----------------------------------------------------------------------
    function Measurements
        
        %compute kinetic energy
        kin = zeros(Np,Dim);
        for d=1:Dim
            kin(:,d) = 0.5*mass.*vel(:,d).^2;
        end
        kinE = sum(sum(kin));
        
        %save potential energy
        potE = potentialE;
        
        %save total energy (this should be conserved)
        totE = kinE +potE;
        fprintf(energyFileID, '%13.6E %13.6E %13.6E %13.6E\n', ...
            time, totE, kinE, potE);
    end
%----------------------------------------------------------------------
%----- End measurements -----------------------------------------------
%----------------------------------------------------------------------


%----------------------------------------------------------------------
%||||| MovieWrite |||||------------------------------------------------
%----------------------------------------------------------------------
    function MovieWrite
        
        %record data for movie software Ovito
        fprintf(movieFileID, '%u\n', Np);
        fprintf(movieFileID, '%1s\n', ' ');
        for ii=1:N(1)
            fprintf(movieFileID, '%3s  %13.6E %13.6E %13.6E\n',...
                '1', pos(ii,1), pos(ii,2), 1/2*mass(ii)*sum(vel(ii,:).^2));
        end
        for ii=N(1)+1:Np
            fprintf(movieFileID, '%3s  %13.6E %13.6E %13.6E\n',...
                '2', pos(ii,1), pos(ii,2), 1/2*mass(ii)*sum(vel(ii,:).^2));
        end
    end
%----------------------------------------------------------------------
%----- End MovieWrite -------------------------------------------------
%----------------------------------------------------------------------



%----------------------------------------------------------------------
%||||| Evolution |||||-------------------------------------------------
%----------------------------------------------------------------------
    function Evolution
        % Velocity Verlet time-integration
        vel = vel +0.5*acc*dt;
        pos = pos +vel*dt;
        %
        % Periodic boundary conditions
        for d=1:Dim
            Filter = abs(pos(:,d)) > halfBox(d);
            pos(Filter,d) = pos(Filter,d) -sign(pos(Filter,d))*boxSize(d);
        end
        %
        start0 = toc;
        
        %rebin particles, compute forces again
        BinParticles
        Forces
        forcesTime = forcesTime + toc-start0;
        %
        vel = vel +0.5*acc*dt;
    end
%----------------------------------------------------------------------
%----- End evolution --------------------------------------------------
%----------------------------------------------------------------------

end