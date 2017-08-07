function MD_NearestNeighbor(~)
    clear all, close all
    tic
    start = toc;

% Molecular Dynamics part of the Multi-Scale model
% Uses atomic units
% Only 2 different species are currently allowed (but easy to upgrade)
% The displayed distribution function hax x-position in horizontal axis and
% Vx-velocity in vertical axis (increasing downward).
% Nearest neighbor potential calculation for ~O(N) time

%--------------------------------------------------------------------------
%||||| Constants |||||-----------------------------------------------------
%--------------------------------------------------------------------------
    kB = 3.16499e-06;           % Boltzmann constant (H.K-1)

%--------------------------------------------------------------------------
%||||| Parameters |||||----------------------------------------------------
%--------------------------------------------------------------------------
% Numerical
                                            % Dimension and Size of the box
    boxSize = [ 1 , 1 ];                                    % 2D version
    %boxSize = [ 1 , 1 , 1 ];                               % 3D version
    dt = 1e-05;                             % Time-step
    Ndt = uint16( 5000 );                   % Number of integration steps
    Njumps = uint16( 10 );                  % Jumps between measurements
    
% Particles
                                            % Number of particles in each species
    Nsp = [ 10000 ];                        % 1 species
    %Nsp = uint16( [ 1000 , 2000 ] );           % 2 species
    Ns = size(Nsp,2);                       % Total number of species
    Np = sum(Nsp);                          % Total number of particles
                                            % Mass and charge of each species
    mSp =  [1.0];                           % 1 species
    chSp = [1.0]; 
%     mSp =  [1.0 , 2.0];                                       % 2 species
    chSp = [1.0 , 1.0]; 
    lambda = 1.0e-2;                        % Debye length
    
% System

    Temp = 1.0e4;                           % Temperature (K)
    MakeIniCondMD = 1;                      % Set to 1 to create initial conditions
    gap = 4.0e-3;                           % To prepare Initial Conditions
    
% Multi-Scale

    nbSpacePts = 32;
    nbVelPts = 60;
    deltaVel = 1.0;
    
% Nearest neighbor
    rcut = 4*lambda;                        % cutoff radius for comparisons
    Ndiv = 1;                               % Number bin subdivisions
    safetyfactor = 10;                      % particles to allow for in bin

    
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This is the distribution function which is supposed to be given from the
% KT part :
% 
distFunc = zeros(nbSpacePts,nbVelPts);
for x=1:nbSpacePts
    for v=1:nbVelPts
        distFunc(x,v) = exp(-0.005*((x-(nbSpacePts+1)/2.0)^2+4.0*(v-15.0)^2));
    end
end
distFunc = distFunc/sum(sum(distFunc));
% figure('name','Theoretical distribution function');
% imagesc(transpose(distFunc), [0 4.0e-3]);
velRange = size(distFunc(2));
for v=1:nbVelPts
    velRange(v) = (double(v) -double(nbVelPts+1)/2.0)*deltaVel;
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%--------------------------------------------------------------------------
%----- End parameters -----------------------------------------------------
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%||||| Initialize |||||----------------------------------------------------
%--------------------------------------------------------------------------
% System
    halfBox = 0.5*boxSize;
    Dim = uint16( size(boxSize,2) );        % Dimension of the system
%
% Vectors for particles
    mass   = zeros(Np,1);                   % mass
    charge = zeros(Np,1);                   % charge
    pos = zeros(Np,Dim);                    % position r(x,y,z)
    vel = zeros(Np,Dim);                    % velocity v(x,y,z)
    acc = zeros(Np,Dim);                    % acceleration
%
    if(size(Nsp,2) == 1)
        mass(1:Nsp(1))      = mSp(1);
        charge(1:Nsp(1))    = chSp(1);
    else
        mass(1:Nsp(1))      = mSp(1);
        mass(Nsp(1)+1:Np)   = mSp(2);
        charge(1:Nsp(1))    = chSp(1);
        charge(Nsp(1)+1:Np) = chSp(2);
    end
%
% Open files
    energyFileId = fopen('energy2.dat','w'); 
    movieFileId = fopen('movieMD.dat','w');
%
% Observables
    potentialE = 0.0;
    kinE = 0.0;
    potE = 0.0;
    totE = 0.0;

% Multi-Scale
    cellSize = boxSize(1)./nbSpacePts;
%     cellCenters = zeros(nbSpacePts,1);
%     for i=1:nbSpacePts
%        cellCenters(i) = (i-1)*cellSize +0.5*cellSize -0.5*boxSize(1);
%     end
    cellCenters = (0:nbSpacePts-1).*cellSize + 0.5*cellSize - ...
        0.5*boxSize(1);
    density = zeros(nbSpacePts,1);
    velocity = zeros(nbSpacePts,1);
    temperature = zeros(nbSpacePts,1);

% Nearest neighbors
    nbBins = floor(Ndiv*boxSize/rcut);      % Number bins in each direction
    nbBins(nbBins < 3) = 1;                 % must have 3x3 or more to work
    Nbins = prod(nbBins);                   % total bins
    lbin = boxSize./nbBins;                 % length of bin sides
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
    % Prepare the new initial state    
    SetNewInitialState  

    % Initial call to forces
    forcesTime = 0;
    BinParticles
    Forces
    % Time evolution
%     figure('name','Time evolution in x-y plane')
%     Display
    time = 0;
    fprintf('End of setup phase: %f sec\n', toc - start)
    forcesTime = 0;
    for i=1:Ndt
       start = toc;
       time = time +1;
       Evolution;
       Measurements;
       % Display particles
       if(mod(time,Njumps) == 0)
           %fprintf('Time: %d', time)
           %Display;
           MovieWrite;
       end
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






%__________________________________________________________________________
%--------------------------------------------------------------------------
%-------------------------- NESTED FUNCTIONS ------------------------------
%--------------------------------------------------------------------------
%__________________________________________________________________________


    %----------------------------------------------------------------------
    %||||| Display |||||-------------------------------------------------
    %----------------------------------------------------------------------
    function Display
        hold off;
        scatter(pos(1:Np/3,1),pos(1:Np/3,2),'filled','r')
        hold on;
        scatter(pos(Np/3+1:Np,1),pos(Np/3+1:Np,2),'filled','b')
        axis([-halfBox(1) halfBox(1) -halfBox(2) halfBox(2)])
        pause(0.1);  
    end
    %----------------------------------------------------------------------
    %----- End display --------------------------------------------------
    %----------------------------------------------------------------------
    
    
        %----------------------------------------------------------------------
    %||||| Evolution |||||-------------------------------------------------
    %----------------------------------------------------------------------
    function Evolution
    % Velocity Verlet time-integration
        vel = vel +0.5*acc*dt;
        pos = pos +vel*dt;
    %
    % Periodic boundary conditions on d-1 dimension
        for d=1:Dim %-1
            Filter = abs(pos(:,d)) > halfBox(d);
            pos(Filter,d) = pos(Filter,d) -sign(pos(Filter,d))*boxSize(d);
        end
%   % Reflective boundary conditions on the remaining dimension
%       Mask = abs(pos(:,Dim)) > halfBox(Dim);
%       vel(Mask,Dim) = -vel(Mask,Dim);
	%
        start0 = toc;
        BinParticles
        Forces
        forcesTime = forcesTime + toc-start0;
	%
        vel = vel +0.5*acc*dt;    
    end
    %----------------------------------------------------------------------
    %----- End evolution --------------------------------------------------
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    %||||| Forces |||||----------------------------------------------------
    %----------------------------------------------------------------------
    function Forces
        start0 = toc;
    % Initialize variables
        potentialE = 0.0;
        acc = zeros(Np,Dim);            % acceleration
    % Loop over the bins (skip empty bins)
        for bin = find(nInBin ~= 0)
            nCompare = sum(nInBin(connectivity(:,bin)));
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
%                 Mask = abs(dr(:,:,d)) > halfBox(d);
%                 dr(Mask) = dr(Mask) - boxSize(d)*sign(dr(Mask));
                for ix = 1:length(p1)
                    Mask = abs(dr(ix,:,d)) > halfBox(d);
                    dr(ix,Mask,d) = dr(ix,Mask,d) - ...
                        sign(pos(p1(ix),d))*boxSize(d);
                end
            end
            % norm of relative distancces
            dd = sqrt(sum(dr.^2,3));
            dd(dd > rcut) = 0;
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
            acc(p1,:) = reshape(sum(Fr,2),size(acc(p1,:))) ...
                ./ repmat(mass(p1),1,2);
        end
    end
    %----------------------------------------------------------------------
    %----- End force ------------------------------------------------------
    %----------------------------------------------------------------------

    
    %----------------------------------------------------------------------
    %||||| Measurements |||||----------------------------------------------
    %----------------------------------------------------------------------
    function Measurements
        kin = zeros(Np,Dim);
        for d=1:Dim
            kin(:,d) = 0.5*mass.*vel(:,d).^2;
        end
        kinE = sum(sum(kin));
        potE = potentialE;
        totE = kinE +potE;
        fprintf(energyFileId, '%13.6E %13.6E %13.6E %13.6E\n', ...
            time, totE, kinE, potE);
    end
    %----------------------------------------------------------------------
    %----- End measurements -----------------------------------------------
    %----------------------------------------------------------------------

    
    %----------------------------------------------------------------------
    %||||| MovieWrite |||||------------------------------------------------
    %----------------------------------------------------------------------
    function MovieWrite
        fprintf(movieFileId, '%u\n', Np);
        fprintf(movieFileId, '%1s\n', ' ');
        for i=1:floor(Np/3)
            fprintf(movieFileId, '%3s  %13.6E %13.6E %13.6E %13.6E\n',...
                '1', pos(i,1), pos(i,2), vel(i,1), vel(i,2));
        end
        for i=floor(Np/3)+1:floor(2*Np/3)
            fprintf(movieFileId, '%3s  %13.6E %13.6E %13.6E %13.6E\n',...
                '2', pos(i,1), pos(i,2), vel(i,1), vel(i,2));
        end    
        for i=floor(2*Np/3)+1:Np
            fprintf(movieFileId, '%3s  %13.6E %13.6E %13.6E %13.6E\n',...
                '3', pos(i,1), pos(i,2), vel(i,1), vel(i,2));
        end
    end
    %----------------------------------------------------------------------
    %----- End movieWrite -----------------------------------------------
    %----------------------------------------------------------------------

    
    %----------------------------------------------------------------------
    %||||| SetNewInitialState |||||----------------------------------------
    %----------------------------------------------------------------------
    function SetNewInitialState
        % Initialize vectors
        density = zeros(nbSpacePts,1);
        velocity = zeros(nbSpacePts,1);
        temperature = zeros(nbSpacePts,1);
        
        % Compute densities and velocity moments in the CM frame
        for x=1:nbSpacePts
            density(x) = sum(distFunc(x,:));
            velocity(x) = distFunc(x,:)*velRange(:)./density(x);
            temperature(x) = distFunc(x,:)*(velRange(:)-velocity(x)).^2./density(x);
        end
        temperature = temperature*mass(1)/kB;
        
        % Determine the number of particles to put in each cell
        ratioOfPart = zeros(nbSpacePts,1);
        ratioOfPart = density.*double(Np);
        nbOfPart = round(ratioOfPart);
        nbOfPartTemp = nbOfPart;
        while(sum(nbOfPart) ~= Np)
            [val idx] = max(nbOfPartTemp);
            if(sum(nbOfPart) < Np)
                nbOfPartTemp(idx) = 0;
                nbOfPart(idx) = nbOfPart(idx) +1;
            elseif(sum(nbOfPart) > Np)
                nbOfPartTemp(idx) = 0;                
                nbOfPart(idx) = nbOfPart(idx) -1;                
            end
        end
        clear nbOfPartTemp;
        
        % Set new positions and velocities
        pos = rand(size(pos))-0.5;              % Uniform distribution
        pos(:,1) = pos(:,1).*cellSize(1);       % extended over a cell along x
        pos(:,2) = pos(:,2).*boxSize(2);        % extended over the box along y
        vel = randn(size(vel));                 % Normal velocity distribution
        partIdx = 1;
        for x=1:nbSpacePts
            increment = 0;
            while(increment < nbOfPart(x))
                pos(partIdx,1) = pos(partIdx,1) + cellCenters(x);
                vel(partIdx,1) = vel(partIdx,1).*sqrt(kB*temperature(x)./mass(1)) +velocity(x);
                vel(partIdx,2) = vel(partIdx,2).*sqrt(kB*temperature(x)./mass(1));
                partIdx = partIdx +1;
                increment = increment +1;
            end
        end
        
        % Display the generated distribution function
        finalDist = zeros(nbSpacePts,nbVelPts);
        posIdx = zeros(Np,1);
        velIdx = zeros(Np,1);
        
        posIdx = int64(pos(:,1)/cellSize(1)) +nbSpacePts/2;
        velIdx = int64(vel(:,1)/deltaVel)    +nbVelPts/2;
        
        for i=1:Np
            if(posIdx(i) > nbSpacePts) posIdx(i) = nbSpacePts; end
            if(posIdx(i) < 1         ) posIdx(i) = 1         ; end
            if(velIdx(i) > nbVelPts  ) velIdx(i) = nbVelPts  ; end
            if(velIdx(i) < 1         ) velIdx(i) = 1         ; end
            finalDist(posIdx(i),velIdx(i)) = finalDist(posIdx(i),velIdx(i)) +1.0;
        end
        finalDist = finalDist./sum(sum(finalDist));
%         figure('name', 'Generated initial distribution function');
%         imagesc(transpose(finalDist), [0 4.0e-3]);
        
    end
    %----------------------------------------------------------------------
    %----- End SetNewInitialState -----------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %||||| GetFinalDistribution |||||--------------------------------------
    %----------------------------------------------------------------------
    function GetFinalDistribution
        finalDist = zeros(nbSpacePts,nbVelPts);
        posIdx = zeros(Np,1);
        velIdx = zeros(Np,1);

        posIdx = int64(pos(:,1)/cellSize(1)) +nbSpacePts/2;
        velIdx = int64(vel(:,1)/deltaVel)    +nbVelPts/2;
        
        % if the velocity v is higher than the maximum value Vmax allowed by the
        % velocity grid, v contributes to Vmax.
        for i=1:Np
            if(posIdx(i) > nbSpacePts) posIdx(i) = nbSpacePts; end
            if(posIdx(i) < 1         ) posIdx(i) = 1         ; end
            if(velIdx(i) > nbVelPts  ) velIdx(i) = nbVelPts  ; end
            if(velIdx(i) < 1         ) velIdx(i) = 1         ; end
            
            finalDist(posIdx(i),velIdx(i)) = finalDist(posIdx(i),velIdx(i)) +1.0;
        end
        finalDist = finalDist./sum(sum(finalDist));
%         figure('name','Final distribution function');
%         imagesc(transpose(finalDist));
    end
    %----------------------------------------------------------------------
    %----- End GetFinalDistribution ---------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %||||| NNconnectivity |||||--------------------------------------------
    %----------------------------------------------------------------------
    function [Nneighbors, connectivity] = NNconnectivity
    %%% build nearest neighbor connectivity %%%
        % prevent unnecessary evaluations if lbin is factor of rcut
        fix = 1e-10 * (mod(rcut./lbin,1) == 0);
        % get max y-offset for each x-offset
        xLen = floor(rcut/lbin(1) + 1 - fix(1));
        yLen = floor(sqrt(rcut^2 - ((0:xLen-1)*lbin(1)).^2) / lbin(2) + ...
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
    %||||| BinParticles |||||----------------------------------------------
    %----------------------------------------------------------------------
    function BinParticles
        nInBin(:) = 0;
        binToParticles(:,:) = 0;
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
    
    
end
