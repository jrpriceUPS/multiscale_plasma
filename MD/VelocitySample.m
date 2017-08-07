function vel = VelocitySample(f,Np,Lv,factor,Nmoment,abserr,relerr)
%%%---------------------------------------------------------------------%%%
%
% Function to randomly sample velocities from a distribution function. Uses
% linear interpolation to refine the grid and then creates a CDF from which
% to sample the 2D velocities.
%
% INPUTS:
%
%   f       -   the distribution function to sample (1 species, 1 cell),
%               must be an Nv x Nv array.
%   Np      -   number of particles to assign velocities
%   Lv      -   maximum (absolute value) velocity represented by f
%   factor  -   refinement factor: results in 2^factor points between each
%               point on the original f surface.
%   Nmoment -   number of moments to match (2-4)
%   abserr  -   absolute error tolerance in the moments (moments x 1)
%   relerr  -   relative error tolerance in the moments (moments x 1)
%
% OUTPUTS:
%
%   vel     -   2D velocity for the particles (Np x 2 array)
%
%%%---------------------------------------------------------------------%%%
if nargin < 5
    Nmoment = 2;
end
if nargin < 6
    abserr  = 0.05 * ones(Nmoment,1);
    relerr  = 0.05 * ones(Nmoment,1);
end
if Nmoment > 4
    warning('Only first 4 moments supported. Using 4 moments.')
    Nmoment = 4;
end

% linearly interpolate f
if factor > 0
    f   = interp2(f,factor);
end

% set up velocity arrays
Nv  = size(f,1);
v   = linspace(-Lv,Lv,Nv)';
vel = zeros(Np,2);

% set up moment variables
nDist   =   0;
uDist   =   zeros(1,2);
TDist   =   0;
qDist   =   zeros(1,2);
M4Dist  =   0;
getDistributionMoments

uPart   =   zeros(1,2);
TPart   =   0;
qPart   =   zeros(1,2);
M4Part  =   0;

% create cumulative distribution
distributionX = cumsum(sum(f,2)) / sum(f(:));
distributionY = cumsum(f,2) ./ repmat(sum(f,2),[1,Nv]);

% sample the distribution for each particle while error too big
success = 0;
count = 0;
while success < Nmoment
    success = 0;
    count = count + 1;
    
    if mod(count,100)==0
        fprintf('sampling distribution, attempt %d\n',count)
    end
    
    if Np==0
        print('oh no!')
    end
    for p = 1:Np
        ix = find(distributionX >= rand(1), 1, 'first');
        jx = find(distributionY(ix,:) >= rand(1), 1, 'first');
        vel(p,:) = [v(ix) v(jx)];
    end
    getParticleMoments
    
    % check moments
    success = success + ...
        (sum(abs(uDist-uPart)./abs(uDist) < relerr(1)) == 2 || ...
        sum(abs(uDist-uPart) < abserr(1)) == 2);
    success = success + (abs(TDist-TPart)/abs(TDist) < relerr(2) || ...
        abs(TDist-TPart) < abserr(2));
    if Nmoment >= 3
        success = success + ...
            (sum(abs(qDist-qPart)./abs(qDist) < relerr(3)) == 2 || ...
            sum(abs(qDist-qPart) < abserr(3)) == 2);
    end
    if Nmoment >= 4
        success = success + (abs(M4Dist-M4Part)/abs(M4Dist) < relerr(4) || ...
            abs(M4Dist-M4Part) < abserr(4));
    end
end

% nested functions
    function getDistributionMoments
        dv          =   v(2) - v(1);
        wtN         =   dv*ones(Nv,1);
        wtN(1)      =   0.5*wtN(1);
        wtN(end)    =   0.5*wtN(end);
        W           =   wtN*wtN';
        v1          =   repmat(v,[1,Nv]);
        v2          =   v1';
        
        nDist       =   sum(sum(W .* f,2),1);
        uDist(1)    =   sum(sum(W .* v1 .* f,2),1) / nDist;
        uDist(2)    =   sum(sum(W .* v2 .* f,2),1) / nDist;
        vmu1        =   v1 - uDist(1);
        vmu2        =   v2 - uDist(2);
        TDist       =   sum(sum(W .* (vmu1.^2 + vmu2.^2) .* f,2),1) / ...
            (2*nDist);
        qDist(1)    =   1/2 * sum(sum(W .* (vmu1.^2 + vmu2.^2) .* ...
            vmu1 .* f,2),1) / nDist;
        qDist(2)    =   1/2 * sum(sum(W .* (vmu1.^2 + vmu2.^2) .* ...
            vmu2 .* f,2),1) / nDist;
        M4Dist      =   sum(sum(W .* (vmu1.^2 + vmu2.^2).^2 .* f,2),1) /...
            nDist;
        
    end

    function getParticleMoments
        uPart   =   sum(vel,1) / Np;
        vmu1    =   vel(:,1) - uPart(1);
        vmu2    =   vel(:,2) - uPart(2);
        TPart   =   sum(vmu1.^2 + vmu2.^2) / (2*Np);
        qPart   =   1/2 * sum(repmat(vmu1.^2 + vmu2.^2,[1,2]) .* ...
            repmat(uPart,[Np,1]),1) / Np;
        M4Part  =   sum((vmu1.^2 + vmu2.^2).^2) / Np;
    end
end