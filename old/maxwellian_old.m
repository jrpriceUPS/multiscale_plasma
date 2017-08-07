function fmax=maxwellian(baseparams,repeatedparams,Nx,Nv,nSp)
%
%fmax=maxwellian(baseparams,repeatedparams,Nx,Nv,nSp)
%
%Creates a maxwellian distribution array with dimension nSp x Nx x Nv x Nv
%using data from either baseparams or repeatedparams. Repeatedparams is the
%structure used in the actual computation, but if it is not provided, it
%can be computed from the easier to understand baseparams.
%
%%%%%%%%
%INPUTS%
%%%%%%%%
%
%baseparams contains information about our desired distribution
%
%baseparams.vrange is the range of velocity values vmin:dv:vmax
%
%baseparams.density is a nSp x Nx array of the density of each particle
%                       type at each x position
%
%baseparams.temperature is a nSp x Nx array of the temperature of each
%                       particle type at each x position
%
%baseparams.velocity is a nSp x Nx x 2 array of the velocity of each
%                       particle at each x position
%
%baseparams.mass is the a 1 x nSp row vector of the mass of each particle
%                       type
%
%
%repeatedparams is a similar structure with the elements of baseparams
%repeated to properly execute the final computation
%
%repeatedparams.V1 is a nSp x Nx x Nv x Nv array where each nSp x Nx
%                       position contains the same Nv x Nv array of vrange
%                       repeated Nv times
%
%repeatedparams.V2 is a nSp x Nx x Nv x Nv array where each nSp x Nx
%                       position contains the same Nv x Nv array of vrange.'
%                       repeated Nv times
%
%repeatedparams.mass is a nSp x Nx x Nv x Nv array where each entry is the
%                       mass corresponding to the correct species
%
%repeatedparams.density is a nSp x Nx x Nv x Nv array where each nSp x Nx
%                       position contains the density
%
%repeatedparams.temperature is a nSp x Nx x Nv x Nv array where each nSp x Nx
%                       position contains the temperature
%
%repeatedparams.velocity is a nSp x 2 xNx x Nv x Nv array where each
%                       nSp x Nx x 2 position contains the velocity
%
%in cases where there is only one species, all these quantities have the
%nSp=1 dimension squeezed away
%
%
%%%%%%%%%
%OUTPUTS%
%%%%%%%%%
%
%fmax is a nSp x Nx x Nv x Nv maxwellian distribution with the desired
%densities, temperatures, and velocities for each species type
%
%
%%%%%%%%%%%%%%%
%EXAMPLE INPUT%
%%%%%%%%%%%%%%%
%
%Copy and paste this to see an example of the function's behavior:
%
%  Nx=32; Nv=60; Lx=800e-9; Lv=2e-6; dx=Lx/Nx;
%  mass=[1.67377e-27 , 2*1.67377e-27];
%
%  density=zeros(2,Nx);
%
%  A = 0.2*1e28;
%
%  density(1,:) = A*sin(2*pi*(0.5*dx:dx:(Lx-0.5*dx))/Lx)+A*10;
%  density(2,1:Nx/2) = 3e28*ones(1,Nx/2);density(2,1+Nx/2:Nx)=1e28*ones(1,Nx/2);
%
%  T1 = 100*1.60217653e-19*sin(2*pi*(0.5*dx:dx:(Lx-0.5*dx))/Lx)+200*1.60217653e-19;
%  T2 = 200*1.60217653e-19*ones(1,Nx); T = [T1;T2];
%
%  u = [1.0 , 2.0;3.0 , 4.0];
%
%  field1 = 'vrange'; value1 = linspace(-Lv,Lv,Nv);
%  field2 = 'density'; value2 = density;
%  field3 = 'temperature'; value3 = T;
%  field4 = 'velocity'; value4 = permute(repmat(u,[1,1,Nx]),[1,3,2]);
%  field5 = 'mass'; value5 = mass;
%
%  maxwellstruct=struct(field1,value1,field2,value2,field3,value3,...
%                       field4,value4,field5,value5);
%
%  f=maxwellian(maxwellstruct,[],Nx,Nv,2);

%if no repeatedparams structure was provided, construct one to fill from
%baseparams
if isempty(repeatedparams)
    repeatedparams=struct('V1',[],'V2',[],'mass',[],'density',[],'temperature',[],'velocity',[]);
end

%if there are ANY empty elements of repeated params, use baseparams to fill
%the missing ones according to the definitions in the header
if isempty(repeatedparams.V1)+isempty(repeatedparams.V2)+isempty(repeatedparams.mass)+isempty(repeatedparams.density)+isempty(repeatedparams.temperature)+isempty(repeatedparams.velocity)~=0
    
    if isempty(repeatedparams.mass)
        if isempty(baseparams.mass)
            disp('No mass provided')
            return;
        end
        repeatedparams.mass=squeeze(repmat(baseparams.mass,[1,1,Nx,Nv,Nv]));
    end
    if isempty(repeatedparams.V1)
        if isempty(baseparams.vrange)
            disp('No velocity range provided')
            return;
        end
        V1=squeeze(repmat(baseparams.vrange,[1,1,Nv]));
        if nSp==1
            repeatedparams.V1=permute(repmat(V1,[1,1,Nx]),[3,1,2]);
        else
            repeatedparams.V1 = permute(repmat(V1,[1,1,nSp,Nx]),[3,4,1,2]);
        end
    end
    if isempty(repeatedparams.V2)
        if isempty(baseparams.vrange)
            disp('No velocity range provided')
            return;
        end
        V2=squeeze(repmat(baseparams.vrange,[1,1,Nv])).';
        if nSp==1
            repeatedparams.V2=permute(repmat(V2,[1,1,Nx]),[3,1,2]);
        else
            repeatedparams.V2 = permute(repmat(V2,[1,1,nSp,Nx]),[3,4,1,2]);
        end
    end
    if isempty(repeatedparams.density)
        if isempty(baseparams.density)
            disp('No density provided')
            return;
        end
        if nSp==1
            repeatedparams.density=squeeze(repmat(baseparams.density,[1,1,Nv,Nv]));
        else
            repeatedparams.density=repmat(baseparams.density,[1,1,Nv,Nv]);
        end
    end
    if isempty(repeatedparams.velocity)
        if isempty(baseparams.velocity)
            disp('No velocity provided')
            return;
        end
        if nSp==1
            repeatedparams.velocity=squeeze(repmat(baseparams.velocity,[1,1,1,Nv,Nv]));
        else
            repeatedparams.velocity=repmat(baseparams.velocity,[1,1,1,Nv,Nv]);
        end
    end
    if isempty(repeatedparams.temperature)
        if isempty(baseparams.temperature)
            disp('No temperature provided')
            return;
        end
        repeatedparams.temperature=squeeze(repmat(baseparams.temperature,[1,1,Nv,Nv]));
    end
end

%grab the variables from repeatedparams
V1r=repeatedparams.V1;
V2r=repeatedparams.V2;
ur=repeatedparams.velocity;
nr=repeatedparams.density;
massr=repeatedparams.mass;
Tr=repeatedparams.temperature;

%compute the distance between the given velocities and the velocity grids
if nSp==1   
    V1a=V1r-squeeze(ur(:,1,:,:));
    V2b=V2r-squeeze(ur(:,2,:,:));  
else  
    V1a=V1r-squeeze(ur(:,:,1,:,:));
    V2b=V2r-squeeze(ur(:,:,2,:,:));
end

%construct the maxwellian
fmax=nr.*massr./(2*pi*Tr).*exp(-.5*massr.*...
    (V1a.^2+V2b.^2)./Tr);