% clear variables and set up paths
%clear all, close all
addpath('../BGK/')
addpath('../MD/')
addpath('../dependencies/')


i=2;


if i==1
    %All the below is just input data for a BGK distribution that has
    %notable density, bulk velocity, and temperature (to make sure
    %equilibrate step is working)
    lambda=1.0e-2*5.2917721092e-11; %1/1000 the size of the domain
    Ndt=0; %this will give us just the density, temperature, velocity, and initial distribution from BGK
    njumps=1; %record data ever 10 steps
    legends=[]; %no legends BGK
    plotting=0; %no plotting BGK
    timing=0; %no timing BGK
    EV_TO_J = 1.60217653e-19; %conversion
    Lx = lambda*1000; %domain size
    Nx = 20; %number of cells
    dx = Lx/Nx;
    Nv = 60; %number of velocity points
    
    
    dt=1e-15; %timestep for BGK
    mass=1.67377e-27; %mass of particle
    charge=1; %charge of particle
    
    %sinusoidal densities, temperatures, and velocities, to make sure
    %equilibrate is working correctly
    A=0.2*1e28;
    x = 0.5*dx:dx:(Lx-0.5*dx);
    density=A*sin(2*pi*x/Lx)+A*10;
    u=[A*cos(2*pi*x/Lx) ; A*sin(2*pi*x/Lx)].'*10^-27;
    %T=10*EV_TO_J*cos(2*pi*x/Lx)+100*EV_TO_J;
    
    
    %below is a step function for temperature, it is handled roughly the
    %same way as the other one
     T=zeros(1,Nx);
     T(1:Nx/2)=100*EV_TO_J;
     T(Nx/2+1:Nx)=200*EV_TO_J;
    
    %other inputs
    therm_speed_max = max(sqrt(T/mass));
    Lv = 10.0*therm_speed_max;
    v = linspace(-Lv,Lv,Nv);
    dv = v(2) - v(1);
    
    %set up struct for maxwellian
    species.Nsp = 1;
    species.mass = mass;
    species.charge = charge;
    moments.n = density;
    moments.u = u';
    moments.T = T;
    
    %BGK parameters
    params.Nx = Nx;
    params.Nv = Nv;
    params.dt = dt;
    params.Ndt = Ndt;
    params.Lx = Lx;
    params.Lv = Lv;
    params.Njumps = njumps;
    params.plotting = plotting;
    params.timing = timing;
    params.axissize = [];
    params.legend = legends;    
    
    %run maxwellian
    f0=zeros(1,Nx,Nv,Nv);
    f0fill=maxwellian([],species,moments,params);
    f0(1,:,:,:)=f0fill;
    
    %BGK func run
    [f,n,u,T]=BGK_Func(mass,charge,f0,spatialparams);
    
    
    
    %number of particles
    N = 1000;
    
    %aspect ratio
    AR=1;
    
    %timestep for MD
    dt = 1e-16;%2.418884326505e-22;
    
    %number of timesteps
    Ndt=500;
    
    %how often to record data
    Njumps=1;
    
    %cutoff for force calculation
    cutoff=4*lambda;
    
    %minimum particle spacing
    spacing=1;
    
    %number of subdivisions in bins
    Ndiv=1;
    
    %amount of timing data to display
    timing=2;
    
    %safety factor for binning
    safetyfactor=10;
    
    %equilibration constants
    gamma=1;
    dteq=1e-16;
    Ndteq=1000;
    
    
    
elseif i==2
    
    %haven't tried to do a 2-particle type demo yet
    
   
    mass=[1.67377e-27 , 1000*1.67377e-27];
    charge=[1.0 , 1.0];
    
    %All the below is just input data for a BGK distribution that has
    %notable density, bulk velocity, and temperature (to make sure
    %equilibrate step is working)
    lambda=1.0e-2*5.2917721092e-11; %1/1000 the size of the domain
    Ndt=0; %this will give us just the density, temperature, velocity, and initial distribution from BGK
    njumps=10; %record data ever 10 steps
    legends=[]; %no legends BGK
    plotting=0; %no plotting BGK
    timing=0; %no timing BGK
    EV_TO_J = 1.60217653e-19; %conversion
    Lx = lambda*1000; %domain size
    Nx = 20; %number of cells
    dx = Lx/Nx;
    Nv = 60; %number of velocity points
    
    
    dt=1e-15; %timestep for BGK
    
    %sinusoidal densities, temperatures, and velocities, to make sure
    %equilibrate is working correctly
    A=0.2*1e28;
    x = 0.5*dx:dx:(Lx-0.5*dx);
    density1=A*sin(2*pi*x/Lx)+A*10;
    density2=A*cos(2*pi*x/Lx)+A*10;
    density=[density1;density2];
    
    u1=[A*cos(2*pi*x/Lx) ; A*sin(2*pi*x/Lx)].'*10^-27;
    
    u2=[A*sin(2*pi*x/Lx) ; A*cos(2*pi*x/Lx)].'*10^-27;
    u=zeros(2,Nx,2);
    
    u(1,:,:)=u1;
    u(2,:,:)=u2;
    u(:) = 0;
    
    T1=10*EV_TO_J*cos(2*pi*x/Lx)+1000*EV_TO_J;
    T2=10*EV_TO_J*sin(2*pi*x/Lx)+1000*EV_TO_J;
    
    T=[T1;T2];
    %below is a step function for temperature, it is handled roughly the
    %same way as the other one
     %T=zeros(1,Nx);
     %T(1:Nx/2)=100*EV_TO_J;
     %T(Nx/2+1:Nx)=200*EV_TO_J;
    
    %other inputs
    therm_speed_max = max(max(sqrt(T./squeeze(repmat(mass,[1,1,Nx])))));
    Lv = 10.0*therm_speed_max;
    v = linspace(-Lv,Lv,Nv);
    dv = v(2) - v(1);
    
    %set up struct for maxwellian
    species.Nsp = 2;
    species.mass = mass';
    species.charge = charge';
    moments.n = density;
    moments.u = u;
    moments.T = T;
    
    %BGK parameters
    params.Nx = Nx;
    params.Nv = Nv;
    params.dt = dt;
    params.Ndt = Ndt;
    params.Lx = Lx;
    params.Lv = Lv;
    params.Njumps = njumps;
    params.plotting = plotting;
    params.timing = timing;
    params.axissize = [];
    params.legends = {'Species 1','Species 2'};
    
    %run maxwellian
    f0=maxwellian([],species,moments,params);
    
    %BGK func run
    [f,n,u,T]=BGK_Func(species,f0,params);
    
    
    
    %number of particles
    N = [250 , 500];
    
    %aspect ratio
    AR=1;
    
    %timestep for MD
    dt=4e-12;
    
    %number of timesteps
    Ndt=500;
    
    %how often to record data
    Njumps=1;
    
    %cutoff for force calculation
    cutoff=4*lambda;
    
    %minimum particle spacing
    spacing=1;
    
    %number of subdivisions in bins
    Ndiv=1;
    
    %amount of timing data to display
    timing=2;
    
    %safety factor for binning
    safetyfactor=10;
    
    %equilibration constants
    gamma=1;
    dteq=4e-12;
    Ndteq=0;
    
    
end

%set up species struct
species.N = N;
species.mass = mass;
species.charge = charge;

%set up distribution struct
distribution.f = f;
distribution.n = n;
distribution.u = u;
distribution.T = T;

%set up domain struct
domain.Lx = Lx;
domain.AR = AR;
domain.lambda = lambda;
domain.dt = dt;
domain.Ndt = Ndt;
domain.Nx = Nx;
domain.Nv = Nv;
domain.dv = dv;

%set up params struct
params.cutoff = cutoff;
params.spacing = spacing;
params.Ndiv = Ndiv;
params.Njumps = Njumps;
params.timing = timing;
params.safetyfactor = safetyfactor;

%equillibration params
equil.gamma = gamma;
equil.dteq = dteq;
equil.Ndteq = Ndteq;

%run MD
MD_Func(species,distribution,domain,params,equil);