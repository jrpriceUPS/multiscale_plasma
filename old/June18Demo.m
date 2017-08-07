%A demo of the functionality of the BGK_func, using interface and
%maxwellian. See documentation for these for more information.


clear all; close all;

%changing i will change which demo is computed
%
%i = 1 is a one species interface problem
%i = 2 is a two species interface problem
%i = 3 is a three species interface problem
%i = 4 is a one species sine function in density initial condition
%    NOTE: This one does not work yet due to our faulty poisson solver
i=4;

%timestep size, 1e-15 works for i=2,3, while 1e-14 is acceptable for i=1
%(the code makes this change automatically)
dt=1e-14;

%the number of timesteps, this can be changed as desired
%our steps are so small, the solution remains stable even for large Ndt (at
%least 10000)
Ndt=10000;

%if you are plotting, njumps is how often the results are plotted
njumps=10;

%legends is where you can name your particles, the default is to call them
%'Species 1', 'Species 2', etc. legends=[]; makes us have the default
legends=[];

%set this to 0 if you don't want the results plotted
plotting=1;

%set this to 1 if you want more detailed timing data
timing=0;

%converts eV to J
EV_TO_J = 1.60217653e-19;

%spatial domain data, you can change as you wish, but all my testing has
%been with these parameters
Lx = 800e-9;
Nx = 32;
dx = Lx/Nx;
Nv = 60;


if i==1
    %as discussed, one species can handle a larger timestep
    dt=1e-14;
    
    %below we set up an initial interface distribution
    
    %mass and charge of particle type
    mass=1.67377e-27;
    charge=1;
    
    %constant density (n), velocity (u), and temperature (T) on the left (L)
    %and right (R) sides of the interval
    %you may change any of these
    n_L=5.0e+28;
    u_L=[0.0 , 0.0];
    T_L=100*EV_TO_J;
    n_R=1e+28;
    u_R=[0.0 , 0.0];
    T_R=200*EV_TO_J;
    
    %this axis size is great for the current parameters
    %to allow the axis to change automatically while plotting, set
    %axissize=NaN
    axissize=[0 6e28
        0 400
        -5e4 5e4];
    
    
    
elseif i==2
    %below we set up an initial interface distribution
    
    %mass and charge of each particle type
    mass=[1.67377e-27 , 2*1.67377e-27];
    charge=[1.0 , 2.0];
    
    %constant density (n), velocity (u), and temperature (T) on the left (L)
    %and right (R) sides of the interval for each particle type
    %you may change any of these
    n_L=[5.0e+28 , 1e10];
    u_L=[0.0 , 0.0
        0.0 , 0.0];
    T_L=[100*EV_TO_J , 200*EV_TO_J];
    n_R=[1e+10 , 1.0e+28];
    u_R=[0.0 , 0.0
        0.0 , 0.0];
    T_R=[100*EV_TO_J , 200*EV_TO_J];
    
    %this axis size is great for the current parameters
    %to allow the axis to change automatically while plotting, set
    %axissize=NaN
    axissize=[  0  6e28
        0   500
        -3e5  3e5];
    
elseif i==3
    %below we set up an initial interface distribution
    
    %mass and charge of each particle type
    mass=[1.67377e-27 , 2*1.67377e-27, 3*1.67377e-27];
    charge=[1.0 , 2.0 , 3.0];
    
    %constant density (n), velocity (u), and temperature (T) on the left (L)
    %and right (R) sides of the interval for each particle type
    %you may change any of these
    n_L=[5.0e+28 , 1e10, 1e18];
    u_L=[0.0 , 0.0
        0.0 , 0.0
        0.0 , 0.0];
    T_L=[100*EV_TO_J , 200*EV_TO_J , 300*EV_TO_J];
    n_R=[1e+10 , 1.0e+28 , 5e+28];
    u_R=[0.0 , 0.0
        0.0 , 0.0
        0.0 , 0.0];
    T_R=[400*EV_TO_J , 500*EV_TO_J , 600*EV_TO_J];
    
    %this axis size is great for the current parameters
    %to allow the axis to change automatically while plotting, set
    %axissize=NaN
    axissize=[0 6e28
        0 1400
        -4e5 4e5];
    
    
elseif i==4
    %An initial distribution with one particle that has constant
    %temperature and initial velocity, but sinusoidal density with
    %amplitude A
    
    dt=1e-14;
    
    %mass and charge of particle
    mass=1.67377e-27;
    charge=1;
    
    %A = amplitude of density sine wave can be changed as desired, make
    %sure it does not lead to negative initial densities
    A=0.2*1e28;
    
    %x grid
    x = 0.5*dx:dx:(Lx-0.5*dx);
    
    %density assigned as sine wave periodic on domain with amplitude A
    density=A*sin(2*pi*x/Lx)+A*10;
    
    %constant velocity and temperature
    u=[0 , 0];
    T=100*EV_TO_J;
    
    %calculate maximum thermal speed and maximum velocity we consider
    therm_speed_max = sqrt(T/mass);
    Lv = 10.0*therm_speed_max;
    
    %set the parameters to make our maxwellian
    field1 = 'vrange'; value1 = linspace(-Lv,Lv,Nv);
    field2 = 'density'; value2 = density;
    field3 = 'temperature'; value3 = squeeze(repmat(T,[1,1,Nx]));
    field4 = 'velocity'; value4 = permute(repmat(u,[1,1,Nx]),[1,3,2]);
    field5 = 'mass'; value5 = mass;
    s1 = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
    
    %initialize initial distribution
    f0=zeros(1,Nx,Nv,Nv);
    
    %find maxwellian
    f0fill=maxwellian(s1,[],Nx,Nv,1);
    
    %fill initial distribution
    f0(1,:,:,:)=f0fill;
    
    %this axis size is great for the current parameters
    %to allow the axis to change automatically while plotting, set
    %axissize=NaN
    axissize=[  0  3e28
                0   150
              -2e4  2e4];
end

%fill spatial parameters structure

field1 = 'Nx'; value1 = Nx;
field2 = 'Nv'; value2 = Nv;
field3 = 'timestep'; value3 = dt;
field4 = 'timestepnum'; value4 = Ndt;
field5 = 'Lx'; value5 = Lx;
field6 = 'Lv'; value6 = 0;
field7 = 'njumps'; value7 = njumps;
field8 = 'plotting'; value8 = plotting;
field9 = 'timing'; value9 = timing;
field10 = 'axissize'; value10 = axissize;
field11 = 'legend'; value11=legends;

spatialparams = struct(field1,value1,field2,value2,field3,value3,...
    field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,...
    field9,value9,field10,value10,field11,value11);

%if we are using an interface problem, construct the maxwellians using
%given parameters
if i<4
    
    field1 = 'density';     value1 = n_L;
    field2 = 'temperature'; value2 = T_L;
    field3 = 'velocity';    value3 = u_L;
    
    leftparams = struct(field1,value1,field2,value2,field3,value3);
    
    field1 = 'density';     value1 = n_R;
    field2 = 'temperature'; value2 = T_R;
    field3 = 'velocity';    value3 = u_R;
    
    rightparams = struct(field1,value1,field2,value2,field3,value3);
    
    [f0,Lv]=interface(mass,spatialparams,leftparams,rightparams);
    
end

spatialparams.Lv = Lv;

%run our BGK_func on the given initial distribution

[f,n,u,T]=BGK_func(mass,charge,f0,spatialparams);



