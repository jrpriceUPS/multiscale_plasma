function [f,Lv]=interface(mass,spatialparams,leftparams,rightparams)
%
%[f,Lv]=interface(mass,spatialparams,leftparams,rightparams)
%
%A function to set up an initial condition with an interface. The
%densities, temperatures, and bulk velocities of each particle type is
%constant on each side. We select maxwellians for each side that have the
%proper densities, temperatures, and velocities.
%
%
%%%%%%%%
%INPUTS%
%%%%%%%%
%
%mass is a 1xn row vector of the masses of each species
%
%spatialparams contains information about our domain and other parameters
%if any elements are empty, they are filled with default values.
%
%spatialparams.Nx is the number of spatial points in our discretization
%                        default value = 32
%
%spatialparams.Nv is the number of velocity points in our discretization
%                        default value = 60
%
%
%leftparams and rightparams contain information about the distributions on
%each side
%
%leftparams.density is an nSp x Nx/2 matrix of the density of each particle
%                   at each x position left of the interval
%
%leftparams.temperature is an nSp x Nx/2 matrix of the temperature of each
%                       particle at each x position left of the interval
%
%leftparams.velocity is an nSp x Nx/2 x 2 matrix of the velocity of each
%                    particle at each x position left of the interval
%
%
%rightparams contains similar information for the right side of the interval
%
%
%%%%%%%%%
%OUTPUTS%
%%%%%%%%%
%
%f is a nSp x Nx x Nv x Nv probability distribution such that each particle
%is maxwellian in space and has the correct density, temperature, and bulk
%velocity on each side of the interval
%
%Lv is the largest velocity that will be used in probability calculations
%(the cutoff velocity). It is selected based on the highest initial
%temperature
%
%
%%%%%%%%%%%%%%%
%EXAMPLE INPUT%
%%%%%%%%%%%%%%%
%
%Copy and paste this to see an example of the function's behavior:
%
%     EV_TO_J = 1.60217653e-19;
%
%     mass=[1.67377e-27 , 2*1.67377e-27];
%
%     n_L=[5.0e+28 , 1e10];
%     u_L=[100.0 , 200.0
%         300.0 , 400.0];     
%     T_L=[100*EV_TO_J , 200*EV_TO_J];
%
%     n_R=[1e+10 , 1.0e+28];
%     u_R=[500.0 , 600.0
%         700.0 , 800.0];
%     T_R=[100*EV_TO_J , 200*EV_TO_J];
%
%     spatialparams = struct('Nx',32,'Nv',60);
%     leftparams = struct('density',n_L,'temperature',T_L,'velocity',u_L);
%     rightparams = struct('density',n_R,'temperature',T_R,'velocity',u_R);
%
%     [f,Lv]=interface(mass,spatialparams,leftparams,rightparams);


%extract parameters from input structures
Nx=spatialparams.Nx;
Nv=spatialparams.Nv;
n_L=leftparams.density;
T_L=leftparams.temperature;
u_L=rightparams.velocity;
n_R=rightparams.density;
T_R=rightparams.temperature;
u_R=rightparams.velocity;
nSp=length(mass);

%compute maximum thermal speed, and use this to select the maximum velocity
%we consider
therm_speed_max = sqrt(max([T_L./mass T_R./mass]));
Lv = 10.0*therm_speed_max;

%generate the data structures for the left and right sides of the interval
%that will be used to construct the maxwellian on either side

%left side
field1 = 'vrange'; value1 = linspace(-Lv,Lv,Nv);
field2 = 'density'; value2 = squeeze(repmat(n_L,[1,1,Nx/2]));
if length(mass)==1
    value2=value2.';
end
field3 = 'temperature'; value3 = squeeze(repmat(T_L,[1,1,Nx/2]));
field4 = 'velocity'; value4 = permute(repmat(u_L,[1,1,Nx/2]),[1,3,2]);
field5 = 'mass'; value5 = mass;
sL = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);

%right side
field1 = 'vrange'; value1 = linspace(-Lv,Lv,Nv);
field2 = 'density'; value2 = squeeze(repmat(n_R,[1,1,Nx/2]));
if length(mass)==1
    value2=value2.';
end
field3 = 'temperature'; value3 = squeeze(repmat(T_R,[1,1,Nx/2]));
field4 = 'velocity'; value4 = permute(repmat(u_R,[1,1,Nx/2]),[1,3,2]);
field5 = 'mass'; value5 = mass;
sR = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);

%compute the maxwellians
fL=maxwellian(sL,[],Nx/2,Nv,nSp);
fR=maxwellian(sR,[],Nx/2,Nv,nSp);

%fill the output array with the proper maxwellian for each side of the
%interval
if length(mass)==1
    
    f=zeros(1,Nx,Nv,Nv);
    f(1,1:Nx/2,:,:)=fL;
    f(1,1+Nx/2:Nx,:,:)=fR;
    
else
    
    f=zeros(length(mass),Nx,Nv,Nv);
    f(:,1:Nx/2,:,:)=fL;
    f(:,1+Nx/2:Nx,:,:)=fR;
    
end