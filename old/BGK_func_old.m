function [f,density,u,temperature]=BGK_func(mass,charge,f0,spatialparams)
%
%[f,density,temperature,u]=BGK_func(mass,charge,f0,spatialparams)
%
%A function to solve the 1D-2V BGK equation for the probability
%distributions of n species of particles on a one-dimensional periodic
%domain with two dimensional velocities.
%
%%%%%%%%
%INPUTS%
%%%%%%%%
%
%mass is a 1xn row vector of the masses of each species              [kg]
%charge is a 1xn row vector of the charge level of each species      [e_0]
%f0 is an n x Nx x Nv x Nv array of the initial probability distribution
%        i.e. f0(i,x,v1,v2) is the probability that a particle of species i
%             is at position x with x velocity v1 and y velocity v2
%
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
%spatialparams.dt is the time step (seconds) of our integration, it must be
%                 very small
%                        default value = 1e-15
%
%spatialparams.Ndt is the number of time steps we will take
%                        default value = 1000
%
%spatialparams.Lx is the width of our spatial domain (m)
%                        default value = 800e-9
%
%spatialparams.Lv is the maximum velocity we integrate [m/s]
%                        default value = 2e6
%
%spatialparams.Njumps is how often (in timesteps) we plot our results if
%                     we are plotting
%                        default value = 1
%
%spatialparams.Nx is the number of spatial points in our discretization
%                        default value = 32
%
%spatialparams.plotting is a logical indicating whether we plot results
%                        default value = 0
%
%spatialparams.timing is a logical indicating whether detailed timing is
%                     undergone. If 0, the full time is printed, if 1, the
%                     average time of each computation is printed
%                        default value = 32
%
%spatialparams.axissize is a matrix of the y-axis for the density plots
%                       (row 1), the temperature plots (row 2),
%                       and the velocity plots (row 3). If this is NaN, no
%                       axis is specified
%                       default value = NaN
%
%spatialparams.legends is a cell array of the string names of each species
%                      used for plotting purposes
%                       default value = {'Species 1', 'Species 2',...}
%
%
%%%%%%%%%
%OUTPUTS%
%%%%%%%%%
%
%f is the final probability distribution with the same shape and
%  interpretation as f0
%
%density is the n x Nx particle density of this distribution [m^-3]
%
%temperature is the n x Nx temperature of this distribution  [kT in J]
%
%u is the n x Nx x 2 bulk velocities of this distribution    [m/s]
%
%
%%%%%%%%%%%%%%%
%EXAMPLE INPUT%
%%%%%%%%%%%%%%%
%
%Copy and paste this to see an example of the function's behavior:
%
%   field1  = 'Nx'; value1 = 32;
%   field2  = 'Nv'; value2 = 60;
%   field3  = 'timestep'; value3 = 1e-15;
%   field4  = 'timestepnum'; value4 = 1000;
%   field5  = 'Lx'; value5 = 800e-9;
%   field6  = 'Lv'; value6 = 2e6;
%   field7  = 'njumps'; value7 = 10;
%   field8 = 'plotting'; value8 = 1;
%   field9 = 'timing'; value9 = 0;
%   field10 = 'axissize'; value10 = [0 4e28; 0 400; -5e4 5e4];
%   field11 = 'legend'; value11=[];
%
%  spatialparams = struct(field1,value1,field2,value2,field3,value3,...
%      field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,...
%      field9,value9,field10,value10,field11,value11);
%
%  mass   = [1.67377e-27 , 2*1.67377e-27];
%  charge = [1.0 , 2.0];
%
%  density=zeros(2,spatialparams.Nx);
%
%  A = 0.2*1e28;
%
%  Lx = spatialparams.Lx; Nx=spatialparams.Nx; dx=Lx/Nx;
%  Lv = spatialparams.Lv; Nv=spatialparams.Nv;
%
%  density(1,:) = A*sin(2*pi*(0.5*dx:dx:(Lx-0.5*dx))/Lx)+A*10;
%
%  density(2,1:Nx/2) = 3e28*ones(1,Nx/2);density(2,1+Nx/2:Nx)=1e28*ones(1,Nx/2);
%
%  T1 = 100*1.60217653e-19*sin(2*pi*(0.5*dx:dx:(Lx-0.5*dx))/Lx)+200*1.60217653e-19;
%  T2 = 200*1.60217653e-19*ones(1,Nx); T = [T1;T2];
%
%  u = [1.0 , 2.0;3.0 , 4.0];
%
%  field14 = 'vrange'; value14 = linspace(-Lv,Lv,Nv);
%  field15 = 'density'; value15 = density;
%  field16 = 'temperature'; value16 = T;
%  field17 = 'velocity'; value17 = permute(repmat(u,[1,1,Nx]),[1,3,2]);
%  field18 = 'mass'; value18 = mass;
%
%  maxwellstruct=struct(field14,value14,field15,value15,field16,value16,...
%                       field17,value17,field18,value18);
%
%  f0=maxwellian(maxwellstruct,[],Nx,Nv,2);
%
%  [f, density, temperature, u] = BGK_func(mass,charge,f0,spatialparams);


%begin recording time
tic;
start = toc;



%%%%%%%%%%%%%%%%
%INITIALIZATION%
%%%%%%%%%%%%%%%%

%if not enough inputs are provided, quit the function

if nargin<2
    disp('Not enough inputs');
    return
end

%if spatialparams is not provided, create an empty one to be filled with
%default values

if nargin<3
    f1 = 'Nx';f2 = 'Nv';f3 = 'timestep';f4 = 'timestepnum';f5 = 'Lx';f6 = 'Lv';
    f7 = 'njumps';f8 = 'plotting';f9 = 'timing';
    f10 = 'axissize';f11 = 'legend';
    spatialparams=struct(f1,[],f2,[],f3,[],f4,[],f5,[],f6,[],f7,[],f8,[],f9,[],...
        f10,[],f11,[]);
end


%set parameters (see header for explanation of each)

Nx = spatialparams.Nx;
Nv = spatialparams.Nv;
dt = spatialparams.timestep;
Ndt = spatialparams.timestepnum;
Lx = spatialparams.Lx;
Lv = spatialparams.Lv;
Njumps = spatialparams.njumps;
plotting = spatialparams.plotting;
timing = spatialparams.timing;
axissize = spatialparams.axissize;
legends = spatialparams.legend;



%set a color for each particle type such that they are visually distinct
%from one another

colors=distinguishable_colors(length(mass));


%initialize marker for assigning default values

flag=0;


%if any parameters were not filled, assign default values and set flag to 1

if isempty(Nx)
    Nx=32;
    flag=1;
end

if isempty(Nv)
    Nv=60;
    flag=1;
end

if isempty(dt)
    dt=1e-15;
    flag=1;
end

if isempty(Ndt)
    Ndt=1000;
    flag=1;
end

if isempty(Lx)
    Lx=800e-9;
    flag=1;
end

if isempty(Lv)
    Lv=2e6;
    flag=1;
end

if isempty(Njumps)
    Njumps=1;
    flag=1;
end

if isempty(plotting)
    plotting=0;
    flag=1;
end

if isempty(timing)
    timing=0;
    flag=1;
end

if isempty(axissize)
    axissize=NaN;
    flag=1;
end

if isempty(legends)
    for i=1:length(mass)
        legends{i}=sprintf('Species %d',i);
    end
    flag=1;
end


%Inform the user that default values were used

if flag==1
    disp('WARNING: Default values for some variables were used!')
end


%if the number of species, nSp, is different between the different inputs,
%return

if (length(mass)~=length(charge)|| length(mass)~=size(f0,1))
    disp('Number of species unclear')
    return;
end


%initialize plot windows

if plotting==1
    figure('OuterPosition',[1,385,560,493]);
    figure('OuterPosition',[881,385,560,493]);
    figure('OuterPosition',[440,49,560,363]);
end



%converts electron volts to joules (plasma convention is to use T when they
%mean kT)
EV_TO_J = 1.60217653e-19;

%background permittivity in F/m = s^2 A^2 / m^3 kg
EPS_0 = 8.854187817e-12;

%Elementary charge
E_0 = 1.602176565e-19;


%number of species
nSp=length(mass);

%set spatial spacing and grid
dx = Lx/Nx;
x = 0.5*dx:dx:(Lx-0.5*dx);

%set velocity grid and velocity spacing
v = linspace(-Lv,Lv,Nv);
dv = v(2) - v(1);

%construct vector of weights for integration
wtN = dv*ones(Nv,1);
wtN(1) = 0.5*wtN(1);
wtN(end) = 0.5*wtN(end);

%permute this vector in several ways for use in computing density,
%temperature, and bulk velocity with vectorization in later steps
W=permute(repmat(wtN*wtN.',[1,1,nSp,Nx]),[3,4,1,2]);
W2=permute(repmat(wtN*(wtN.*v.').',[1,1,nSp,Nx]),[3,4,2,1]);
W3=permute(repmat((wtN.*v.')*wtN.',[1,1,nSp,Nx]),[3,4,2,1]);

%initialize array of forces from poisson
F = ones(nSp,Nx);

%initialize phi and ne0 from poisson
phi = zeros(Nx,1);
ne0=0;

%set initial distribution
f=f0;

%initialize density, temperature, and velocity arrays
density = zeros(nSp,Nx);

temperature = zeros(nSp,Nx);

if nSp==1
    u = zeros(Nx,2);
else
    u = zeros(nSp,Nx,2);
end

%construct matrices of velocities for comparisons in Maxwellian codes
V1=squeeze(repmat(v,[1,1,Nv]));
V2=V1.';


%repeated data for vectorized use in later calculations

charger=squeeze(repmat(charge,[1,1,Nx]));   %same dimension as density
massr=squeeze(repmat(mass,[1,1,Nx]));      %same dimension as density
V1adv=permute(repmat(V1,[1,1,Nx]),[3,1,2]); %same dimension as one row of f


if nSp==1
    
    %same dimension as one row of f
    V1r=permute(repmat(V1,[1,1,Nx]),[3,1,2]);
    V2r=permute(repmat(V2,[1,1,Nx]),[3,1,2]);
    
else
    
    %same dimension as f
    V1r=permute(repmat(V1,[1,1,nSp,Nx]),[3,4,1,2]);
    V2r=permute(repmat(V2,[1,1,nSp,Nx]),[3,4,1,2]);
    
end

%initialize maxwellian structures for collision terms
maxw = struct('V1',[],'V2',[],'mass',[],'density',[],'temperature',[],'velocity',[]);




maxw.V1=permute(repmat(V1,[1,1,Nx]),[3,1,2]);
maxw.V2=permute(repmat(squeeze(repmat(v,[1,1,Nv])).',[1,1,Nx]),[3,1,2]);

if nSp==1
    maxw.mass=squeeze(repmat(mass,[1,1,Nx,Nv,Nv]));
end





%initialize timing variables
if timing==1
    avgdenstime=0;
    avgbulkvtime=0;
    avgtemptime=0;
    avgadvect=0;
    avgpoisson=0;
    avgcollide=0;
end


%%%%%%%%%%%
%MAIN LOOP%
%%%%%%%%%%%

%iterate as many time steps as Ndt
for t = 1:Ndt
    
    %start a timer for each step
    if timing==1
        start0 = toc;
    end
    
    %compute the density with the integration weights in W
    density=sum(sum(W.*f,4),3);
    
    %record time (this comment will be omitted for future timing steps
    if timing==1
        avgdenstime=avgdenstime+(toc-start0)/(Ndt+1);
        start0=toc;
    end
    
    %compute the bulk velocity with integration weights
    if nSp==1
        u(:,1)=sum(sum(W2.*f,4),3);
        u(:,2)=sum(sum(W3.*f,4),3);
        u=u./squeeze(repmat(density,[1,1,2]));
    else
        u(:,:,1)=sum(sum(W2.*f,4),3);
        u(:,:,2)=sum(sum(W3.*f,4),3);
        u=u./repmat(density,[1,1,2]);
    end
    
    if timing==1
        avgbulkvtime=avgbulkvtime+(toc-start0)/(Ndt+1);
        start0=toc;
    end
    
    %compute the temperature
    if nSp==1
        
        temperature=sum(sum(squeeze(W).*(...
            (V1r-squeeze(repmat(u(:,1),[1,1,Nv,Nv]))).^2+...
            (V2r-squeeze(repmat(u(:,2),[1,1,Nv,Nv]))).^2).*...
            squeeze(f),3),2).*massr./(2*density.');
        
    else
        
        temperature=sum(sum(W.*(...
            (V1r-squeeze(repmat(u(:,:,1),[1,1,1,Nv,Nv]))).^2+...
            (V2r-squeeze(repmat(u(:,:,2),[1,1,1,Nv,Nv]))).^2).*...
            f,4),3).*massr./(2*density);
    end
    
    temperature(density==0)=0;
    
    if timing==1
        avgtemptime=avgtemptime+(toc-start0)/(Ndt+1);
    end
    
    %if we are plotting, plot every Njumps
    if plotting==1
        if(mod(t-1,Njumps) == 0 || t==1)
            figure(1)
            hold off
            
            if nSp==1
                plot(x,density,'Color',colors(1,:))
            else
                
                for i=1:nSp
                    plot(x,density(i,:),'Color',colors(i,:))
                    hold on
                end
            end
            xlabel('x')
            ylabel('n [m^{-3}]')
            if isnan(axissize)==0
                axis([0 Lx axissize(1,:)]);
            end
            legend(legends{1:nSp})
            title(sprintf('Timestep %d', t-1))
            
            
            figure(2)
            hold off
            if nSp==1
                plot(x,temperature/EV_TO_J,'Color',colors(1,:))
            else
                
                for i=1:nSp
                    plot(x,temperature(i,:)/EV_TO_J,'Color',colors(i,:))
                    hold on
                end
            end
            xlabel('x')
            ylabel('T [eV]')
            legend(legends{1:nSp})
            if isnan(axissize)==0
                axis([0 Lx axissize(2,:)]);
            end
            title(sprintf('Timestep %d', t-1))
            
            figure(3)
            hold off
            
            if nSp==1
                plot(x,u(:,1),'Color',colors(1,:))
            else
                
                for i=1:nSp
                    plot(x,u(i,:,1),'Color',colors(i,:))
                    hold on
                end
            end
            xlabel('x')
            ylabel('u [m/s]')
            if isnan(axissize)==0
                axis([0 Lx axissize(3,:)]);
            end
            legend(legends{1:nSp})
            title(sprintf('Timestep %d', t-1))
        end
    end
    
    %compute the transport step
    
    if timing==1
        start0=toc;
    end
    
    %compute the field from the poisson equation (see code below)
    nonlinPoisson;
    
    %compute the forces on each particle from this
    if nSp==1
        F=E_0*charge*mass*F;
    else
        F=(E_0*charger./massr).*squeeze(permute(repmat(F,[1,1,nSp]),[3,1,2]));
    end
    
    if timing==1
        avgpoisson=avgpoisson+(toc-start0)/Ndt;
        start0=toc;
    end
    
    %compute the advection step using F and the velocities (see code below)
    AdvectOne;
    %this is first order, second order code can be derived and implemented
    
    
    %if at any time our code encounters NaN, abort
    
    if( sum(isnan(fadv))  ~= 0)
        disp('NaN detected');
        break;
    end
    
    %compute the collision step
    BGKCollide;
    fcoll=0;
    
    if timing==1
        avgcollide=avgcollide+(toc-start0)/Ndt;
    end
    
    %advance the distribution using both the collision and advection terms
    f=f+fadv+fcoll;    
end

%compute the final density, temperature, and velocity
density=sum(sum(W.*f,4),3);

if timing==1
    avgdenstime=avgdenstime+(toc-start0)/(Ndt+1);
    start0=toc;
end

if nSp==1
    u(:,1)=sum(sum(W2.*f,4),3);
    u(:,2)=sum(sum(W3.*f,4),3);
    u=u./squeeze(repmat(density,[1,1,2]));
else
    u(:,:,1)=sum(sum(W2.*f,4),3);
    u(:,:,2)=sum(sum(W3.*f,4),3);
    u=u./repmat(density,[1,1,2]);
end

if timing==1
    avgbulkvtime=avgbulkvtime+(toc-start0)/(Ndt+1);
    start0=toc;
end

if nSp==1
    
    temperature=sum(sum(squeeze(W).*(...
        (V1r-squeeze(repmat(u(:,1),[1,1,Nv,Nv]))).^2+...
        (V2r-squeeze(repmat(u(:,2),[1,1,Nv,Nv]))).^2).*...
        squeeze(f),3),2).*massr./(2*density.');
    
else
    
    temperature=sum(sum(W.*(...
        (V1r-squeeze(repmat(u(:,:,1),[1,1,1,Nv,Nv]))).^2+...
        (V2r-squeeze(repmat(u(:,:,2),[1,1,1,Nv,Nv]))).^2).*...
        f,4),3).*massr./(2*density);
end

temperature(density==0)=0;

if timing==1
    avgtemptime=avgtemptime+(toc-start0)/(Ndt+1);
end


%print the total elapsed time
fprintf('Total time: %f secs\n', toc - start);

%if we are timing, print the average time spent on each step
if timing==1
    fprintf('Average density time: %f sec\n', avgdenstime)
    fprintf('Average bulk velocity time: %f sec\n', avgbulkvtime)
    fprintf('Average temperature time: %f sec\n', avgtemptime)
    fprintf('Average advection time: %f sec\n', avgadvect)
    fprintf('Average Poisson time: %f sec\n', avgpoisson)
    fprintf('Average collision time: %f sec\n', avgcollide)
end

%if we are plotting, plot the final results
if plotting==1
    
    figure(1)
    hold off
    
    if nSp==1
        plot(x,density,'Color',colors(1,:))
    else
        
        for i=1:nSp
            plot(x,density(i,:),'Color',colors(i,:))
            hold on
        end
    end
    xlabel('x')
    ylabel('n [m^{-3}]')
    if isnan(axissize)==0
        axis([0 Lx axissize(1,:)]);
    end
    legend(legends{1:nSp})
    title(sprintf('Timestep %d', t))
    
    
    figure(2)
    hold off
    if nSp==1
        plot(x,temperature/EV_TO_J,'Color',colors(1,:))
    else
        
        for i=1:nSp
            plot(x,temperature(i,:)/EV_TO_J,'Color',colors(i,:))
            hold on
        end
    end
    xlabel('x')
    ylabel('T [eV]')
    legend(legends{1:nSp})
    if isnan(axissize)==0
        axis([0 Lx axissize(2,:)]);
    end
    title(sprintf('Timestep %d', t))
    
    figure(3)
    hold off
    
    if nSp==1
        plot(x,u(:,1),'Color',colors(1,:))
    else
        
        for i=1:nSp
            plot(x,u(i,:,1),'Color',colors(i,:))
            hold on
        end
    end
    xlabel('x')
    ylabel('u [m/s]')
    if isnan(axissize)==0
        axis([0 Lx axissize(3,:)]);
    end
    legend(legends{1:nSp})
    title(sprintf('Timestep %d', t))
end






%%%%%%%%%%%%%%%%%%%%%%
%ADVECTION SUBROUTINE%
%%%%%%%%%%%%%%%%%%%%%%

%first order in space upwind advection scheme
    function AdvectOne
        
        %initialize advection array
        fadv=zeros(size(f));
        
        %update advection using upwinding and downwinding
        for i=1:nSp
            fadv(i,:,1:Nv/2,:)=-dt.*V1adv(:,1:Nv/2,:)/dx.*...
                squeeze(f(i,[2:Nx 1],1:Nv/2,:)-f(i,:,1:Nv/2,:));
            fadv(i,:,Nv/2+1:Nv,:)=-dt.*V1adv(:,Nv/2+1:Nv,:)/dx.*...
                squeeze(f(i,:,Nv/2+1:Nv,:)-f(i,[Nx 1:Nx-1],Nv/2+1:Nv,:));
        end
        
        %the same result is found with the below code, but it is slower for
        %small numbers of particle species:
        %    fadv(:,:,1:Nv/2,:)=-dt.*V1r(:,:,1:Nv/2,:)/dx.*...
        %        (f(:,[2:Nx 1],1:Nv/2,:)-f(:,:,1:Nv/2,:));
        %    fadv(:,:,Nv/2+1:Nv,:)=-dt.*V1r(:,:,Nv/2+1:Nv,:)/dx.*...
        %        (f(:,:,Nv/2+1:Nv,:)-f(:,[Nx 1:Nx-1],Nv/2+1:Nv,:));
        
        if timing==1
            avgadvect=avgadvect+(toc-start0)/Ndt;
            start0=toc;
        end
        
        %using the forces computed in the poisson subroutine, advance the
        %advection term using upstreaming and downstreaming
        Fpos=F>0;
        Fneg=F<0;
        Frep=repmat(F,[1,1,Nv,Nv]);
        
        for i=1:nSp
            fadv(i,Fpos(i,:),2:Nv,:)=fadv(i,Fpos(i,:),2:Nv,:)+...
                dt/dv*Frep(i,Fpos(i,:),2:Nv,:).*(f(i,Fpos(i,:),2:Nv,:)-f(i,Fpos(i,:),1:Nv-1,:));
            fadv(i,Fpos(i,:),1,:)=(1+dt/dv*Frep(i,Fpos(i,:),1,:)).*...
                fadv(i,Fpos(i,:),1,:);
            fadv(i,Fneg(i,:),1:Nv-1,:)=fadv(i,Fneg(i,:),1:Nv-1,:)+...
                dt/dv*Frep(i,Fneg(i,:),1:Nv-1,:).*(f(i,Fneg(i,:),2:Nv,:)-f(i,Fneg(i,:),1:Nv-1,:));
            fadv(i,Fneg(i,:),Nv,:)=(1-dt/dv*Frep(i,Fneg(i,:),Nv,:)).*...
                fadv(i,Fneg(i,:),Nv,:);
        end
    end

%%%%%%%%%%%%%%%%%%%%
%POISSON SUBROUTINE%
%%%%%%%%%%%%%%%%%%%%
%
% Solves nonlinear 1D Poisson equation
%
% - phi_xx + g(phi) = s
%
% with periodic BCs
% where g is the nonlinear function
%
% n_{e0} exp(E_0 phi / T)
%
% and s is a source from the ion terms
%
% The discretization is
% A_ij phi_j + g(phi_i) = s_i
% where A is the usual three point stencil for Laplacian

 function nonlinPoisson
        
        %set source
        if nSp==1
            source=charge*density;
        else
            source=sum(charger.*density,1);
        end
        source=source.';
 
        %find ne0 if this is the first time through
        if t==1
            Ci = sum(source)*dx;
            ne0 = (Ci/Lx)*ones(Nx,1);
        end
        
        %construct the second derivative discretization matrix
        A = (-diag(ones(Nx-1,1),-1) - diag(ones(Nx-1,1),1) + 2*diag(ones(Nx,1))) /dx/dx ;
        
        %set periodic boundary conditions
        A(1,end) = A(1,2);
        A(end,1) = A(end,end-1);
        
        %set electron temperature (presently arbitrary)
        T_e = 200*EV_TO_J;
        
        %Add the linear term
        A = A + diag(4*pi*E_0*ne0/T_e/EPS_0);
        
        %solve the linearized problem
        RHS = (4*pi/EPS_0)*(source - ne0);       
        phi = A\RHS;
        
        %the force is the x derivative of phi, use a second order
        %discretization to find this
        F=-(phi([2:Nx 1])-phi([Nx 1:Nx-1]))/(2*dx);
    end



%%%%%%%%%%%%%%%%%%%%%%
%COLLISION SUBROUTINE%
%%%%%%%%%%%%%%%%%%%%%%

    function BGKCollide
        
        %initialize collision array
        fcoll=zeros(size(f));
        fcoll2=zeros(size(f));
        
        %if the number of species types is one, we only have self
        %collisions
        if nSp==1
            
            %compute bmax^-2
            bmax2_inv= (charge^2*E_0*E_0*density.')./(EPS_0*temperature);
            
            %compute bmax^2
            bmax2=1.0./bmax2_inv;
            
            %compute bmin
            bmin=charge^2*E_0*E_0./temperature;
            
            %compute Lambda
            Lambda=bmax2./(bmin.*bmin);
            
            %compute relaxation parameter tau
            tau=3.0*sqrt(2*mass)*(2*pi*temperature).^(1.5)*EPS_0*EPS_0./(...
                (charge*charge*E_0*E_0)^2*0.5*log(1+Lambda.*Lambda));
            
            %compute lambda_ii
            lambda_ii=(density.')./tau;
            
            %update new maxwellian parameters
            maxw.density=squeeze(repmat(squeeze(density),[1,1,Nv,Nv]));
            maxw.velocity=squeeze(repmat(u,[1,1,1,Nv,Nv]));
            maxw.temperature=squeeze(repmat(squeeze(temperature),[1,1,Nv,Nv]));
            
            %compute self-collisions
            fcoll(i,:,:,:)=squeeze(fcoll2(i,:,:,:))+dt*squeeze(repmat(lambda_ii,[1,1,Nv,Nv]))...
                .*(maxwellian([],maxw,Nx,Nv,1)-squeeze(f(i,:,:,:)));
            
            %This version is more reader-friendly, but slower
            %
            %field1 = 'vrange'; value1 = v;
            %field2 = 'density'; value2 = squeeze(density);
            %field3 = 'temperature'; value3 = squeeze(temperature);
            %field4 = 'velocity'; value4 = u;
            %field5 = 'mass'; value5 = mass;
            %
            %s=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
            %
            %fcoll(i,:,:,:)=squeeze(fcoll(i,:,:,:))+dt*squeeze(repmat(lambda_ii,[1,1,Nv,Nv]))...
            %    .*(maxwellian(s)-squeeze(f(i,:,:,:)));
            
            
            
        else
            
            %compute bmax^-2
            bmax2_inv = sum(charger.^2*E_0*E_0.*density./(EPS_0*temperature));
            
            %compute bmax^2
            bmax2=1.0./bmax2_inv;
            
            %loop over all pairs of species types (order does not matter)
            for i=1:nSp
                for j=i:nSp
                    if i==j
                        %self-collisions
                        
                        %compute bmin
                        bmin=charge(i)^2*E_0*E_0./temperature(i,:);
                        
                        %compute Lambda
                        Lambda=bmax2./(bmin.*bmin);
                        
                        %compute relaxation parameter tau
                        tau=3.0*sqrt(2*mass(i))*(2*pi*temperature(i,:)).^(1.5)*EPS_0*EPS_0./(...
                            (charge(i)*charge(i)*E_0*E_0)^2*0.5*log(1+Lambda.*Lambda));
                        
                        %compute lambda_ii
                        lambda_ii=density(i,:)./tau;
              
                        %update new maxwellian parameters
                        maxw.density=squeeze(repmat(squeeze(density(i,:)),[1,1,Nv,Nv]));
                        maxw.velocity=squeeze(repmat(squeeze(u(i,:,:)),[1,1,1,Nv,Nv]));
                        maxw.temperature=squeeze(repmat(squeeze(temperature(i,:)),[1,1,Nv,Nv]));
                        maxw.mass=squeeze(repmat(mass(i),[1,1,Nx,Nv,Nv]));
                        
                        %compute self-collisions
                        fcoll(i,:,:,:)=squeeze(fcoll(i,:,:,:))+dt*squeeze(repmat(lambda_ii,[1,1,Nv,Nv]))...
                            .*(maxwellian([],maxw,Nx,Nv,1)-squeeze(f(i,:,:,:)));
                
                        %This version is more reader-friendly, but slower
                        %
                        %field1 = 'vrange'; value1 = linspace(-Lv,Lv,Nv);
                        %field2 = 'density'; value2 = squeeze(density(i,:));
                        %field3 = 'temperature'; value3 = squeeze(temperature(i,:));
                        %field4 = 'velocity'; value4 = squeeze(u(i,:,:));
                        %field5 = 'mass'; value5 = mass(i);
                        %
                        %si=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
                        %
                        %fcoll(i,:,:,:)=squeeze(fcoll(i,:,:,:))+dt*squeeze(repmat(lambda_ii,[1,1,Nv,Nv]))...
                        %    .*(maxwellian(si,[],Nx,Nv,1)-squeeze(f(i,:,:,:)));
         
                    else
                        %interparticle collisions
                        
                        %compute mass
                        mu=mass(i)*mass(j)/(mass(i)+mass(j));
                        
                        %compute mixT
                        mixT=(density(i,:).*temperature(i,:)+density(j,:).*temperature(j,:))...
                            ./(density(i,:)+density(j,:));
                        
                        %compute bmin
                        bmin=charge(i)*charge(j)*E_0*E_0./mixT;
                        
                        %compute Lambda
                        Lambda=bmax2./(bmin.*bmin);
                        
                        %compute relaxation parameter tau
                        tau=6.0*sqrt(2)*(pi*mixT./mu).^(1.5)*EPS_0*EPS_0.*mu*mass(i)./(...
                            (charge(i)*charge(j)*E_0*E_0)^2.*0.5.*log(1+Lambda.*Lambda));
                        
                        %compute lambda_ij
                        lambda_ij=density(j,:)./tau;
                        
                        %compute lambda_ji
                        lambda_ji=lambda_ij*mass(i)/mass(j).*density(i,:)./density(j,:);
                        
                        %repeate density for mixU_max calculation
                        densityr=repmat(density,[1,1,2]);
                        
                        %compute mixU_max
                        mixU_max=(mass(i)*densityr(i,:,:).*repmat(lambda_ij,[1,1,2])...
                            .*u(i,:,:)+mass(j)*densityr(j,:,:).*repmat(lambda_ji,[1,1,2]).*u(j,:,:))...
                            ./(mass(i)*densityr(i,:,:).*repmat(lambda_ij,[1,1,2])...
                            +mass(j)*densityr(j,:,:).*repmat(lambda_ji,[1,1,2]));
                        
                        %compute mixT_max
                        mixT_max = (lambda_ij.*density(i,:).*temperature(i,:)...
                            +lambda_ji.*density(j,:).*temperature(j,:)...
                            -(mass(i)*density(i,:).*lambda_ij.*sum(mixU_max.^2-u(i,:,:).^2,3)...
                            +mass(j)*density(j,:).*lambda_ji.*sum(mixU_max.^2-u(j,:,:).^2,3))/2.0)...
                            ./(density(i,:).*lambda_ij+density(j,:).*lambda_ji);
                        
                        maxw.density=squeeze(repmat(squeeze(density(i,:)),[1,1,Nv,Nv]));
                        maxw.velocity=squeeze(repmat(squeeze(mixU_max),[1,1,1,Nv,Nv]));
                        maxw.temperature=squeeze(repmat(squeeze(mixT_max),[1,1,Nv,Nv]));
                        maxw.mass=squeeze(repmat(mass(i),[1,1,Nx,Nv,Nv]));
                        
                        
                        fcoll(i,:,:,:)=squeeze(fcoll(i,:,:,:))+dt*squeeze(repmat(lambda_ij,[1,1,Nv,Nv]))...
                            .*(maxwellian([],maxw,Nx,Nv,1)-squeeze(f(i,:,:,:)));
                        
                        maxw.density=squeeze(repmat(squeeze(density(j,:)),[1,1,Nv,Nv]));
                        maxw.mass=squeeze(repmat(mass(j),[1,1,Nx,Nv,Nv]));
                        
                        fcoll(j,:,:,:)=squeeze(fcoll(j,:,:,:))+dt*squeeze(repmat(lambda_ji,[1,1,Nv,Nv]))...
                            .*(maxwellian([],maxw,Nx,Nv,1)-squeeze(f(j,:,:,:)));
                        
                        %This version is more reader-friendly, but slower
                        %
                        %field1 = 'vrange'; value1 = v;
                        %field2 = 'density'; value2 = squeeze(density(i,:));
                        %field3 = 'temperature'; value3 = squeeze(mixT_max);
                        %field4 = 'velocity'; value4 = squeeze(mixU_max);
                        %field5 = 'mass'; value5 = mass(i);
                        %
                        %sij=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
                        %
                        %field1 = 'vrange'; value1 = v;
                        %field2 = 'density'; value2 = squeeze(density(j,:));
                        %field3 = 'temperature'; value3 = squeeze(mixT_max);
                        %field4 = 'velocity'; value4 = squeeze(mixU_max);
                        %field5 = 'mass'; value5 = mass(j); 
                        %
                        %sji=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5); 
                        % 
                        %fcoll(i,:,:,:)=squeeze(fcoll(i,:,:,:))+dt*...
                        %    squeeze(repmat(lambda_ij,[1,1,Nv,Nv]))...
                        %    .*(maxwellian(sij,[],Nx,Nv,1)-squeeze(f(i,:,:,:)));      
                        %fcoll(j,:,:,:)=squeeze(fcoll(j,:,:,:))+dt*...
                        %    squeeze(repmat(lambda_ji,[1,1,Nv,Nv]))...
                        %    .*(maxwellian(sji,[],Nx,Nv,1)-squeeze(f(j,:,:,:)));
                        
                    end
                end
            end
        end
    end
end