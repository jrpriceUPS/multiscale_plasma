function BGK_1spec(~)

%Matlab version of 1D-2V BGK code
%Uses SI Units 
%Assumes ONE species
%Avg Atom model, second order to be included later
    
%%%%%%%%%%%%%%%%%%
%useful constants%
%%%%%%%%%%%%%%%%%%
    
    %converts electron volts to joules (since plasma convention is to use T
    %then they mean kT)
    EV_TO_J = 1.60217653e-19;
        
    %background permittivity in F/m = s^2 A^2 / m^3 kg
    EPS_0 = 8.854187817e-12;
        
    %Elementary charge 
    E_0 = 1.602176565e-19;
        


%%%%%%%%%%%%%%%%%%%%
%Initial conditions%    
%%%%%%%%%%%%%%%%%%%%

%Species info
    m1 = 1.67377e-27;                  %mass of species 1 [kg]

    Z1 = 1.0;                          %charge level of sp 1


%For now, let's assume an interface problem

%left side ICs

    n1_l = 5.0e+28;                      %number density [m^-3]
    u1_l = [0.0,0.0];                   %bulk velocity [m/s]
    T1_l = 100*EV_TO_J;                 %kT [J]

%right side ICs

    n1_r = 0.0;                         %number density [m^-3]
    u1_r = [0.0,0.0];                   %bulk velocity [m/s]
    T1_r = 0.0;                         %kT [J]

    therm_speed_max = sqrt(max([T1_l/m1 T1_r/m1]));

%%%%%%%%%%%%    
%Parameters%
%%%%%%%%%%%%

%Time parameters
    dt = 1.0e-15;                       %timestep [seconds]
    Ndt = 1000;                        %Total timesteps in sim
    Njumps = 1;                       %Number of timesteps between
                                       %dumping data

%Physical space parameters
    Lx = 800e-9;                      %Width of computational box [m]
    Nx = 32;                          %Number of spatial cells
    dx = Lx/Nx;
    x = 0.5*dx:dx:(Lx-0.5*dx);

%Velocity space parameters
    Lv = 10.0*therm_speed_max;          %Max velocity [m/s]
    Nv = 60;                           %Number of velocity space
                                       %grid points in one
                                       %direction
    
    v = linspace(-Lv,Lv,Nv);           %assuming this is centered
                                       %at zero, can be generalized

    dv = v(2) - v(1);                  %velocity spacing
    wtN = dv*ones(1,Nv);               %integration weights (trapezoid)
    wtN(1) = 0.5*wtN(1);
    wtN(end) = 0.5*wtN(end);
    
    %thus the velocity vector is given by [v(i), v(j)]

%Poisson solver iteration parameters + setup
    relTol = 1e-8;
    absTol = 1e-8;
    phi = ones(Nx,1);
    F = ones(Nx,1);

%Other initializations
    dens1 = zeros(Nx,1);
    u1 = zeros(Nx,2);
    Temperature1 = zeros(Nx,1);
    
%set up distribution functions

    f1     = zeros(Nx,Nv,Nv);
    f1_adv = zeros(Nx,Nv,Nv);

    %used if doing second order
    f1_old = zeros(Nx,Nv,Nv);


    for i = 1:Nx/2
        %fill left side Maxwellians
        if(T1_l ~= 0)
            
            for j = 1:Nv
                for k = 1:Nv
                    f1(i,j,k) = n1_l*(m1/(2.0*pi*T1_l)) * ...
                        exp(-0.5*m1*norm([v(j),v(k)] - u1_l,2)^2 / T1_l);
                end
            end
        end
        
    end
    
    for i = (Nx/2 + 1):Nx
        %fill right side Maxwellians
        if(T1_r ~= 0)
            
            for j = 1:Nv
                for k = 1:Nv
                    f1(i,j,k) = n1_r*(m1/(2.0*pi*T1_r)) * ...
                        exp(-0.5*m1*norm([v(j),v(k)] - u1_r,2)^2 / T1_r);
                end
            end
        end
        
    end

    %initial conditions set

    
%%%%%%%%%%%
%MAIN LOOP%
%%%%%%%%%%%

    for t = 1:Ndt
        
        GetDens;
        GetBulkV;
        GetTemperature;

        t

        if(mod(t,Njumps) == 0)
            figure(1)
            plot(x,dens1,'-b')
            figure(2)
            plot(x,Temperature1/EV_TO_J,'-b')
            pause(0.5);
        end

        %transport step

        %-------------------------%
        %first order in space/time%
        %-------------------------%

        nonlinPoisson;

        f1_adv = AdvectOne(f1);
        
        if( sum(isnan(f1_adv)) ~= 0)
            disp('NaN detected');
            break;
        end

        %collision step
            
        BGKCollide;
                
        %update distribution function        
                        
        %NOTE - if collision rate is not assumed to be time dependent, this can be implicitly updated
        
        for i = 1:Nx
            f1(i,:,:) = f1(i,:,:) + f1_adv(i,:,:); + reshape(Q11,[1, Nv, Nv]);
        end
        
        %--------------------------%
        %second order in space/time%
        %--------------------------%

        %f1_adv = AdvectTwo(f1);

        %BGKCollide;

        %f1_old = f1;

        %for i = 1:Nx
            %f1(i,:,:) = f1_old(i,:,:) + 0.5*(f1_adv(i,:,:) + reshape(Q11 + Q12,[1, Nv, Nv]));
        %end

        %GetDens;
        %GetBulkV;
        %GetTemperature;
        
        %f1_adv = AdvectTwo(f1);

        %BGKCollide;
        
        %for i = 1:Nx
            %f1(i,:,:) = f1_old(i,:,:) + f1_adv(i,:,:) + reshape(Q11 + Q12,[1, Nv, Nv]);
        %end

    end
        
    
    GetDens;
    GetBulkV;
    GetTemperature;



%%%%%%%%%%%%%%%%%%
%Nested functions%
%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%
    %Moment subroutines%
    %%%%%%%%%%%%%%%%%%%%
    
    function GetDens
        for i = 1:Nx
            dens1(i) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    dens1(i) = dens1(i) + wtN(j)*wtN(k)*f1(i,j,k);
                end
            end
        end
    end

    function GetBulkV
        for i = 1:Nx
            u1(i,1) = 0.0;
            u1(i,2) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    u1(i,1) = u1(i,1) + wtN(j)*wtN(k)*v(j)*f1(i,j,k);
                    u1(i,2) = u1(i,2) + wtN(j)*wtN(k)*v(k)*f1(i,j,k);
                end
            end
            u1(i,1) = u1(i,1)/dens1(i);
            u1(i,2) = u1(i,2)/dens1(i);
        end
    end

    function GetTemperature
        for i = 1:Nx
            Temperature1(i) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    Temperature1(i) = Temperature1(i) + wtN(j)*wtN(k)*...
                    norm([(v(j) - u1(i,1)), (v(k) - u1(i,2))],2)^2 * f1(i,j,k);
                end
            end

            if(dens1(i) ~= 0)
                Temperature1(i) = Temperature1(i)*m1/(2*dens1(i));
            else
                Temperature1(i) = 0.0;
            end
            
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    %Advection subroutines%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %first order in space upwind advection scheme
    function f_adv = AdvectOne(f)
       
        for i = 2:Nx-1
            for k = 1:Nv
                for j = 1:Nv/2              %stuff moving left
                    %upwinding
                    f_adv(i,j,k) = -(dt*v(j)/dx)*(f(i+1,j,k) - f(i,j,k));
                end
                for j = (Nv/2+1):Nv         %stuff moving right
                    f_adv(i,j,k) = -(dt*v(j)/dx)*(f(i,j,k) - f(i-1,j,k));
                end                
            end
        end

        %periodic BC Handling - i=1
        for k = 1:Nv
            for j = 1:Nv/2
                f_adv(1,j,k) = -(dt*v(j)/dx)*(f(2,j,k) - f(1,j,k));
            end
            for j = (Nv/2+1):Nv
                f_adv(1,j,k) = -(dt*v(j)/dx)*(f(1,j,k) - f(Nx,j,k));
            end            
        end

        %periodic BC Handling - i=Nx
        for k = 1:Nv
            for j = 1:Nv/2
                f_adv(Nx,j,k) = -(dt*v(j)/dx)*(f(1,j,k) - f(Nx,j,k));
            end
            for j = (Nv/2+1):Nv
                f_adv(Nx,j,k) = -(dt*v(j)/dx)*(f(Nx,j,k) - f(Nx-1,j,k));
            end
        end
        
        %poisson terms
        for i= 1:Nx
            for j = 1:Nv
                for k = 1:Nv
                    if(F(i) > 0)
                        if(j ~= 1)
                            f_adv(i,j,k) = f_adv(i,j,k) + (dt*F(i)/dv)*(f(i,j,k) - f(i,j-1,k));
                        else
                            f_adv(i,j,k) = f_adv(i,j,k) + (dt*F(i)/dv)*(f(i,j,k));
                        end
                    else
                        if(j ~= Nv)
                            f_adv(i,j,k) = f_adv(i,j,k) + (dt*F(i)/dv)*(f(i,j+1,k) - f(i,j,k));
                        else
                            f_adv(i,j,k) = f_adv(i,j,k) + (dt*F(i)/dv)*(-f(i,j,k));
                        end
                    end
                end
            end
        end
    end


    %NOTE - poisson terms need to be added here...
    
    function f_adv = AdvectTwo(f)
        slope = zeros(1,3);

        for i = 3:Nx-2
            for k = 1:Nv
                for j = 1:Nv/2  %stuff moving left
                    slope(2) = minmod((f(i,j,k)   - f(i-1,j,k))/dx, ...
                                      (f(i+1,j,k) - f(i,j,k))/dx, ...
                                      (f(i+1,j,k) - f(i-1,j,k))/(2*dx));
                    slope(1) = minmod((f(i-1,j,k) - f(i-2,j,k))/dx, ...
                                      (f(i,j,k)   - f(i-1,j,k))/dx, ...
                                      (f(i,j,k)   - f(i-2,j,k))/(2*dx));
                    f_adv(i,j,k) = f(i,j,k) - (dt*v(j)/dx)*( ...
                        f(i,j,k) + 0.5*dx*slope(2) - f(i-1,j,k) - ...
                        0.5*dx*slope(1));
                end
                for j = 1:Nv/2  %stuff moving right
                    slope(3) = minmod((f(i+1,j,k) - f(i,j,k))/dx, ...
                                      (f(i+2,j,k) - f(i+1,j,k))/dx, ...
                                      (f(i+2,j,k) - f(i,j,k))/(2*dx));
                    slope(2) = minmod((f(i,j,k)   - f(i-1,j,k))/dx, ...
                                      (f(i+1,j,k) - f(i,j,k))/dx, ...
                                      (f(i+1,j,k) - f(i-1,j,k))/(2*dx));
                    f_adv(i,j,k) = f(i,j,k) - (dt*v(j)/dx)*( ...
                        f(i+1,j,k) + 0.5*dx*slope(3) - f(i,j,k) - ...
                        0.5*dx*slope(2));
                end
            end
        end

        %boundary cells
        %i=1
        for k = 1:Nv
            for j = 1:Nv/2  %stuff moving left
                slope(2) = minmod((f(1,j,k)   - f(Nx,j,k))/dx, ...
                                  (f(2,j,k)   - f(1,j,k))/dx, ...
                                  (f(2,j,k)   - f(Nx,j,k))/(2*dx));
                slope(1) = minmod((f(Nx,j,k)  - f(Nx-1,j,k))/dx, ...
                                  (f(1,j,k)   - f(Nx,j,k))/dx, ...
                                  (f(1,j,k)   - f(Nx-1,j,k))/(2*dx));
                f_adv(1,j,k) = f(1,j,k) - (dt*v(j)/dx)*( ...
                    f(1,j,k) + 0.5*dx*slope(2) - f(Nx,j,k) - ...
                    0.5*dx*slope(1));
            end
            for j = 1:Nv/2  %stuff moving right
                slope(3) = minmod((f(2,j,k) - f(1,j,k))/dx, ...
                                  (f(3,j,k) - f(2,j,k))/dx, ...
                                  (f(3,j,k) - f(1,j,k))/(2*dx));
                slope(2) = minmod((f(1,j,k) - f(Nx,j,k))/dx, ...
                                  (f(2,j,k) - f(1,j,k))/dx, ...
                                  (f(2,j,k) - f(Nx,j,k))/(2*dx));
                f_adv(1,j,k) = f(1,j,k) - (dt*v(j)/dx)*( ...
                    f(2,j,k) + 0.5*dx*slope(3) - f(1,j,k) - ...
                    0.5*dx*slope(2));
            end
        end
        %i=2
        for k = 1:Nv
            for j = 1:Nv/2  %stuff moving left
                slope(2) = minmod((f(2,j,k)   - f(1,j,k))/dx, ...
                                  (f(3,j,k)   - f(2,j,k))/dx, ...
                                  (f(3,j,k)   - f(1,j,k))/(2*dx));
                slope(1) = minmod((f(1,j,k)   - f(Nx,j,k))/dx, ...
                                  (f(2,j,k)   - f(1,j,k))/dx, ...
                                  (f(2,j,k)   - f(Nx,j,k))/(2*dx));
                f_adv(2,j,k) = f(2,j,k) - (dt*v(j)/dx)*( ...
                    f(2,j,k) + 0.5*dx*slope(2) - f(1,j,k) - ...
                    0.5*dx*slope(1));
            end
            for j = 1:Nv/2  %stuff moving right
                slope(3) = minmod((f(3,j,k) - f(2,j,k))/dx, ...
                                  (f(4,j,k) - f(3,j,k))/dx, ...
                                  (f(4,j,k) - f(2,j,k))/(2*dx));
                slope(2) = minmod((f(2,j,k) - f(1,j,k))/dx, ...
                                  (f(3,j,k) - f(2,j,k))/dx, ...
                                  (f(3,j,k) - f(1,j,k))/(2*dx));
                f_adv(2,j,k) = f(2,j,k) - (dt*v(j)/dx)*( ...
                    f(3,j,k) + 0.5*dx*slope(3) - f(2,j,k) - ...
                    0.5*dx*slope(2));
            end
        end

        %i=Nx
        for k = 1:Nv
            for j = 1:Nv/2  %stuff moving left
                slope(2) = minmod((f(Nx,j,k)   - f(Nx-1,j,k))/dx, ...
                                  (f(1,j,k)    - f(Nx,j,k))/dx, ...
                                  (f(1,j,k)    - f(Nx-1,j,k))/(2*dx));
                slope(1) = minmod((f(Nx-1,j,k) - f(Nx-2,j,k))/dx, ...
                                  (f(Nx,j,k)   - f(Nx-1,j,k))/dx, ...
                                  (f(Nx,j,k)   - f(Nx-2,j,k))/(2*dx));
                f_adv(Nx,j,k) = f(Nx,j,k) - (dt*v(j)/dx)*( ...
                    f(Nx,j,k) + 0.5*dx*slope(2) - f(Nx-1,j,k) - ...
                    0.5*dx*slope(1));
            end
            for j = 1:Nv/2  %stuff moving right
                slope(3) = minmod((f(1,j,k)  - f(Nx,j,k))/dx, ...
                                  (f(2,j,k)  - f(1,j,k))/dx, ...
                                  (f(2,j,k)  - f(Nx,j,k))/(2*dx));
                slope(2) = minmod((f(Nx,j,k) - f(Nx-1,j,k))/dx, ...
                                  (f(1,j,k)  - f(Nx,j,k))/dx, ...
                                  (f(1,j,k)  - f(Nx-1,j,k))/(2*dx));
                f_adv(Nx,j,k) = f(Nx,j,k) - (dt*v(j)/dx)*( ...
                    f(1,j,k) + 0.5*dx*slope(3) - f(Nx,j,k) - ...
                    0.5*dx*slope(2));
            end
        end
        %i=Nx-1
        for k = 1:Nv
            for j = 1:Nv/2  %stuff moving left
                slope(2) = minmod((f(Nx-1,j,k)   - f(Nx-2,j,k))/dx, ...
                                  (f(Nx,j,k)    - f(Nx-1,j,k))/dx, ...
                                  (f(Nx,j,k)    - f(Nx-2,j,k))/(2*dx));
                slope(1) = minmod((f(Nx-2,j,k) - f(Nx-3,j,k))/dx, ...
                                  (f(Nx-1,j,k)   - f(Nx-2,j,k))/dx, ...
                                  (f(Nx-1,j,k)   - f(Nx-3,j,k))/(2*dx));
                f_adv(Nx-1,j,k) = f(Nx-1,j,k) - (dt*v(j)/dx)*( ...
                    f(Nx-1,j,k) + 0.5*dx*slope(2) - f(Nx-2,j,k) - ...
                    0.5*dx*slope(1));
            end
            for j = 1:Nv/2  %stuff moving right
                slope(3) = minmod((f(Nx,j,k)  - f(Nx-1,j,k))/dx, ...
                                  (f(1,j,k)  - f(Nx,j,k))/dx, ...
                                  (f(1,j,k)  - f(Nx-1,j,k))/(2*dx));
                slope(2) = minmod((f(Nx-1,j,k) - f(Nx-2,j,k))/dx, ...
                                  (f(Nx,j,k)  - f(Nx-1,j,k))/dx, ...
                                  (f(Nx,j,k)  - f(Nx-2,j,k))/(2*dx));
                f_adv(Nx-1,j,k) = f(Nx-1,j,k) - (dt*v(j)/dx)*( ...
                    f(Nx,j,k) + 0.5*dx*slope(3) - f(Nx-1,j,k) - ...
                    0.5*dx*slope(2));
            end
        end

    end

    function minmod_ret = minmod(a,b,c)
        if( (a > 0) && (b > 0) && (c > 0))
            minmod_ret = min([a,b,c]);
        elseif ( (a < 0) && (b < 0) && (c < 0))
            minmod_ret = max([a,b,c]);
        else
            minmod_ret = 0;
        end
            
    end

    %%%%%%%%%%%%%%%%%%%%%
    %Poisson solve stuff%
    %%%%%%%%%%%%%%%%%%%%%

    % solves nonlinear 1D Poisson equation
    %
    % - phi_xx + g(phi) = s
    %
    % with periodic BCs
    % where g is the nonlinear function 
    %
    % n_e0 exp(E_0 phi / T)
    %
    % and s is a source from the ion terms.
    % 
    % The discretization is
    % A_ij phi_j + g(phi_i) = s_i
    % where A is the usual three point stencil for Laplacian

    function nonlinPoisson
            
        %set source
        source = (4*pi/EPS_0)*(Z1*dens1);
            
        A = (-diag(ones(Nx-1,1),-1) - diag(ones(Nx-1,1),1) + 2*diag(ones(Nx,1))) /dx/dx ;
        
        relErr = relTol + 1.0;
        absErr = absTol + 1.0;
        
        %set A's BCs for periodic
        A(1,end) = A(1,2);
        A(end,1) = A(end,end);
        
        
        while ( ( relErr > relTol ) || ( absErr > absTol ) )
            [g,gPrime] = electronSource(phi,source);
            B = A + diag(gPrime);
            rhs = source - g + gPrime.*phi;
            phiNext = B\rhs;
            dphi = phiNext - phi;
            relErr = norm(dphi)/norm(phi); 
            absErr = norm(dphi);
            phi = phiNext;
        end

        %we want to return the force for the acceleration
        % F = (Z1*E_0/m1) d_x phi
        for i = 2:Nx-1
            F(i) = (E_0*Z1/m1)*(phi(i+1) - phi(i-1))/(2*dv);
        end
        F(1)  = (E_0*Z1/m1)*(phi(2) - phi(Nx)  )/(2*dv);
        F(Nx) = (E_0*Z1/m1)*(phi(1) - phi(Nx-2))/(2*dv);
        
    end
        

    function [g , gPrime] = electronSource(phi, ne0)
            
        Te  = 200*EV_TO_J;            %electron temperature - this should be improved later...
        
        g =          ne0   .*exp(E_0*phi/Te);
        gPrime = E_0*ne0/Te.*exp(E_0*phi/Te);        
        
    end    



    %%%%%%%%%%%%%%%%%%%%%%%
    %Collision subroutines%
    %%%%%%%%%%%%%%%%%%%%%%%

    function M = GetMaxwell(m,n,u,T)
        M = zeros(Nv,Nv);
        for j = 1:Nv
            for k = 1:Nv
                M(i,j) = n*(m/(2*pi*T)) * exp(-0.5*m/T * norm([v(j),v(k)] - u,2)^2);   
            end
        end
    end

    function BGKCollide
        for i = 1:Nx            
            %placeholder stuff in here for now - MD informed collision rates will replace this
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      
            bmax2_inv = Z1*Z1*E_0*E_0*dens1(i)/(EPS_0*Temperature1(i));   %Debye
                                                                          %length,
                                                                          %kinda
            bmax2 = 1.0/bmax2_inv;
            
            %species 1 self collisions
            bmin = Z1*Z1*E_0*E_0/Temperature1(i);
            Lambda = bmax2/(bmin*bmin);
            tau11 = 12.0*pi*Temperature1(i)*EPS_0*EPS_0*sqrt(m1)/( ...
                (Z1*Z1*E_0*E_0)^2 * 0.5*log(1+Lambda*Lambda));
            lambda11 = dens1(i)/tau11;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %calculate cross-species Maxwellian parameters
                        
            %update with collision terms
            Q11 = dt*lambda11*(GetMaxwell(m1,dens1(i),u1(i,:),Temperature1(i)) - squeeze(f1(i,:,:)));

        end
    end
end
    