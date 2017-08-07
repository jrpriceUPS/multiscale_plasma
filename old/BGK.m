function BGK(~)

%Matlab version of 1D-2V BGK code
%Uses SI Units 
%Assumes two species
%avg atom, second order to be added later
    
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
    m2 = 2.0*m1;                       %mass of species 2 [kg]

    Z1 = 1.0;                         %charge level of sp 1
    Z2 = 2.0;                         %charge level of sp 2


%For now, let's assume an interface problem

%left side ICs

    n1_l = 5.0e+28;                      %number density [m^-3]
    n2_l = 0.0;
    u1_l = [0.0,0.0];                   %bulk velocity [m/s]
    u2_l = [0.0,0.0];
    T1_l = 100*EV_TO_J;                 %kT [J]
    T2_l = 0;

%right side ICs

    n1_r = 0.0;                         %number density [m^-3]
    n2_r = 5.0e+28;
    u1_r = [0.0,0.0];                   %bulk velocity [m/s]
    u2_r = [0.0,0.0];
    T1_r = 0.0;                         %kT [J]
    T2_r = 200*EV_TO_J;

    therm_speed_max = sqrt(max([T1_l/m1 T2_l/m2 T1_r/m1 T2_r/m2]));

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
    dens2 = zeros(Nx,1);
    u1 = zeros(Nx,2);
    u2 = zeros(Nx,2);
    Temperature1 = zeros(Nx,1);
    Temperature2 = zeros(Nx,1);
    
%set up distribution functions

    f1     = zeros(Nx,Nv,Nv);
    f2     = zeros(Nx,Nv,Nv);
    f1_adv = zeros(Nx,Nv,Nv);
    f2_adv = zeros(Nx,Nv,Nv);

    %used if doing second order
    f1_old = zeros(Nx,Nv,Nv);
    f2_old = zeros(Nx,Nv,Nv);


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
        
        if(T2_l ~= 0)
            for j = 1:Nv
                for k = 1:Nv            
                    f2(i,j,k) = n2_l*(m2/(2.0*pi*T2_l)) * ...
                        exp(-0.5*m2*norm([v(j),v(k)] - u2_l,2)^2 / T2_l);
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
        
        if(T2_r ~= 0)
            for j = 1:Nv
                for k = 1:Nv            
                    f2(i,j,k) = n2_r*(m2/(2.0*pi*T2_r)) * ...
                        exp(-0.5*m2*norm([v(j),v(k)] - u2_r,2)^2 / T2_r);
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
            plot(x,dens1,'-b',x,dens2,'-g')
            figure(2)
            plot(x,Temperature1/EV_TO_J,'-b',x,Temperature2/EV_TO_J,'-g')
            pause(0.5);
        end

        %transport step

        %-------------------------%
        %first order in space/time%
        %-------------------------%
        
        nonlinPoisson;

        f1_adv = AdvectOne(f1,(E_0*Z1/m1)*F);
        f2_adv = AdvectOne(f2,(E_0*Z2/m2)*F);
        
        if( sum(isnan(f1_adv)) + sum(isnan(f2_adv)) ~= 0)
            disp('NaN detected');
            break;
        end

        %collision step
            
        BGKCollide;
                
        %update distribution function        
                        
        %NOTE - if collision rate is not assumed to be time dependent, this can be implicitly updated
        
        for i = 1:Nx
            f1(i,:,:) = f1(i,:,:) + f1_adv(i,:,:); + reshape(Q11 + Q12,[1, Nv, Nv]);
            f2(i,:,:) = f2(i,:,:) + f2_adv(i,:,:); + reshape(Q22 + Q21,[1, Nv, Nv]);   
        end
        
        %--------------------------%
        %second order in space/time%
        %--------------------------%

        %f1_adv = AdvectTwo(f1);
        %f2_adv = AdvectTwo(f2);

        %BGKCollide;

        %f1_old = f1;
        %f2_old = f2;

        %for i = 1:Nx
            %f1(i,:,:) = f1_old(i,:,:) + 0.5*(f1_adv(i,:,:) + reshape(Q11 + Q12,[1, Nv, Nv]));
            %f2(i,:,:) = f2_old(i,:,:) + 0.5*(f2_adv(i,:,:) + reshape(Q22 + Q21,[1, Nv, Nv]));   
        %end

        %GetDens;
        %GetBulkV;
        %GetTemperature;
        
        %f1_adv = AdvectTwo(f1);
        %f2_adv = AdvectTwo(f2);

        %BGKCollide;
        
        %for i = 1:Nx
            %f1(i,:,:) = f1_old(i,:,:) + f1_adv(i,:,:) + reshape(Q11 + Q12,[1, Nv, Nv]);
            %f2(i,:,:) = f2_old(i,:,:) + f2_adv(i,:,:) + reshape(Q22 + Q21,[1, Nv, Nv]);   
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
            dens2(i) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    dens1(i) = dens1(i) + wtN(j)*wtN(k)*f1(i,j,k);
                    dens2(i) = dens2(i) + wtN(j)*wtN(k)*f2(i,j,k);
                end
            end
        end
    end

    function GetBulkV
        for i = 1:Nx
            u1(i,1) = 0.0;
            u1(i,2) = 0.0;
            u2(i,1) = 0.0;
            u2(i,2) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    u1(i,1) = u1(i,1) + wtN(j)*wtN(k)*v(j)*f1(i,j,k);
                    u1(i,2) = u1(i,2) + wtN(j)*wtN(k)*v(k)*f1(i,j,k);
                    u2(i,1) = u2(i,1) + wtN(j)*wtN(k)*v(j)*f2(i,j,k);
                    u2(i,2) = u2(i,2) + wtN(j)*wtN(k)*v(k)*f2(i,j,k);
                end
            end
            u1(i,1) = u1(i,1)/dens1(i);
            u1(i,2) = u1(i,2)/dens1(i);
            u2(i,1) = u2(i,1)/dens2(i);
            u2(i,2) = u2(i,2)/dens2(i);
        end
    end

    function GetTemperature
        for i = 1:Nx
            Temperature1(i) = 0.0;
            Temperature2(i) = 0.0;
            for j = 1:Nv
                for k = 1:Nv
                    Temperature1(i) = Temperature1(i) + wtN(j)*wtN(k)*...
                    norm([(v(j) - u1(i,1)), (v(k) - u1(i,2))],2)^2 * f1(i,j,k);
                    Temperature2(i) = Temperature2(i) + wtN(j)*wtN(k)*...
                    norm([(v(j) - u2(i,1)), (v(k) - u2(i,2))],2)^2 * f2(i,j,k);
                end
            end

            if(dens1(i) ~= 0)
                Temperature1(i) = Temperature1(i)*m1/(2*dens1(i));
            else
                Temperature1(i) = 0.0;
            end
            
            if(dens2(i) ~= 0)
                Temperature2(i) = Temperature2(i)*m2/(2*dens2(i));
            else
                Temperature2(i) = 0.0;
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%
    %Advection subroutines%
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %first order in space upwind advection scheme
    function f_adv = AdvectOne(f,F)
       
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
        source = (4*pi/EPS_0)*(Z1*dens1 + Z2*dens2);
            
        A = (-diag(ones(Nx-1,1),-1) - diag(ones(Nx-1,1),1) + 2*diag(ones(Nx,1))) /dx/dx ;
        
        relErr = relTol + 1.0;
        absErr = absTol + 1.0;
        
        %set A's BCs for periodic
        A(1,end) = A(1,2);
        A(end,1) = A(end,end-1);
        
        
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
        % F = d_x phi
        for i = 2:Nx-1
            F(i) = (phi(i+1) - phi(i-1))/(2*dv);
        end
        F(1)  = (phi(2) - phi(Nx)  )/(2*dv);
        F(Nx) = (phi(1) - phi(Nx-2))/(2*dv);
        
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
            
            mu = m1*m2/(m1 + m2);             %reduced mass
            
            mixT = (dens1(i)*Temperature1(i) + dens2(i)* ...
                    Temperature2(i))/(dens1(i) + dens2(i));  %Mixture
                                                             %temperature
            bmax2_inv = Z1*Z1*E_0*E_0*dens1(i)/(EPS_0*Temperature1(i)) + Z2*Z2*E_0* ...
                E_0*dens2(i)/(EPS_0*Temperature2(i));                     %Debye
                                                                          %length,
                                                                          %kinda
            bmax2 = 1.0/bmax2_inv;
            
            %species 1 self collisions
            bmin = Z1*Z1*E_0*E_0/Temperature1(i);
            Lambda = bmax2/(bmin*bmin);
            tau11 = 12.0*pi*Temperature1(i)*EPS_0*EPS_0*sqrt(m1)/( ...
                (Z1*Z1*E_0*E_0)^2 * 0.5*log(1+Lambda*Lambda));
            lambda11 = dens1(i)/tau11;
            
            %species 1 collisions with species 2
            bmin = Z1*Z2*E_0*E_0/mixT;
            Lambda = bmax2/(bmin*bmin);
            tau12 = 6.0*pi*mixT*EPS_0*EPS_0*m1/( ...
                (Z1*Z2*E_0*E_0)^2 * sqrt(mu) * 0.5*log(1+Lambda*Lambda));
            lambda12 = dens2(i)/tau12;
            
            %species 2 collisions with species 1
            bmin = Z1*Z2*E_0*E_0/mixT;
            Lambda = bmax2/(bmin*bmin);
            tau21 = 6.0*pi*mixT*EPS_0*EPS_0*m2/( ...
                (Z1*Z2*E_0*E_0)^2 * sqrt(mu) * 0.5*log(1+Lambda*Lambda));
            lambda21 = dens1(i)/tau21;
            
            %species 2 self collisions
            bmin = Z2*Z2*E_0*E_0/Temperature2(i);
            Lambda = bmax2/(bmin*bmin);
            tau22 = 12.0*pi*Temperature2(i)*EPS_0*EPS_0*sqrt(m2)/( ...
                (Z2*Z2*E_0*E_0)^2 * 0.5*log(1+Lambda*Lambda));
            lambda22 = dens2(i)/tau22;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %calculate cross-species Maxwellian parameters
            
            mixU_Max = (m1*lambda12*u1(i,:) + m2*lambda21*u2(i,:))/(m1*lambda12 + m2*lambda21);
            
            
            mixT_Max = (lambda12*Temperature1(i) + lambda21*Temperature2(i) ...
                        - (m1*lambda12*(mixU_Max*mixU_Max' - u1(i,:)*u1(i,:)') ...
                        +  m2*lambda21*(mixU_Max*mixU_Max' - u2(i,:)*u2(i,:)'))/3.0 ) ...
                        /(lambda12 + lambda21);
            
            %update with collision terms
            Q11 = dt*lambda11*(GetMaxwell(m1,dens1(i),u1(i,:),Temperature1(i)) - squeeze(f1(i,:,:)));
            Q22 = dt*lambda22*(GetMaxwell(m2,dens2(i),u2(i,:),Temperature2(i)) - squeeze(f2(i,:,:)));
            Q12 = dt*lambda12*(GetMaxwell(m1,dens1(i),mixU_Max,mixT_Max)       - squeeze(f1(i,:,:)));
            Q21 = dt*lambda21*(GetMaxwell(m2,dens2(i),mixU_Max,mixT_Max)       - squeeze(f2(i,:,:)));
        end
    end
end
    