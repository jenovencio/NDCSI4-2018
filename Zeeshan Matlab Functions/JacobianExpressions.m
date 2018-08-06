function [f, ur, Jac] = Fric3Dcoup_Jac( x, w, Jac_r, nH, nt, Nt)

% this function evaluates the coupled tangential friction force and 
% Jacobian at one time step for all nodes in the displacement vector x.
% The function's inputs and outputs are:

% INPUTS: 



% sample x
x = [1; 1; 1; 2; 2; 2];

Nt = 100;               % no. of time steps - outerloop for Jacobian
nH = 2;                 % no. of harmonics
tH = 2*nH+1;            % total no. of harmonic coefficients
% create a vector of harmonic coefficient indices - for innerloop
nh = HarmIndexVector(nH);
% nh = 0:1:nH;

CM = 2*pi/Nt;           % constant multiplier
dpn = 3;                % degrees of freedom per node
N_node = length(x)/dpn; % total number of nodes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables from FricElem3D function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kn  = 50;
ktx = 50;           % tangential stiffness in x-direction
kty = 50;           % tangential stiffness in y-direction
N0 = 50;
mu = 0.2;

ux =  zeros(1, N_node);
uy =  zeros(1, N_node);
v  =  zeros(1, N_node);
uxr = zeros(1, N_node);
uyr = zeros(1, N_node);
Coul = zeros(1, N_node);
Tx =  zeros(1, N_node);
Ty =  zeros(1, N_node);
N  =  zeros(1, N_node);
F = zeros(dpn, N_node);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jac = zeros( dpn*N_node , dpn*(2*nH+1) );

for n = 1:N_node
    % dof's at the current time step
    ux(n) = x((dpn*n-2), 1);          % 1st dof in x-tangential dir
    uy(n) = x((dpn*n-1), 1);          % 2nd dof in y-tangential dir
    v(n) = x(dpn*n,1);
    % dof's at the previous time step
    uxr(n) = w((dpn*n-2), 1);
    uyr(n) = w((dpn*n-1), 1);
    
    % normal load variable. v = 0 makes it constant load case
    N(n) = max (N0 + kn * v(n), 0);

    if N(n) <= 0
        Tx(n) = 0;
        Ty(n) = 0;
        uxr(n) = ux(n);
        uyr(n) = uy(n);
        ur(:,n) = [uxr(n); uyr(n)];
        ID(1,n) = 0;         % separation 
        ID(2,n) = 0;         % separation
        
%         if JacOpt == 1
%             for j = 1:tH   
%                 % for x-force
%                 Jac(dpn*n-2, j) = 0;                  % dfxdX
%                 Jac(dpn*n-2, 1*tH+j) = 0;             % dfxdY
%                 Jac(dpn*n-2, 2*tH+j) = 0;             % dfxdX
%                 
%                 % for y-force
%                 Jac(dpn*n-1, j) = 0;                  % dfydX
%                 Jac(dpn*n-1, 1*tH+j) = 0;             % dfydY     
%                 Jac(dpn*n-1, 2*tH+j) = 0;             % dfydZ
%
%                 % for z-force
%                 Jac(dpn*n, j) = 0;                    % dfzdX       
%                 Jac(dpn*n, j) = 0;                    % dfzdY
%                 Jac(dpn*n, j) = 0;                    % dfzdZ
%             end
%         end
        
    else
        Coul(n) = mu*N(n); 
        Tx(n) = ktx *( ux(n) - uxr(n) );
        Ty(n) = kty *( uy(n) - uyr(n) );
        ID(1,n) = 1;         % stick phase
        ID(1,n) = 1;
        
        T(n) = sqrt( Tx(n)^2 + Ty(n)^2 );       % coupled tangential force
        d(n) = sqrt( (ux(n)-uxr(n))^2 + (uy(n)-uyr(n))^2 );        % direction cosine
        
        if JacOpt == 1
                        
            for j = 1:tH  
                % in stick condition
                if mod(j,2) == 0        % compute cosines (or real coefficients)
                    
                    % x-force
                    Jac(dpn*n-2, j) = kt * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) )...
                                      + Jac_r(dpn*n-2, j);          % dfxdX
                    Jac(dpn*n-2, 1*tH+j) = Jac_r(dpn*n-2, 1*tH+j);  % dfxdY
                    Jac(dpn*n-2, 2*tH+j) = Jac_r(dpn*n-2, 2*tH+j);  % dfxdZ
                    % y-force
                    Jac(dpn*n-1, j) = Jac_r(n,j);                   % dfydX
                    Jac(dpn*n-1, 1*tH+j) = kt * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) )...
                                           + dfydYr(n,j);           % dfydY
                    Jac(dpn*n-1, 2*tH+j) = Jac_r(dpn*n-1, 1*tH+j);  % dfydZ
                    
                    %z-force
                    Jac(dpn*n, j) = 0;                              % dfzdX
                    Jac(dpn*n, 1*tH+j) = 0;                         % dfzdY
                    Jac(dpn*n, 2*tH+j) = kn * cos(nh(j)*(nt)*CM);   % dfzdZ
                    

                else                    % compute sines (or imaginary coefficients)
                    % x-force
                    Jac(dpn*n-2, j) = kt * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) )...
                                      + Jac_r(dpn*n-2, j);          % dfxdX
                    Jac(dpn*n-2, 1*tH+j) = Jac_r(dpn*n-2,1*tH+j);   % dfxdY
                    Jac(dpn*n-2, 2*tH+j) = Jac_r(dpn*n-2, 2*tH+j);  % dfxdZ
                    
                    % y-force
                    Jac(dpn*n-1, j) = Jac_r(dpn*n-1, j);            % dfydX    
                    Jac(dpn*n-1, 1*tH+j) = kt * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) )...
                                           + Jac_r(dpn*n-1, 1*tH+j);% dfydZ
                    Jac(dpn*n-1, 2*tH+j) = Jac_r(dpn*n-1, 2*tH+j);  % dfydZ
                    
                    % z-force
                    Jac(dpn*n, j) = 0;                              % dfzdX
                    Jac(dpn*n, 1*tH+j) = 0;                         % dfzdY
                    Jac(dpn*n, 2*tH+j) = kn * -1*sin(nh(j)*(nt)*CM);% dfzdZ
                end
            end
        end
        
        
        
        if T(n) >= Coul(n)
            
            Tx(n) = Coul(n) * (ux(n)-uxr(n))/d(n);
            uxr(n) = ux(n) - Tx(n)/ktx;
            
            Ty(n) = Coul(n) * (uy(n)-uyr(n))/d(n);
            uyr(n) = uy(n) - Ty(n)/kty;
            
            ID(1,n) = 2;     % x slip  
            ID(2,n) = 2;     % x slip  
            
            
            if JacOpt == 1
                for j = 2:tH   
                    % in slip condition
                    if mod(j,2) == 0        % compute cosines (or real coefficients)

                        % x-force
                        Jac(dpn*n-2, j) = mu*N(n)* ( 1/d(n) * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) + 1/kt*Jac_r(dpn*n-2, j) ) ...
                            - ( (ux(n)-uxr(n))/d(n)^3 * (ux(n)-uxr(n)) * ( cos(nh(j)*(nt)*CM) - cos( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-2, j) )...
                            + (uy(n)-uyr(n))/kt * Jac_r(dpn*n-1,j) ) ); % dfxdX
                        
                         
                        Jac(dpn*n-2, 1*tH+j) = mu*N(n)* ( (uy(n)-uyr(n))/d(n)^3 * ((uy(n)-uyr(n))/kt * Jac_r(dpn*n-2, 1*tH+j) ...
                            - (ux(n)-uxr(n))* ( cos(nh(j)*(nt)*CM) - cos( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) )));  % dfxdY
                        
                        Jac(dpn*n-2, 2*tH+j) = mu*N(n)* ( (uy(n)-uyr(n))/d(n)^3 * ( (uy(n)-uyr(n))/kt * Jac_r(dpn*n-2, 2*tH+j) ...
                            - (ux(n)-uxr(n))/kt*Jac_r(dpn*n-1, 2*tH+j) )) ;  % dfxdZ
                        
                        % y-force
                        Jac(dpn*n-1, j) = mu*N(n)* ( (ux(n)-uxr(n))/d(n)^3 * ( (ux(n)-uxr(n))/kt * Jac_r(dpn*n-1, j) ...
                            - (uy(n)-uyr(n)) * ( cos(nh(j)*(nt)*CM) - cos( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-2, j)))) ;      %dfydX                % dfydX

                        Jac(dpn*n-1, 1*tH+j) = mu*N(n)* ( 1/d(n) * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) ) ...
                            - ( (uy(n)-uyr(n))/d(n)^3 * (uy(n)-uyr(n)) * ( cos(nh(j)*(nt)*CM) - cos( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) )...
                            + (ux(n)-uxr(n))/kt * Jac_r(dpn*n-2,1*tH+j) ) );           % dfydY
                  
                        Jac(dpn*n-1, 2*tH+j) = mu*N(n)* ( (ux(n)-uxr(n))/d(n)^3 * ( (ux(n)-uxr(n))/kt * Jac_r(dpn*n-1, 2*tH+j) ...
                            - (uy(n)-uyr(n))/kt*Jac_r(dpn*n-2, 2*tH+j) )) ;  % dfydZ
                        
                        
                        

                        %z-force
                        Jac(dpn*n, j) = 0;                              % dfzdX
                        Jac(dpn*n, 1*tH+j) = 0;                         % dfzdY
                        Jac(dpn*n, 2*tH+j) = kn * cos(nh(j)*nt*CM);     % dfzdZ


                    else                    % compute sines (or imaginary coefficients)
                        % x-force
                        Jac(dpn*n-2, j) = mu*N(n)* ( 1/d(n) * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) + 1/kt*Jac_r(dpn*n-2, j) ) ...
                            - ( (ux(n)-uxr(n))/d(nt)^3 * (ux(n)-uxr(n)) * ( -sin(nh(j)*(nt)*CM) + sin( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-2, j) )...
                            + (uy(n)-uyr(n))/kt * Jac_r(dpn*n-1,j) ) ); % dfxdX
                        
                        
                        Jac(dpn*n-2, 1*tH+j) = mu*N(n)* ( (uy(n)-uyr(n))/d(n)^3 * ( (uy(n)-uyr(n))/kt * Jac_r(dpn*n-2, 1*tH+j) ...
                            - (ux(n)-uxr(n))* ( -sin(nh(j)*(nt)*CM) + sin( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) )));   % dfxdY
                        
                        Jac(dpn*n-2, 2*tH+j) = mu*N(n)* ( (uy(n)-uyr(n))/d(n)^3 * ( (uy(n)-uyr(n))/kt * Jac_r(dpn*n-2, 2*tH+j) ...
                            - (ux(n)-uxr(n))/kt*Jac_r(dpn*n-1, 2*tH+j) )) ;  % dfxdZ

                        % y-force
                        Jac(dpn*n-1, j) = mu*N(n)* ( (ux(n)-uxr(n))/d(n)^3 * ( (ux(n)-uxr(n))/kt * Jac_r(dpn*n-1, j) ...
                            - (uy(n)-uyr(n)) * ( cos(nh(j)*(nt)*CM) - cos( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-2, j)))) ;      % dfydX 
                        
                           
                        Jac(dpn*n-1, 1*tH+j) = mu*N(n)* ( 1/d(n) * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) ) ...
                            - ( (uy(n)-uyr(n))/d(n)^3 * (uy(n)-uyr(n)) * ( -sin(nh(j)*(nt)*CM) + sin( nh(j)*(nt-1)*CM ) + 1/kt*Jac_r(dpn*n-1, 1*tH+j) )...
                            + (ux(n)-uxr(n))/kt * Jac_r(dpn*n-2,1*tH+j) ) );           % dfydY
                        
                        Jac(dpn*n-1, 2*tH+j) = mu*N(n)* ( (ux(n)-uxr(n))/d(n)^3 * ( (ux(n)-uxr(n))/kt * Jac_r(dpn*n-1, 2*tH+j) ...
                            - (uy(n)-uyr(n))/kt*Jac_r(dpn*n-2, 2*tH+j) )) ;  % dfydZ

                        
                        
                        % z-force
                        Jac(dpn*n, j) = 0;                              % dfzdX
                        Jac(dpn*n, 1*tH+j) = 0;                         % dfzdY
                        Jac(dpn*n, 2*tH+j) = kn * -1*sin(nh(j)*nt*CM);     % dfzdZ
                    end
                end
            end
            
        end
        

        f(:,n) = [Tx(n); Ty(n); N(n)-N0];       
        ur(:,n) = [uxr(n); uyr(n)];
    end
end
    
    
 