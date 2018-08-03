
% sample x
x = [1; 1; 1; 2; 2; 2];

Nt = 100;               % no. of time steps - outerloop for Jacobian
nH = 2;                 % no. of harmonics
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

Jaco_time = zeros( dpn*N_node , dpn*(2*nH+1) );

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
%             for j = 1:(2*nH+1)    
%                 dfxdX(n,j) = 0;     dfxdY(n,j) = 0;     dfxdZ(n,j) = 0;
%                 dfydX(n,j) = 0;     dfydY(n,j) = 0;     dfydZ(n,j) = 0;
%                 dfzdX(n,j) = 0;     dfzdY(n,j) = 0;     dfzdZ(n,j) = 0;
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
            
            dfxdX(n,1) = kt * ( cos(nh(1)*(nt)*CM) - cos(nh(1)*(nt-1)*CM) ) + dfxdXr(n,1); 
            Jac(n,1) = kt + Jac_r(n,1); 
%             dfxdY(n,1) = dfxdYr(n,1);
%             dfxdZ(n,(2*nH+1)+1) = dfxdZr(n,(2*nH+1)+1);             
%                          
%             dfydX(n,(2*nH+1)+1) = kt * ( cos(nh(1)*(nt)*CM) - cos(nh(1)*(nt-1)*CM) )...
%                              + dfxdXr(n,1); 
            for j = 2:(2*nH+1)   
                % in stick condition
                if mod(j,2) == 0
                     dfxdX(n,j) = kt * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) ) + dfxdXr(n,j); 
                     dfxdY(n,j) = dfxdYr(n,j);
                     dfxdZ(n,j) = dfxdZr(n,j);
                     
                     dfydX(n,j) = dfydXr(n,j);
                     dfydY(n,j) = kt * ( cos(nh(j)*(nt)*CM) - cos(nh(j)*(nt-1)*CM) ) + dfydYr(n,j);
                     dfydZ(n,j) = dfydZr(n,j);
                
                else
                     dfxdX(n,j) = kt * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) ) + dfxdXr(nt); 
                     dfxdY(n,j) = dfxdYr(n,j);
                     dfxdZ(n,j) = dfxdZr(n,j);
                     
                     dfydX(n,j) = dfydXr(n,j);
                     dfydY(n,j) = kt * ( -sin(nh(j)*(nt)*CM) + sin(nh(j)*(nt-1)*CM) ) + dfydYr(n,j);
                     dfydZ(n,j) = dfydZr(n,j);
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
        end
        

        F(:,n) = [Tx(n); Ty(n); N(n)-N0];       
        ur(:,n) = [uxr(n); uyr(n)];
    end
end
    
    
    
    
for j = 1

    % if separation
    dfxdX(j) = 0;
    
    dfxdY(j) = 0;


    % if stick
    

    dfxdY(j) = dfxdY(nt-1);



    % if slip
    dfxdX(j) = mu*N(nt)* ( 1/d(nt) * ( cos( nh*(nt)*CM) - cos( nh*(nt-1)*CM) + 1/kt* dfxdX(nt-1))...
        - Xr(nt)/d(nt)^3 * ( Xr(nt) * ...
        ( cos( nh*(nt)*CM) - cos( nh*(nt-1)*CM) + 1/kt* dfxdX(nt-1)) + Yr(nt)/kt * dfydX(nt-1)));
    
    
    dfxdY(j) = mu*N(nt) * ( Y   )
