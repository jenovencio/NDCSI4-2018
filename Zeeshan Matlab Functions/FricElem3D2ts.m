function [F,ur] = FricElem3D2ts(x,w)

% ***INPUT***
%       x = vector of displacements in time at t = tn for all nodes
%       w = vector of displacements in time at t = tn-1 for all nodes
% ***OUTPUT***
%       F = friction force in time
%       ur = updated displacements in time at t = tn which become
%       displacement at t = tn-1 in the next time step
%       ID = vector of contact state 
%           ID = 0, --------> separation
%           ID = 1, --------> stick
%           ID = 2, --------> slip

% --------- test data ----------------

%current time step
% x = [1; 2; 3; 4; 5; 6; 7; 8; 9];
% w = [1; 1; 1; 2; 2; 2; 3; 3; 3];

% ------------------------------------
dpn = 3;                % degrees of freedom per node
N_node = length(x)/dpn; % total number of nodes

kn  = 50;
ktx = 50;           % tangential stiffness in x-direction
kty = 50;           % tangential stiffness in y-direction

N0 = 20;
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
ID = zeros(2, N_node);          % contact condition identification

% arrange dof's in tangential plane 'u' and normal direction 'v'
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
    
    % Currently, the loop is computing the tangential force even when 
    % N = 0. It waits fot Tx and Ty to exceed Coulomb limit and only then 
    % the tangential forces becomes zero. Another aspect to consider is to
    % define sense of the normal load that defines separation. Say N > 0.
    % There should be another loop to check this so that the tangential 
    % forces  be returned as Tx = 0 and Ty = 0. 
    if N(n) <= 0
        Tx(n) = 0;
        Ty(n) = 0;
        ID(1,n) = 0;         % separation 
        ID(2,n) = 0;         % separation
        uxr(n) = ux(n);
        uyr(n) = uy(n);
        ur(:,n) = [uxr(n); uyr(n)];
    else
        Coul(n) = mu*N(n); 
        Tx(n) = ktx *( ux(n) - uxr(n) );
        Ty(n) = kty *( uy(n) - uyr(n) );
        ID(1,n) = 1;         % stick phase
        ID(2,n) = 1;
       
        if abs(Tx(n)) > Coul(n)
            Tx(n) = sign(Tx(n)).*Coul(n);
            uxr(n) = ux(n) - Tx(n)/ktx;
            ID(1,n) = 2;     % x slip         
        end
        
        if abs(Ty(n)) > Coul(n)
            Ty(n) = sign(Ty(n)).*Coul(n);
            uyr(n) = uy(n) - Ty(n)/kty;
            ID(2,n) = 2;     % y slip 
        end
        % updated to remove the static normal load
        %F(:,n) = [Tx(n); Ty(n); N(n)-N0]; 
        F(:,n) = [Tx(n); Ty(n); N(n)];
        ur(:,n) = [uxr(n); uyr(n)];
        id = ID(:,n);
        
    end
end
end


