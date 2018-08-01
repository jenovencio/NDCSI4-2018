function [F,ur,ID] = FricElem3D2ts(x,w)

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
% dof's per node are only 2, there is only one dof in tangential direction
% if dpn == 2
%     kt = kt(1,1);       
% end

N0 = 70;
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
        IDx(n) = 0;         % separation phase
        IDy(n) = 0;         
        uxr(n) = ux(n);
        uyr(n) = uy(n);
    else
        Coul(n) = mu*N(n); 
        Tx(n) = ktx *( ux(n) - uxr(n) );
        Ty(n) = kty *( uy(n) - uyr(n) );
        IDx(n) = 1;         % stick phase
        IDy(n) = 1;
        if abs(Tx(n)) > Coul(n)
            Tx(n) = sign(Tx(n)).*Coul(n);
            uxr(n) = ux(n) - Tx(n)/ktx;
            IDx(n) = 2;     % slip phase         
        end
        if abs(Ty(n)) > Coul(n)
            Ty(n) = sign(Ty(n)).*Coul(n);
            uyr(n) = uy(n) - Ty(n)/kty;
            IDy(n) = 2;     % slip phase  
        end
        % updated to remove the static normal load
        %F(:,n) = [Tx(n); Ty(n); N(n)-N0]; 
        F(:,n) = [Tx(n); Ty(n); N(n)];
        ur(:,n) = [uxr(n); uyr(n)];
        ID = [IDx(n); IDy(n)];
    end
end
end


