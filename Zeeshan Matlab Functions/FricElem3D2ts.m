function [F,ur] = FricElem3D2ts(x,w)

% --------- test data ----------------

%current time step
% x = [1; 2; 3; 4; 5; 6; 7; 8; 9];
% w = [1; 1; 1; 2; 2; 2; 3; 3; 3];

% ------------------------------------
dpn = 3;                % degrees of freedom per node
N_node = length(x)/dpn; % total number of nodes

kn  = 50;
ktx = 50;          % tangential stiffness in x-direction
kty = 50;           % tangential stiffness in y-direction
% dof's per node are only 2, there is only one dof in tangential direction
% if dpn == 2
%     kt = kt(1,1);       
% end

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
    Coul(n) = mu*N(n); 
    Tx(n) = ktx *( ux(n) - uxr(n) );
    Ty(n) = kty *( uy(n) - uyr(n) );

    if abs(Tx(n)) > Coul(n)
        Tx(n) = sign(Tx(n)).*Coul(n);
        uxr(n) = ux(n) - Tx(n)/ktx;
    end
    if abs(Ty(n)) > Coul(n)
        Ty(n) = sign(Ty(n)).*Coul(n);
        uyr(n) = uy(n) - Ty(n)/kty;
    end
    F(:,n) = [Tx(n); Ty(n); N(n)];
    ur(:,n) = [uxr(n); uyr(n)];
end

end


