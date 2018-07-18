
%function [T,N] = FricElem2D(U,Om,mu,kt,kn)

clear
% ------------------------------------------------------------------------
% Inputs: 
% ------------------------------------------------------------------------
% U = displacements in frequency domain 
% Om = frequency - continuation parameter - not explicitly required 
% mu = friction coefficient
% kt = tangential stiffness in 1D
% kn = normal stiffness

nStep = 100;
UX = [1; 1.0; 0.5; 0.05] ;
% uy = [1.0; 0.2; 0.05] ;
UZ = [0.5; 0.5; 0.2; 0.1] ;

ux = myInvFFT(UX,nStep);
ux = [ux; ux];
v = myInvFFT(UZ,nStep);
v = [v; v];

wex = 1;                % excitation freq
% nNode = 1;              % no. of nodes
% nH = 3;                 % no. of harmonics
% dpn = 3;                % no. of DOF per node
% nDOF = dpn * nNode;     % no. of total DOF in the system        

mu = 0.2;
N0 = 50;
kn = 50;
ktx = 50; 
% ktx = kt;
kty = 40;

N = max (N0 + kn * v, 0);

w = zeros(length(ux),1);
t = linspace(0,1,nStep*2);
figure(2000); 
subplot(1,3,1)
title('U')
subplot(1,3,2)
% legend('tangential','normal')
title('forces evolution')
subplot(1,3,3)
title(['N_0 = ' num2str(N0) ' N'] )
for nt = 1:length(ux)
    Coul(nt) = mu*N(nt); 
    if nt == 1
        % Predictor
        w(nt) = ux(nt);
    else
        % Corrector
        w(nt) = w(nt-1);
    end
    
    
    Tx(nt) = ktx*(ux(nt) - w(nt));
    
    if abs(Tx(nt)) > Coul(nt)
        Tx(nt) = sign(Tx(nt))*Coul(nt);
        w(nt) = ux(nt) - Tx(nt)/ktx;
    end
    
    
    figure(2000); hold on
    subplot(1,3,1)
    plot(t(nt),ux(nt),'ro'), hold on
    
    subplot(1,3,2), 
    plot(t(nt),Tx(nt),'ro'), hold on
    plot(t(nt),N(nt),'bo'), 
    
    subplot(1,3,3), hold on
    if nt < length(ux)/2
        plot(ux(nt),Tx(nt),'ro')
    else
        plot(ux(nt),Tx(nt),'bo')
    end

end
figure(2000)
subplot(1,3,2)
legend('tangential','normal')