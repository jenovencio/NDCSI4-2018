
%function [T,N] = FricElem3D(U,Om,mu,ktx,kty,kn,est)
clear

nStep = 100;
t = linspace(0,1,nStep*1);
wex = 1;     
UX = [1; 1.0; 0.5; 0.05] ;
UY = [1; 1.0; 0.5; 0.05] ;
UZ = [0.5; 0.5; 0.2; 0.1] ;

% ux = myInvFFT(UX,nStep);
ux = UX(2)*sin(2*pi*wex*t);
ux = ux';
ux = [ux; ux];

% uy = myInvFFT(UY,nStep);

% uy = UY(1) + UY(2)*cos( 1*2*pi()*wex*t + pi()/2 )+ ...
%              UY(3)*cos( 2*2*pi()*wex*t + pi()/2 )+...
%              UY(4)*cos( 3*2*pi()*wex*t + pi()/2 );

uy = UY(2)*sin(2*pi*wex*t+pi/2);
uy = uy';
uy = [uy; uy];
v = myInvFFT(UZ,nStep);
v = [v; v];


% nNode = 1;              % no. of nodes
% nH = 3;                 % no. of harmonics
% dpn = 3;                % no. of DOF per node
% nDOF = dpn * nNode;     % no. of total DOF in the system        

mu = 0.2;
N0 = 50;
kn = 50;
ktx = 50; 
kty = 50;

kt = [ktx 0; 0 kty];

% N = max (N0 + kn * v, 0);
N = N0*ones(1,length(v));
% w1 = zeros(length(ux),1);
% w2 = zeros(length(uy),1);
w = zeros(2,length(ux));
u = [ux'; uy'];
T = zeros(2,length(ux));
Coul = zeros(2,length(ux));
t = linspace(0,1,nStep*2);
est = 2;

if est == 1
    for nt = 1:length(ux)
        Coul(nt) = mu*N(nt); 
        if nt == 1
            % Predictor
            w1(nt) = ux(nt);
            w2(nt) = uy(nt);
        else
            % Corrector
            w1(nt) = w1(nt-1);
            w2(nt) = w2(nt-1);
        end


        Tx(nt) = ktx*(ux(nt) - w1(nt));
        Ty(nt) = kty*(uy(nt) - w2(nt));

        if abs(Tx(nt)) > Coul(nt)
            Tx(nt) = sign(Tx(nt)).*Coul(nt);
            w1(nt) = ux(nt) - Tx(nt)/ktx;
        end
        if abs(Ty(nt)) > Coul(nt)
            Ty(nt) = sign(Ty(nt)).*Coul(nt);
            w2(nt) = uy(nt) - Ty(nt)/kty;
        end

%         if abs(Ty(nt)) > Coul(nt)
%             Ty(nt) = sign(Ty(nt))*Coul(nt);
%             w2(nt) = uy(nt) - Ty(nt)/kty;
%         end

    figure(est*2000+est);
plot(ux(nt),Tx(nt),'ro'), hold on

figure(est*2001+est);
plot(uy(nt),Ty(nt),'bo'), hold on

figure(est*2002+est);
plot(Tx(nt),Ty(nt),'go'), hold on
    end
  

elseif est == 2
    
    for nt = 1:length(ux)
        Coul(nt) = mu*N(nt); 
        if nt == 1
            % Predictor
            w1(nt) = ux(nt);
            w2(nt) = uy(nt);
        else
            % Corrector
            w1(nt) = w1(nt-1);
            w2(nt) = w2(nt-1);
        end


        Tx(nt) = ktx*(ux(nt) - w1(nt));
        Ty(nt) = kty*(uy(nt) - w2(nt));
        Tnorm(:,nt) = [Tx(nt); Ty(nt)]./sqrt(Tx(nt)^2+Ty(nt)^2);
        
        if sqrt(Tx(nt)^2+Ty(nt)^2) > Coul(nt) %'OR' abs(Tx(nt))>Coul(nt)
            T(:,nt) = (Tnorm(:,nt))*Coul(nt);%*sign(Tnorm(:,nt));
                w1(nt) = ux(nt) - T(1,nt)/ktx;
                w2(nt) = uy(nt) - T(2,nt)/kty;
        
       
        end

%         if T(nt) > Coul(nt) %& abs(Ty(nt)) > Coul(nt)
%             Ty(nt) = sign(Ty(nt))*Coul(nt);
%             w2(nt) = uy(nt) - Ty(nt)/kty;
%         end
    figure(est*3000+est);
plot(ux(nt),T(1,nt),'ro'), hold on

figure(est*3001+est);
plot(uy(nt),T(2,nt),'bo'), hold on

figure(est*3002+est);
plot(T(1,nt),T(2,nt),'go'), hold on
    end
end


