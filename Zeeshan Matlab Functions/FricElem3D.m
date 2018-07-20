
%function [T,N] = FricElem3D(U,Om,mu,ktx,kty,kn,est)

nStep = 100;
wex = 1;     
UX = [1; 1.0; 0.5; 0.05] ;
UY = [1; 1.0; 0.5; 0.05] ;
UZ = [0.5; 0.5; 0.2; 0.1] ;

ux = myInvFFT(UX,nStep);
ux = [ux; ux];
% uy = myInvFFT(UY,nStep);
t = linspace(0,1,nStep*1);
uy = UY(1) + UY(2)*cos( 1*2*pi()*wex*t + pi()/2 )+ ...
             UY(3)*cos( 2*2*pi()*wex*t + pi()/2 )+...
             UY(4)*cos( 3*2*pi()*wex*t + pi()/2 );
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
est = 1;

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
end



%%
if est == 1
    for nt = 1:length(ux)
        Coul(:,nt) = mu*N(nt); 
        if nt == 1
            % Predictor
            w(:,nt) = u(:,nt);
        else
            % Corrector
            w(:,nt) = w(:,nt-1);
        end
        
        T(:,nt) = kt *( u(:,nt) - w(:,nt) );

%         Tx(nt) = ktx*(ux(nt) - w1(nt));
%         Ty(nt) = kty*(uy(nt) - w2(nt));

        if abs(T(:,nt)) > Coul(:,nt)
            T(:,nt) = sign(T(:,nt)).*Coul(:,nt);
            w(:,nt) = u(:,nt) - kt\T(:,nt);
        end

%         if abs(Ty(nt)) > Coul(nt)
%             Ty(nt) = sign(Ty(nt))*Coul(nt);
%             w2(nt) = uy(nt) - Ty(nt)/kty;
%         end

    figure(est*1000+est*10);
plot(u(1,nt),T(1,nt),'ro'), hold on

figure(est*1001+est*10);
plot(u(2,nt),T(2,nt),'bo'), hold on

figure(est*1002+est*10);
plot(T(1,nt),T(2,nt),'go'), hold on
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
        T(nt) = sqrt(Tx(nt)^2+Ty(nt)^2);
        
        if T(nt) > Coul(nt) %'OR' abs(Tx(nt))>Coul(nt)
            T(nt) = Coul(nt);
            if abs(Tx) > Coul(nt)
                Tx(nt) = sign(Tx(nt))*Coul(nt);
                w1(nt) = ux(nt) - Tx(nt)/ktx;
                Ty(nt) = sign(Ty(nt))*sqrt(T(nt)^2-Tx(nt)^2);
                w2(nt) = uy(nt) - Ty(nt)/kty;
            
            else
                Ty(nt) = sign(Ty(nt))*Coul(nt);
                w2(nt) = uy(nt) - Ty(nt)/kty;
                Tx(nt) = sign(Tx(nt))*sqrt(T(nt)^2-Ty(nt)^2);
                w1(nt) = ux(nt) - Tx(nt)/ktx;
            end 
%             
        end

%         if T(nt) > Coul(nt) %& abs(Ty(nt)) > Coul(nt)
%             Ty(nt) = sign(Ty(nt))*Coul(nt);
%             w2(nt) = uy(nt) - Ty(nt)/kty;
%         end
    figure(est*1000+est);
plot(ux(nt),Tx(nt),'ro'), hold on

figure(est*1001+est);
plot(uy(nt),Ty(nt),'bo'), hold on

figure(est*1002+est);
plot(Tx(nt),Ty(nt),'go'), hold on
    end

%Tylp = Coul(
% figure(3001);
% plot(ux(nt),Tx(nt),'ro'), hold on
% 
% figure(3002);
% plot(uy(nt),Ty(nt),'bo'), hold on
% 
% figure(3003);
% plot(Tx(nt),Ty(nt),'go'), hold on
% plot(Tx(nt),Ty(nt),'ro'), hold on
% figure(3001);
% plot(ux(nt),Tx(nt),'ro'), hold on
% 
% figure(3001);
% 
end
