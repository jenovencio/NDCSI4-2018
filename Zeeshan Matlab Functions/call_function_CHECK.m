
clear

nStep = 100;
wex = 1;     
UX = [1; 1.0; 0.2; 0.05] ;
UX = [0; 1.0; 0; 0] ;
UY = [1; 1.0; 0.2; 0.05] ;
UY = [0; 2.0; 0; 0] ;
UZ = [0.5; 0.5; 0.1; 0.1] ;
UZ = [0; 1.0; 0; 0] ;


ux = myInvFFT(UX,nStep);

ux = [ux; ux; ux];
% ux = zeros(length(ux),1);
% uy = myInvFFT(UY,nStep);
t = linspace(0,1,nStep*1);
uy = UY(1) + UY(2)*cos( 1*2*pi()*wex*t + pi()/2 )+ ...
             UY(3)*cos( 2*2*pi()*wex*t + pi()/2 )+...
             UY(4)*cos( 3*2*pi()*wex*t + pi()/2 );
uy = uy';
uy = [uy; uy; uy];
% uy = zeros(length(uy),1);
v = myInvFFT(UZ,nStep);
v = [v; v; v];
% v = zeros(length(v),1);
% N = N0*ones(1,length(v));
% w1 = zeros(length(ux),1);
% w2 = zeros(length(uy),1);
w = zeros(2,length(ux));
u = [ux'; uy'];
% T = zeros(2,length(ux));


for nt = 1:length(ux)
        if nt == 1
            % Predictor
            w(:,nt) = u(:,nt);
        else
            % Corrector
            w(:,nt) = w(:,nt-1);
        end
        
    U = [u(:,nt); v(nt)];
    W = [w(:,nt); 0];
    [f(:,nt),w(:,nt)] = FricElem3D2tsCoup(U,W);
    
%     [Jac] = JacFric3D2tsCoup(

nloop = ceil(nt/nStep);

% res1 = abs(f(1,100)-f(1,1));

switch nloop
    case 1
        
        %         figure(9011);
%         plot(U(1),f(1,nt),'ro'), hold on
% 
%         figure(9012);
%         plot(U(2),f(2,nt),'bo'), hold on
% 
%         figure(9013);
%         plot(f(1,nt),f(2,nt),'go'), hold on
    case 2
%         figure(9011);
%         plot(U(1),f(1,nt),'rs'), hold on
% 
%         figure(9012);
%         plot(U(2),f(2,nt),'bs'), hold on
% 
%         figure(9013);
%         plot(f(1,nt),f(2,nt),'gs'), hold on
    case 3
        figure(9011);
        plot(U(1),f(1,nt),'b.','MarkerSize',10), hold on
        
%         Us1(nt) = U(1);
%         fs1(nt) = f(1,nt);
        figure(9012);
        plot(U(2),f(2,nt),'b^'), hold on

        figure(9013);
        plot(f(1,nt),f(2,nt),'bx'), hold on
end
end


