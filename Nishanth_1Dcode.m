clear;
KT = 0.5;                             %Temperature
n = 200000;                             %Number of e
d = n/100;                              %Length of box
%np=20000;
M = 50;                                %mass of p
mel = 50;                             %mass of e
me = ones(1,n)*0.51;                    %array of mass of e
P = -30:0.01:30;                        %momenta of e
f = @(p)exp(-sqrt(mel^2+p.^2)/KT);      %momentum distribution
n1 = integral(f,-inf,inf);
fn = (1/n1)*exp(-sqrt(mel^2+P.^2)/KT);  %normalized momentum distribution

Mo = randpdf(fn,P,[1,n]);               %Sampled momenta of e
Mo(n/2+1:n) = -Mo(1:n/2);               %Symmetrized momenta for e, to make it net zero
Vel = Mo./sqrt(mel^2+Mo.^2);            %Vel of e
X=d*rand(1,n);                          %Positions of e from 0 to d

%Initializing mass,velocity and position of Proton
me(n+1) = M;
Vel(n+1) = 0;
X(n+1) = d/2;

A = X';                                 %transpose of X,Vel,me
B = Vel';
Me = me';
C(:,1) = A;
C(:,2) = B;
C(:,3) = Me;
C = sortrows(C,1);
X = C(:,1);
Vel = C(:,2);
me = C(:,3);

backward_coll = 0;

for c = 1:5000000
    ind = find(me==M);                  % Index of Proton
    ti = (X(ind)-X)./(Vel-Vel(ind));
    lind = find(ti<0);
    
    for i=1:length(lind)                % Periodic boundary check
        if lind(i)>ind
            xbuf = X(lind(i)) - d;
        else
            xbuf = X(lind(i)) + d;
        end
        ti(lind(i)) = (X(ind)-xbuf)/(Vel(lind(i))-Vel(ind));
    end
%     for i=1:length(lind)                % Periodic boundary check
%         if lind(i)>ind
%             xbuf = -X(lind(i));         % CHECK THIS 
%             ti(lind(i)) = (X(ind)-xbuf)/(Vel(lind(i))-Vel(ind));
%         end
%         if lind(i)<ind && abs(Vel(lind(i)))>abs(Vel(ind))
%             xbuf = d+X(lind(i));
%             ti(lind(i)) = (X(ind)-xbuf)/(Vel(lind(i))-Vel(ind));
%         end
%     end
    
    M2 = min(ti(ti>0));
    ind2 = find(ti==M2);                % Index of colliding electron
    if length(ind2)>1
        disp('here');
    end
    
    if sign(Vel(ind))==sign(Vel(ind2)) && norm(Vel(ind2)) > norm(Vel(ind))
        backward_coll = backward_coll + 1;
    end
    
    Pt = me(ind)*Vel(ind)/sqrt(1-Vel(ind)^2) + me(ind2)*Vel(ind2)/sqrt(1-Vel(ind2)^2);
    Et = me(ind)/sqrt(1-Vel(ind)^2) + me(ind2)/sqrt(1-Vel(ind2)^2);
    
%     X = X+(M2)*Vel;                         %Updating position of all particles 
%     X = mod(X,d);
    V = Pt/Et;
    v1p = (Vel(ind)-V)/(1-Vel(ind)*V);
    v2p = (Vel(ind2)-V)/(1-Vel(ind2)*V);
%     V_e(c) = Vel(ind2);                     % Vel of colliding electron
    Vel(ind) = (-v1p+V)/(1-v1p*V);
    % Vel(ind2) = (-v2p+V)/(1-v2p*V);
    V_n(c) = Vel(ind);
    
    clear Vel;
    
    % Reinitializing electron momentum, velocity, position
    Mo = randpdf(fn,P,[1,n]);               %Sampled momenta of e
    Mo(n/2+1:n) = -Mo(1:n/2);               %Symmetrized momenta for e, to make it net zero
    Vel = Mo./sqrt(mel^2+Mo.^2);            %Vel of e
    X = d*rand(1,n);                        %Positions of e from 0 to d
    
    % Initializing mass, velocity and position of p
    me = ones(1,n)*0.51;                   %array of mass of e
    me(n+1) = M;
    Vel(n+1) = V_n(c);
    X(n+1) = d/2;
    
    A = X';                                 %transpose of X,Vel,me
    B = Vel';
    Me = me';
    C(:,1) = A;
    C(:,2) = B;
    C(:,3) = Me;
    C = sortrows(C,1);
    X = C(:,1);
    Vel = C(:,2);
    me = C(:,3);
end