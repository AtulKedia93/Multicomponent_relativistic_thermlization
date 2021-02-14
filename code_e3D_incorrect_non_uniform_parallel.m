syms thetatheta phiphi aa bb                    %solving equation for incoming angle transform. Used in if condition.
parpool('local',4);
parfor ii = 1:4

% Incorrect Fast Parallel 3-D code
% Natural units
% Energies and masses in MeV
% Neuclei and protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame

n = 1000;                                       %Number of scattering events
KT = 0.1;                                       %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = m1 + KT;                                   %Neuclei energy we start with
c = 1;                                          %speed of light
Ee = m2 : 0.001 : m2 + 10*KT;                   %electron energies
KE = Ee - m2;                                   %Kinetic enegies for e

% MB2 = (1/(KT)).*exp(-KE/(KT));
% MJ2 = (KE+m2)./(KT*(m2+KT)).*exp(-KE/(KT));
MB3 = 2*sqrt(KE/pi)*(1/KT)^(3/2).*exp(-KE/KT);
MJ3 = (KE+m2)/m2^2./(KT*besselk(2,m2/KT)).*sqrt(KE.^2+2*m2*KE).*exp(-(KE+m2)./KT);
% MB in |v| = MB2 * m1*v = (1/(KT))*exp(-0.5*m1*v^2/(KT))*m1*v
% MJ in |v| = MJ2 * m1*v =
% (0.5*m1*v^2+m2)/(KT*(m2+KT))*exp(-0.5*m1*v^2/(KT))*m1*v
% incorrectMJ2 =((KE+m2)./m2).*(1/(m2*KT)).*(1/besselk(2,m2/KT)).*exp(-(KE+m2)./KT);
% Eg = 0.01 : 0.001 : 2;                        %photon energies
% Fg = Eg.^2./(exp(Eg/KT)-1)/0.002404;          %photon distribution in 2-D

e2=randpdf(MJ3,KE,[1,n]);                       %random e- energies selected

ctheta3=-1:0.001:0.999;                         %scattering angle
% dzig=(sqrt(1-ctheta3.^2).*ctheta3./(1-ctheta3).^2);%Rutherford scattering
dzig=1./(1-ctheta3).^2;                         %Rutherford scattering
T=randpdf(dzig,ctheta3,[1,n]);                  %T = cos(rand angle selected)

eqn = tan(phiphi) == bb*sin(thetatheta)/(bb*cos(thetatheta)-aa);
sol = abs(solve(eqn,thetatheta));

Eout=[];Eout2=[];Vout=[];

i=1;
while i < n
    if isnan(e2(i)) || isnan(T(i))
        e2(i) = []; T(i) = [];
        n = n-1;
        continue
    end
    e2(i) = e2(i) + m2;
    
    a = sqrt(1-(m1/e1)^2);
    b = sqrt(1-(m2/e2(i))^2);
    if a < 0.999*b
        theta1= -pi:0.0001:pi;                  %-pi to pi because of looping effect at pi
        phi = atan2(b*sin(theta1),(b*cos(theta1)-a));
        f_phi = (a^2+b^2-2*a*b*cos(theta1))./(b^2-a*b*cos(theta1));
        phi=randpdf(f_phi,phi,[1,1]);
        
        solsol=double(subs(subs(subs(sol,aa,a),bb,b),phiphi,phi));
        if abs(abs(phi)-solsol(1))<abs(abs(phi)-solsol(2))
            theta1 = solsol(1);
        else
            theta1 = solsol(2);
        end
    else
        theta1 = pi;                            %incidence angle
    end
    
    v1 = ([c*a,0]);                             %incident velocities in L-frame
    v2 = ([cos(theta1)*c*b,sin(theta1)*c*b]);
    g1 = 1/sqrt(1-v1(1)^2/c^2);                 %Gamma factor
    
    V1 = ([0,0]);                               %Boosting to particle 1's frame
    V2 = ([(v2(1)-v1(1))/(1-v2(1)*v1(1)/c^2),v2(2)/(1-v2(1)*v1(1)/c^2)/g1]);
    V2norm = norm(V2);
    theta2 = atan2(V2(2),V2(1));                %incidence angle in 1's frame: relativistic
    G2 = 1/sqrt(1-V2norm^2/c^2);
    E2 = G2*m2*c^2;
    
    P2 = ([E2/c*sqrt((G2^2-1)/G2^2),0]);        %Momentum assuming 2 comes along x-axis
	p0 = norm(P2);                              %Momentum of 2
    
    Etot = E2 + m1*c^2;                         %Total energy of both particles before scattering
    
    a1=(Etot^2-(p0*T(i))^2);                    %Transferred energy calculation
    b1=-2*Etot*(m2^2 +m1*E2);
    c1= m1^2*E2^2 + m2^4 + 2*m1*m2^2*E2 + (m2*p0*T(i))^2;
    Ee1 = (-b1+sqrt(b1^2-4*a1*c1))/(2*a1);
    Ee2 = (-b1-sqrt(b1^2-4*a1*c1))/(2*a1);
    
    if (abs(Ee1*Etot - sqrt(Ee1^2-m2^2)*p0*T(i) -(m2^2 + m1*E2)) < abs(Ee2*Etot - sqrt(Ee2^2-m2^2)*p0*T(i) -(m2^2 + m1*E2)))
        Ee = Ee1;                               %Based on the solution of quadratic equation
    else
        Ee = Ee2;
    end
    if imag(Ee)~=0
        Ee = real(Ee);
    end
    if Ee>E2
        Ee = real(Ee2);
    end
    if Ee<m2
        continue
    end
    
    Pe = ([T(i)*sqrt(Ee^2 - m2^2*c^4)/c,sqrt(1-T(i)^2)*sqrt(Ee^2 - m2^2*c^4)/c]);
    Pn = P2 - Pe;
    alpha = atan2(Pn(2),Pn(1));                 %Calculating V of n in COM frame
    Vn = 1/sqrt(1/c^2 + m1^2/norm(Pn)^2);
    Vn = ([cos(theta2 - alpha)*Vn,sin(theta2 - alpha)*Vn]); %Vel of scattered nucleon in COM frame
    
    vn = ([(Vn(1)+v1(1))/(1+Vn(1)*v1(1)/c^2),Vn(2)/(1+Vn(1)*v1(1)/c^2)/g1]); %Boosting back to L-frame
    Vout(i) = norm(vn);
    Eout(i)= m1*c^2/sqrt(1-Vout(i)^2/c^2);
    e1 = Eout(i);
    i = i+1;                                    %loop parameter
end

i=i-1;
Eout2=Eout(1:i)-m1;
% h=histogram(Eout2,'normalization','pdf');
% x=h.BinEdges;y=h.Values;
% x(1)=[];
% hold;
% plot(KE,MJ3)
% plot(KE,MB3)
% plot(KE,MB3Dimensions)

%save e3D_incoming_non_uniform_100thsnd0.mat
parsave(sprintf('e3D_incoming_non_uniform_parallel_100thsnd%d.mat', ii),Eout,Vout);
end
delete(gcp)