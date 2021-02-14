% Complete Long 2-D code
% Natural units
% Energies and masses in MeV
% Neuclei and protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;
n = 20000000;                                   %Number of scattering events
KT = 0.06;                                         %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = m1;                                        %Neuclei energy we start with
c = 1;                                          %speed of light
Ee = m2 + 0.001*KT : 0.005*KT : m2 + 10*KT;                  %electron energies
KE = Ee - m2;                                   %Kinetic enegies for e

MB1 = (pi*KT*KE).^(-1/2).*exp(-KE/(KT));
MJ1 = (KE+m2)./(besselk(1,m2/KT).*sqrt(KE.^2+2*m2*KE)*m2).*exp(-(KE+m2)./KT);
% MB2 = (1/(KT)).*exp(-KE/(KT));
% MJ2 = (KE+m2)./(KT*(m2+KT)).*exp(-KE/(KT));
% MB3 = 2*sqrt(KE/pi)*(1/KT)^(3/2).*exp(-KE/KT);
% MJ3 = (KE+m2)/m2^2./(KT*besselk(2,m2/KT)).*sqrt(KE.^2+2*m2*KE).*exp(-(KE+m2)./KT);
% MB in |v| = MB2 * m1*v = (1/(KT))*exp(-0.5*m1*v^2/(KT))*m1*v
% MJ in |v| = MJ2 * m1*v =
% (0.5*m1*v^2+m2)/(KT*(m2+KT))*exp(-0.5*m1*v^2/(KT))*m1*v
% incorrectMJ2 =((KE+m2)./m2).*(1/(m2*KT)).*(1/besselk(2,m2/KT)).*exp(-(KE+m2)./KT);
% Eg = 0.01 : 0.001 : 2;                        %photon energies
% Fg = Eg.^2./(exp(Eg/KT)-1)/0.002404;          %photon distribution in 2-D


P=-30.0:.01:30.0;                   %momenta of e
f=@(p)exp(-sqrt(m2^2+p.^2)./KT);   %momentum distribution
n1=integral(f,-inf,inf);
fn=(1/n1).*exp(-sqrt(m2^2+P.^2)./KT);%normalized momentum distribution

Mo=randpdf(fn,P,[1,n]);


% e2 = randpdf(MJ1,KE,[1,n]);                     %random e- energies selected
% theta1 = (sign(rand(1,n)-0.5)+1)*pi/2;          %uniform incoming angles
backward_coll = 0;
vn=0;
i=1;
while i < n
%     if isnan(e2(i))
%         e2(i) = [];
%         n = n-1;
%         continue
%     end
%     e2(i) = e2(i) + m2;
%     
%     a = sqrt(1-(m1/e1)^2);
%     b = sqrt(1-(m2/e2(i))^2);
    
    v1 = vn;                                     %incident velocities in L-frame
%     v2 = cos(theta1(i))*b;
    v2 = Mo(i)/sqrt(m2^2+Mo(i)^2);
    if sign(v1)==sign(v2) && norm(v2) > norm(v1)
        backward_coll = backward_coll + 1;
    end
%     V1 = 0;                                     %Boosting to particle 1's frame
%     V2 = (v2-v1)/(1-v2*v1/c^2);
%     G2 = 1/sqrt(1-V2^2/c^2);
%     E2 = G2*m2*c^2;
%     
%     P2 = sign(V2)*E2/c*sqrt((G2^2-1)/G2^2);     %Momentum assuming 2 comes along x-axis
%     
%     Etot = E2 + m1*c^2;                         %Total energy of both particles before scattering
%     
%     Vcom=P2/Etot;
%     Vn = 2*Vcom/(1+Vcom^2);
%     Pn = m1*Vn/sqrt(1-Vn^2);
    
    Ptot = m1*v1/sqrt(1-v1^2) + m2*v2/sqrt(1-v2^2);
    Etot = m1/sqrt(1-v1^2) + m2/sqrt(1-v2^2);
    Vcom = Ptot/Etot;
    
    Vn = (v1-Vcom)/(1-Vcom*v1);
    
    vn = (-Vn+Vcom)/(1-Vn*Vcom/c^2); %Boosting back to L-frame
    
    Vout(i) = vn;
    Eout(i)= m1/sqrt(1-Vout(i)^2);
    e1 = Eout(i);
    i = i+1;                                    %loop parameter
end

i=i-1;
Eout2=Eout(1:i)-m1;
h=histogram(Eout2,'normalization','pdf');
hold;
plot(KE,MJ1,'k')
plot(KE,MB1,'g')
x=h.BinEdges;y=h.Values;
x(1)=[];

% save data_e1D.mat