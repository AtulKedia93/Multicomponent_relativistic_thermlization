%Photon scattering (gamma)(g) with Klein-Nishina scattering
%Natural units 2-D case
%Energies and masses in MeV
%Neutrons and protons are refered to as particle 1 or n(=nucleon)
%and gammas as particle 2 or g(=gamma)
%L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame

%potentially change theta2 to theta1 in doppler shift. Rybicki-Lyghtman.
%radiative procesess in astrophysics

n=10000000;                                     %Number of scattering events
e1=1;                                           %Neutron energy
KT=0.1;                                         %Temperature when BBN gets over
m1=1;                                           %Masses of particles
m2=0;
c = 1;                                          %speed of light
%Ee = 0.51 : 0.001*KT : 20*KT;                  %electron energies
%KE = Ee - m2;                                  %Kinetic enegies for e
%MB = (1/(KT)).*exp(-KE/(KT));
%MJ = (KE+m2)./(KT*(m2+KT)).*exp(-KE/(KT));
%incorrectMJ =((KE+m2)./m2).*(1/(m2*KT)).*(1/besselk(2,m2/KT)).*exp(-(KE+m2)./KT);     %Incorrect MJ distribution
Eg = 0.01*KT : 0.01*KT : 20*KT;                 %photon energies
Fg = Eg./(exp(Eg/KT)-1)/0.01635;                %photon distribution in 2-Dim
% normalization = 1/6 * pi * KT^2
% Fg in |v| = Fg * m1 v = 0.5*939*v^2/(exp(0.5*939*v^2/KT)-1)/0.01635 * m1 * v
Vout(n)=0;                               %Velocities of n
Eout(n)=0;                               %Energies of n

i=0;
while i < n
    e2=randpdf(Fg,Eg,[1,1]);                    %random energy selected
    if isnan(e2)
        continue
    end
    i = i+1;                                    %loop parameter
    
    e2 = e2 + m2;
    
    theta1 = randi([0 36000],1,1)/100;          %incidence angle
    theta1 = deg2rad(theta1);
    v1 = ([c*(1-m1^2*c^4/e1^2)^(1/2),0]);       %incident velocities in L-frame
    v2 = ([cos(theta1)*c*(1-m2^2*c^4/e2^2)^(1/2),sin(theta1)*c*(1-m2^2*c^4/e2^2)^(1/2)]);
    g1 = 1/sqrt(1-v1(1)^2/c^2);                 %Gamma factor
    
    V1 = ([0,0]);                               %Boosting to particle 1's frame
    V2 = ([(v2(1)-v1(1))/(1-v2(1)*v1(1)/c^2),v2(2)/(1-v2(1)*v1(1)/c^2)/g1]);
    theta2 = atan2(V2(2),V2(1));                %incidence angle in 1's frame: relativistic
    
    E2 = e2*(g1*(1-v1(1)*cos(theta2)));         %Relativistic Doppler shift of energy
    
    ctheta3=-1:0.01:1;                          %scattering angles
    Egnew = E2./(1+E2/m1*(1-ctheta3));
    
	dzig=(Egnew/E2).^2.*(Egnew/E2 + E2./Egnew + ctheta3.^2-1).*sqrt(1-ctheta3.^2); %d(sigma)/d(theta) Klein-Nishina scattering
    T=randpdf(dzig,ctheta3,[1,1]);              %T = cos(rand angle selected)
    Egnew = E2/(1+E2/m1*(1-T));
    
    Pn = ([E2-Egnew*T,-Egnew*sqrt(1-T^2)]);
    alpha = atan2(Pn(2),Pn(1));                 %Calculating V of n in COM frame
    
    Vn = 1/sqrt(1/c^2 + m1^2/norm(Pn)^2);
    Vn = ([cos(theta2 - alpha)*Vn,sin(theta2 - alpha)*Vn]); %Vel of scattered nucleon in COM frame
    
    vn = ([(Vn(1)+v1(1))/(1+Vn(1)*v1(1)/c^2),Vn(2)/(1+Vn(1)*v1(1)/c^2)/g1]); %Boosting back to L-frame
    Vout(i) = norm(vn);
    Eout(i)= m1*c^2/sqrt(1-Vout(i)^2/c^2);
    e1 = Eout(i);
end
save photonKN_10million4.mat