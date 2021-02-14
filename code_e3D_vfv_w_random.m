% Complete Long 3-D code
% Natural units
% Energies and masses in MeV
% Protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;
rng('shuffle');

n = 3000000;                                    %Number of scattering events
KT = 1;                                         %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = m1;                                        %Neuclei energy we start with
multiplier = 1;

% intRFD1  = integral(@(KE)(KE+m1).*sqrt(KE.^2+2*m1*KE)./(exp(KE/KT)+1),0,inf);
% intRFD2  = integral(@(KE)(KE+m2).*sqrt(KE.^2+2*m2*KE)./(exp(KE/KT)+1),0,inf);
intRFDM1 = integral(@(KE)(KE+m1).*sqrt(KE.^2+2*m1*KE)./(exp(KE/KT)+exp(-m1/KT)),0,inf);
intRFDM2 = integral(@(KE)(KE+m2).*sqrt(KE.^2+2*m2*KE)./(exp(KE/KT)+exp(-m2/KT)),0,inf);

KE = 0 : 0.01*KT : 10*KT;                        %Kinetic enegies for e

MB3 = 2*sqrt(KE/pi)*(1/KT)^(3/2).*exp(-KE/KT);
MJ31 = (KE+m1)./(m1^2*KT*besselk(2,m1/KT)).*sqrt(KE.^2+2*m1*KE).*exp(-(KE+m1)./KT);
MJ32 = (KE+m2)./(m2^2*KT*besselk(2,m2/KT)).*sqrt(KE.^2+2*m2*KE).*exp(-(KE+m2)./KT);
% RFD1 = 1/intRFD1*(KE+m1).*sqrt(KE.^2+2*m1*KE)./(1+exp(KE/KT));
% RFD2 = 1/intRFD2*(KE+m2).*sqrt(KE.^2+2*m2*KE)./(1+exp(KE/KT));
RFDM1 = 1/intRFDM1*(KE+m1).*sqrt(KE.^2+2*m1*KE)./(exp(-m1/KT)+exp(KE/KT));
RFDM2 = 1/intRFDM2*(KE+m2).*sqrt(KE.^2+2*m2*KE)./(exp(-m2/KT)+exp(KE/KT));

i=1;
while i < n
    a = sqrt(1-(m1/e1)^2);
%     a = 0.5;
    g_a = 1/sqrt(1-a^2);
    
    while a < 1
        rand_vel = rand(1,4);
        v2(1) = rand_vel(1)*2-1;
        v2(2) = rand_vel(2)*2-1;
        v2(3) = rand_vel(3)*2-1;
%         v2(1) = rand_vel(1)*cos(pi*rand_vel(2));
%         v2(2) = rand_vel(1)*sin(pi*rand_vel(2))*cos(2*pi*rand_vel(3));
%         v2(3) = rand_vel(1)*sin(pi*rand_vel(2))*sin(2*pi*rand_vel(3));
        g_v = 1/sqrt(1-v2(1)^2-v2(2)^2-v2(3)^2);
        if v2(1)^2+v2(2)^2+v2(3)^2 > 1
            continue
% chckd KT = 1.0 Mult = 14 | chckd KT = 0.5 Mult = 1.7 | chckd KT = 0.1 Mult = 3*10^-3 | chckd KT = 0.05 Mult = 10^-5 | chckd KT = 0.01 Mult = 10^-23
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^3*exp(-g_a*g_v*(1+v2(1)*a)*m2/KT)*(1+g_v^2*v2(1)^2)*(1+g_v^2*v2(2)^2)*(1+g_v^2*v2(3)^2)
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^5*exp(-g_a*g_v*(1+v2(1)*a)*m2/KT)
        elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^5/(exp((g_a*g_v*(1+v2(1)*a)*m2)/KT)+1)
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^5/(exp((g_a*g_v*(1+v2(1)*a)*m2-m2)/KT)+exp(-m2/KT))
            break
        end
    end
    
    v2(2) = sqrt(v2(2)^2+v2(3)^2);
    v2(3) = [];
    alpha = atan2(v2(2),v2(1));                 %Electron's velocity angle. Nucleus' velocity angle would be the same after collision.
    v2 = norm(v2);
    
    Ptot = m2*v2/sqrt(1-v2^2);
    Etot = m1 + m2/sqrt(1-v2^2);
    Vcom = Ptot/Etot;
                                                %Head-on collisions
    Vn = 2*Vcom/(1+Vcom^2);                     %Speed of Nucleus after collision in Nucleus' pre-collision frame
    vn1 = ([Vn*cos(alpha),Vn*sin(alpha)]);       %Velocity of Nucleus ''  ''  ''
    vn = ([(vn1(1)+a)/(1+a*vn1(1)) , vn1(2)/(g_a*(1+a*vn1(1)))]); %Velocity in Lab frame
    Vout(i) = norm(vn);                         %Speed in Lab frame
    
    Eout(i)= m1/sqrt(1-Vout(i)^2);
    e1 = Eout(i);
    i = i+1;                                    %loop parameter
end

i=i-1;
Eout2=Eout(1:i)-m1;
h=histogram(Eout2,'normalization','pdf');
hold;
% plot(KE,MB3)
% plot(KE,MJ31)
% plot(KE,MJ32)
% plot(KE,RFD1)
% plot(KE,RFD2)
plot(KE,RFDM1,'r','LineWidth',3);
plot(KE,RFDM2,'k','LineWidth',3);
x=h.BinEdges;y=h.Values;
x(1)=[];

% Adding non-relativistic Fermi-Dirac distribution for protons in the plot
% FD = 4*pi*sqrt(2*m1*KE)./(1+exp(KE/KT));

% FD = @(ke) sqrt(ke)./(1+exp(ke/KT));
% inte = integral(FD,0,inf);
% FD1= sqrt(KE)./(1+exp(KE/KT))/inte;
% FD1 = 2*(2+sqrt(2))*(1/KT)^(3/2)/(sqrt(pi)*zeta(3/2))*sqrt(KE)./(1+exp(KE/KT));
% plot(KE,FD1)

% save data_e3D.mat