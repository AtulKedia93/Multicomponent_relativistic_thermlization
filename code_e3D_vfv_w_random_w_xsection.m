% Complete Long 3-D code
% Natural units
% Energies and masses in MeV
% Protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;
rng('shuffle');

n = 5*10^6;                                    %Number of scattering events
KT = 0.1;                                         %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = m1;                                        %Neuclei energy we start with
multiplier = 10;

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
        elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^5/(exp((g_a*g_v*(1+v2(1)*a)*m2-m2)/KT)+exp(-m2/KT))
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^5/(exp((g_a*g_v*(1+v2(1)*a)*m2)/KT)+1)
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^3/(v2(1)^2+v2(2)^2+v2(3)^2)^2/(exp((g_a*g_v*(1+v2(1)*a)*m2)/KT)+1)
%         elseif rand_vel(4)*(multiplier) < sqrt(v2(1)^2+v2(2)^2+v2(3)^2)*g_v^3/(v2(1)^2+v2(2)^2+v2(3)^2)^2/(exp((g_a*g_v*(1+v2(1)*a)*m2 - m2)/KT)+exp(-m2/KT))
            break
        end
    end
    
    v2(2) = sqrt(v2(2)^2+v2(3)^2);
    v2(3) = [];
    theta = atan2(v2(2),v2(1));                 %Electron's velocity angle. Nucleus' velocity angle would be the same after collision.
    v2 = norm(v2);
    e2 = m2/sqrt(1-v2^2);
    p2 = v2*e2;
    
    Etot = m1 + e2;
    
    ctheta3=-0.999:0.001:-0.7;                         %scattering angle
    dzig=(1-v2^2*(1-ctheta3)./2).*sqrt(1-ctheta3.^2)./(1-ctheta3).^2;	%d(sigma)/d(theta) Mott(rutherford) scattering
%     dzig=(1-v2^2*(1-ctheta3)./2).*ctheta3.^2./(1-ctheta3).^2;                       %d(sigma)/d(theta) Mott(rutherford) scattering
    T=randpdf(dzig,ctheta3,[1,1]);              %T = cos(rand angle selected)
%     T = -1;
    
    a1=(Etot^2-(p2*T)^2);                    %Transferred energy calculation
    b1=-2*Etot*(m2^2 +m1*e2);
    c1= m1^2*e2^2 + m2^4 + 2*m1*m2^2*e2 + (m2*p2*T)^2;
    Ee1 = (-b1+sqrt(b1^2-4*a1*c1))/(2*a1);      %calculation based on http://www.np.ph.bham.ac.uk/research_resources/programs/ckin/kinematics.pdf
    Ee2 = (-b1-sqrt(b1^2-4*a1*c1))/(2*a1);      %They have a couple of typos/algebraic mistakes.
    
    if (abs(Ee1*Etot - sqrt(Ee1^2-m2^2)*p2*T -(m2^2 + m1*e2)) < abs(Ee2*Etot - sqrt(Ee2^2-m2^2)*p2*T -(m2^2 + m1*e2)))
        Ee = Ee1;                               %Based on the solution of quadratic equation
    else
        Ee = Ee2;
    end
    if imag(Ee)~=0
        Ee = abs(Ee);
    end
    if Ee>e2
        Ee = abs(Ee2);
    elseif Ee<m2
        Ee = abs(Ee1);
    end
    if (Ee>e2 || Ee<m2)
        continue
    end
    Pe = ([T*sqrt(Ee^2 - m2^2),sqrt(1-T^2)*sqrt(Ee^2 - m2^2)]);
    Pn = ([p2,0]) - Pe;
    
    alpha = atan2(Pn(2),Pn(1));                 %Calculating V of n in COM frame
    alpha = alpha*(2*round(rand)-1);            %Up and down symmetry
    Vn = 1/sqrt(1 + m1^2/norm(Pn)^2);
    
    vn1 = ([cos(theta - alpha)*Vn,sin(theta - alpha)*Vn]);      %Velocity of Nucleus after collision in Nucleus' pre-collision frame
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
% x=h.BinEdges;y=h.Values;
% x(1)=[];

% save data_e3D.mat