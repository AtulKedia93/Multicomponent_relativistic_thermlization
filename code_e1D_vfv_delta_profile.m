% Complete Long 1-D code
% Natural units
% Energies and masses in MeV
% Protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;

n = 10000000;                                    %Number of scattering events
KT = 1;                                         %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = m1;                                        %Neuclei energy we start with

v1 = sqrt(KT^2+2*m2*KT)/(KT+m2);

intRFDM1 = integral(@(KE)(KE+m1)./sqrt(KE.^2+2*m1*KE)./(exp(KE/KT)+exp(-m1/KT)),0,inf);
intRFDM2 = integral(@(KE)(KE+m2)./sqrt(KE.^2+2*m2*KE)./(exp(KE/KT)+exp(-m2/KT)),0,inf);

KE = 0.001*KT : 0.005*KT : 10*KT;

MB1 = (pi*KT*KE).^(-1/2).*exp(-KE/(KT));
MJ11 = (KE+m1)./(besselk(1,m1/KT).*sqrt(KE.^2+2*m1*KE)*m1).*exp(-(KE+m1)./KT);
MJ12 = (KE+m2)./(besselk(1,m2/KT).*sqrt(KE.^2+2*m2*KE)*m2).*exp(-(KE+m2)./KT);

RFDM1 = 1/intRFDM1*(KE+m1)./sqrt(KE.^2+2*m1*KE)./(exp(KE/KT)+exp(-m1/KT));
RFDM2 = 1/intRFDM2*(KE+m2)./sqrt(KE.^2+2*m2*KE)./(exp(KE/KT)+exp(-m2/KT));

i=1;
while i < n
    a = sqrt(1-(m1/e1)^2);
%     a = 0.1;
    
    v2 = rand(1);
    v_1 = ((v1-a)/(1-a*v1));
    v_2 = ((-v1-a)/(1+a*v1));
    
%     g_a = 1/sqrt(1-a^2);
%     g_1 = 1/sqrt(1-v_1^2);
%     g_2 = 1/sqrt(1-v_2^2);
%     
%     prob_1 = abs(v_1)*g_1^3/(exp(g_a*g_1*(1+v_1*a)*m2/KT)+1);
%     prob_2 = abs(v_2)*g_2^3/(exp(g_a*g_2*(1+v_2*a)*m2/KT)+1);
%     
%     if v2 < prob_1/(prob_1 + prob_2)
    
    
    if v2 < abs(v_1)/(abs(v_1)+abs(v_2))
        v2 = v_1;
    else
        v2 = v_2;
    end
    
    Ptot = m2*v2/sqrt(1-v2^2);
    Etot = m1 + m2/sqrt(1-v2^2);
    Vcom = Ptot/Etot;
    
    Vn = 2*Vcom/(1+Vcom^2);                     %velocity of Nucleus after collision in Nucleus' pre-collision frame
    Vout(i) = (a+Vn)/(1+a*Vn);                  %Boosting back to L-frame the V of nucleus after collision
    
    Eout(i)= m1/sqrt(1-Vout(i)^2);
    e1 = Eout(i);
    i = i+1;                                    %loop parameter
end

% i=i-1;
Eout2=Eout-m1;
h=histogram(Eout2,'normalization','pdf');
hold;
Plot2 = plot(KE,RFDM1,'r','LineWidth',3);
Plot1 = plot(KE,RFDM2,'k','LineWidth',3);
x=h.BinEdges;y=h.Values;
x(1)=[];

% save data_e1D.mat