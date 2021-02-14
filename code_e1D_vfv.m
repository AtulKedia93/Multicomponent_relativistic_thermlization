% Complete Long 1-D code
% Natural units
% Energies and masses in MeV
% Protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;

n = 3000000;                                    %Number of scattering events
KT = 1;                                         %Temperature when BBN gets over
m1 = 939;                                       %Masses of particles
m2 = 0.51;
e1 = 2*m1;                                        %Neuclei energy we start with

v = -0.99999:0.0001:0.99999;
g_v = 1./sqrt(1-v.^2);

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
%     a = 0;
    g_a = 1/sqrt(1-a^2);
%     v_MJv = abs(v).*g_v.^3.*(1+v*a).^3.*exp(-g_v*g_a.*(1+v*a)*m2/KT);
%     v_MJv = abs(v).*g_v.^3.*exp(-g_a*g_v.*(1+v*a)*m2/KT)./(g_v.*sqrt(1-((v+a)/(1+v*a)).^2));
%     v_MJv = abs(v).*exp(-m2*(v+a).^2/(2*KT));
%     v_MJv = abs(v).*g_v.^3.*exp(-g_a*g_v.*(1+v*a)*m2/KT);
    v_MJv = abs(v).*g_v.^3./(exp(g_a*g_v.*(1+v*a)*m2/KT)+1);
%     plot(v,v_MJv)
    v2 = randpdf(v_MJv,v,[1,1]);
    if isnan(v2)
        continue
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
plot(KE,MB1)
plot(KE,MJ11)
plot(KE,MJ12)
x=h.BinEdges;y=h.Values;
x(1)=[];

% save data_e1D.mat