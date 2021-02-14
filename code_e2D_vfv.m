% Complete Long 2-D code
% Natural units
% Energies and masses in MeV
% Protons are refered to as particle 1 or n(=nucleon)
% and electrons as particle 2 or e(=electron)
% L-frame = Lab frame, COM-frame = particle 1's frame initial rest frame
clear;

n = 10000000;                                   %Number of scattering events
KT = 1;                                      %Temperature when BBN gets over
m1 = 1;                                       %Masses of particles
m2 = 0.51;

e1 = m1;                                        %Neuclei energy we start with
v = -0.999:0.001:0.999;
l = length(v);
vx = repmat(v,l,1);
v = v';
vy = repmat(v,1,l);
g_v = 1./sqrt(1-vx.^2-vy.^2);
g_v(imag(g_v) ~= 0) = 0;                        %Remove velocities which are greater that c.

i=1;
while i < n
%     a = sqrt(1-(m1/e1)^2);
    a = 0.994;
    g_a = 1/sqrt(1-a^2);
%     v_MJv = sqrt(vx.^2+vy.^2).*g_v.^3.*exp(-g_a*g_v.*(1+vx*a)*m2/KT);
    v_MJv = sqrt(vx.^2+vy.^2).*g_v.^4.*exp(-g_a*g_v.*(1+vx*a)*m2/KT);
%     v_MJv = sqrt(vx.^2+vy.^2).*g_v.^2.*exp(-g_a*g_v.*(1+vx*a)*m2/KT).*(1+g_v.^2*vx.^2).*(1+g_v.^2*vy.^2);
%     plot(v,v_MJv)
    [v2(1),v2(2)] = pinky(vx,vy,v_MJv);
    
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
plot(KE,MJ2)
plot(KE,MB2)
x=h.BinEdges;y=h.Values;
x(1)=[];

save data_e2D.mat