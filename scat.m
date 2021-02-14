%syms E2p
Etot=9.395560111779342e+02;   %E total
m2=.51;   %mass of electron
m1=939;   %mass of proton
E2=Etot-m1; %energy of electron (including mass)in stationary nucleon frame

T =1.472289062331803e-05; %scattering angle of electron in the " " frame 
p0 = sqrt(E2^2-m2^2); %momentum of electron in this stationary nucleon frame
%eqn = (sqrt(E2p^2-m2^2)*2*(p0*T))^2==(m1^2-m2^2+p0^2+(2*Etot*E2p)-Etot^2)^2; %Eqn 31 , since it is much easier to code http://www.np.ph.bham.ac.uk/research_resources/programs/ckin/kinematics.pdf
a1 = (Etot^2-(p0*T)^2);
b1=-2*Etot*(m2^2 +m1*E2);
c1= m1^2*E2^2 + m2^4 + 2*m1*m2^2*E2 + m2^2*p0^2*T^2;
Ee = (-b1+sqrt(b1^2-4*a1*c1))/(2*a1);
if ~(sqrt(Ee^2-m2^2)*p0*T + m2^2 + m1*E2 == Ee*(m1 + E2))
    Ee = (-b1-sqrt(b1^2-4*a1*c1))/(2*a1);
    if Ee<0
        1
    end
end
%bar = solve(eqn,E2p);
%double(E2p)

%EE2p=solve(eqn,E2p); % E2 in the paper, energy of electron(including mass) after scatter in the frame where nucleon "WAS" stationary
%double(EE2p)
% E3p=ET-EE2p;  % E3, energy of nucleon in the same frame above, (recoil energy including mass)
% p2p=(EE2p^2-m2^2)^.5; %momentum of electron after scatter, same frame as above 
% p3p=(E3p^2-m1^2)^.5;%momentum of nucleon after scatter, same frame as above 
% phi=asin((p2p/p3p)*sin(theta)); %angle of recoil after for nucleon.
% double(radtodeg(phi))