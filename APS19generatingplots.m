Eout2=Eout(1:i)-m1;
h=histogram(Eout2,'normalization','pdf');
hold;
plot(KE,MJ2)
plot(KE,MB2)
x=h.BinEdges;y=h.Values;
x(1)=[];

(1/(KT))*exp(-0.5*939*x^2/(KT))*939*x


Distribution f(v/c)
Velocity (v/c)


(1/(KT))*exp(-0.5*m1*v^2/(KT))*m1*v

Fit M-Boltzmann at KT = 1.319 MeV
Simulation data

M-Boltzmann at KT = 1 MeV
Simulation data at Bkg T = 1 MeV