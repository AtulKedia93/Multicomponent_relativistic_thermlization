% 2D plotting
figure('units','normalized','outerposition',[0 0 0.75 1]); hold;

h=histogram(Eout2,'normalization','pdf');
x=h.BinEdges;y=h.Values;
x(1)=[];

% plot(KE,MB1,'g','LineWidth',2);
% plot(KE,MJ11,'r','LineWidth',3);
% plot(KE,MJ12)
Plot2 = plot(KE,RFDM1,'r','LineWidth',3);
Plot1 = plot(KE,RFDM2,'k','LineWidth',3);

xlabel('Energy (MeV)'), ylabel('Distribution f(E)'),
lgd = legend([Plot1 Plot2 h],{'Electron distribution (Fermi-Dirac)', 'Maxwell Boltzmann distribution', 'Proton distribution from simulation'});
title(lgd,'Simulation delta temperature kT = 1 MeV')
lgd.Position(1) = 0.4; 
lgd.Position(2) = 0.6;
lgd.Position(3) = 0.3;
lgd.Position(4) = 0.18;

set(gca,'FontSize',30)
box on
axis ( [-KT*0.1 KT*3 0 max(y)*0.3] )

savefig('Plot_Delta_1D_for_Paper_2020Apr_KT1.fig')
saveas(gcf,'Plot_Delta_1D_for_Paper_2020Apr_KT1.png')