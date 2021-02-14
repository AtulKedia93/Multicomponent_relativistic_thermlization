% 3D plotting
figure('units','normalized','outerposition',[0 0 0.75 1]); hold;

h=histogram(Eout2,'normalization','pdf');
x=h.BinEdges;y=h.Values;
x(1)=[];

% Plot2 = plot(KE,MB3,'r','LineWidth',3);
% plot(KE,MJ31)
% plot(KE,MJ32)
% plot(KE,RFD1)
% plot(KE,RFD2)
Plot2 = plot(KE,RFDM1,'r','LineWidth',3);
Plot1 = plot(KE,RFDM2,'k','LineWidth',3);

xlabel('Energy (MeV)'), ylabel('Distribution f(E)'),
lgd = legend([Plot1 Plot2 h],{'Electron distribution (Fermi-Dirac)', 'Maxwell Boltzmann distribution', 'Proton distribution from simulation'});
title(lgd,'Simulation temperature kT = 0.1 MeV')
lgd.Position(1) = 0.5;
lgd.Position(2) = 0.6;
lgd.Position(3) = 0.3;
lgd.Position(4) = 0.18;

set(gca,'FontSize',30)
box on
axis ( [-KT*0.2 KT*8 0 max(y)*1.1] )

% savefig('Plot_for_Paper_2020Feb_KT001.fig')
% saveas(gcf,'Plot_for_Paper_2020Feb_KT001.png')



% MB_for_residuals = 2*sqrt(x/pi)*(1/KT)^(3/2).*exp(-x/KT);;
% % plot(x, y)
% % plot(x, (y-MB_for_residuals)./(MB_for_residuals))
% x_min = 0.25*KT;
% x_max = 3*KT;
% 
% Numerator = 0;
% Denominator = 0;
% for i =1:1:length(x)
%     if x(i)> x_min && x(i) < x_max
%         Numerator = Numerator + ((y(i)-MB_for_residuals(i))/(MB_for_residuals(i)))^2;
%         Denominator = Denominator + (1/(MB_for_residuals(i)))^2;
%     end
% end
% result1 = sqrt(Numerator/Denominator);
% 
% Numerator = 0;
% Denominator = 0;
% for i =1:1:length(x)
%     if x(i)> x_min && x(i) < x_max
%         Numerator = Numerator + ((y(i)-MB_for_residuals(i))/(MB_for_residuals(i)))^2;
%         Denominator = Denominator + 1;
%     end
% end
% result2 = sqrt(Numerator/Denominator);
% 
% num2 = 0;
% Numerator = 0;
% Denominator = 0;
% for i =1:1:length(x)
%     if x(i)> x_min && x(i) < x_max
%         Numerator = Numerator + ((y(i)-MB_for_residuals(i)))^2;
%         Denominator = Denominator + 1;
%         num2 = num2 + MB_for_residuals(i);
%     end
% end
% num2 = num2/Denominator;
% result3 = sqrt(Numerator/Denominator);