% 2D and 3D regular

% h1 = figure('units','normalized','outerposition',[0 0 0.75 1]);
% hold on
% set(gca,'FontSize',20)
% xlabel('Energy (MeV)'), ylabel('Distribution f(E)'),
% plot(KE,RFDM1,'r','LineWidth',3);
% axis ( [-KT*0.2 KT*8 0 15] )
% v = VideoWriter('Movie_3D_KT_01_inject_GeV.avi');
% v.FrameRate = 10;
% open(v);
% n = i;
% i=1;
% while i < n
%     cla
%     plot(KE,RFDM1,'r','LineWidth',3);
%     histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.003)
%     strng = sprintf('Scatterings %.1E ',i);
%     lgd = legend('MB dist',strng);
%     title(lgd,'3D Simulation temperature kT = 0.1 MeV')
%     Frame = getframe(h1);
%     writeVideo(v,Frame);
%     i = i + 6^(floor(log10(i)));
% end
%     cla
%     i = n;
%     plot(KE,RFDM1,'r','LineWidth',3);
%     histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.03)
%     strng = sprintf('Scatterings %.1E ',i);
%     lgd = legend('MB dist',strng);
%     title(lgd,'3D Simulation temperature kT = 0.1 MeV')
%     Frame = getframe(h1); 
%     writeVideo(v,Frame);
% close(v);

%1D regular

% h1 = figure('units','normalized','outerposition',[0 0 0.75 1]);
% hold on
% set(gca,'FontSize',20)
% xlabel('Energy (MeV)'), ylabel('Distribution f(E)'),
% plot(KE,RFDM1,'r','LineWidth',3);
% axis ( [-KT*0.2 KT*3 0 5] )
% v = VideoWriter('Movie_1D_KT_1_Injection.avi');
% v.FrameRate = 10;
% open(v);
% n = i;
% i=1;
% while i < 2999999
%     cla
%     plot(KE,RFDM1,'r','LineWidth',3);
%     histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.03)
%     strng = sprintf('Scatterings %.1E ',i);
%     lgd = legend('MB dist',strng);
%     title(lgd,'1D Simulation temperature kT = 1 MeV')
%     Frame = getframe(h1);
%     writeVideo(v,Frame);
%     i = i + 6^(floor(log10(i)));
% end
%     cla
%     i = n-1;
%     plot(KE,RFDM1,'r','LineWidth',3);
%     histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.03)
%     strng = sprintf('Scatterings %.1E ',i);
%     lgd = legend('MB dist',strng);
%     title(lgd,'1D Simulation temperature kT = 1 MeV')
%     Frame = getframe(h1); 
%     writeVideo(v,Frame);
% close(v);

% 2D and 3D injecting

h1 = figure('units','normalized','outerposition',[0 0 0.75 1]);
hold on
set(gca,'FontSize',20,'XScale', 'log','YScale', 'log')
xlabel('Energy (MeV)'), ylabel('Distribution f(E)'),
axis ( [KT*0.1 20000 0.0001 25] )
% v = VideoWriter('Movie_3D_KT_01_inject_GeV.avi');
% v.FrameRate = 10;
% open(v);
% n = i;
i=10;
% while i < 2999999
%     cla
    plot(KE,RFDM1,'r','LineWidth',3);
    histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.03, 'EdgeColor', [0,0.5070,0.7410])
    strng = sprintf('No. of Scatterings = %.1E',i);
    lgd = legend('MB dist',strng);
    title(lgd,'3D Simulation temperature kT = 0.025 MeV')
    saveas(gcf,'Plot_for_Method_Paper_2020Oct_KT0025_inject_10GeV_10.png')
%     Frame = getframe(h1);
%     writeVideo(v,Frame);
%     i = i + 6^(floor(log10(i)));
% end
%     cla
%     i = n;
%     plot(KE,RFDM1,'r','LineWidth',3);
%     histogram(Eout2(1:i),'normalization','pdf', 'BinWidth', KT*0.03)
%     strng = sprintf('Scatterings %.1E ',i);
%     lgd = legend('MB dist',strng);
%     title(lgd,'3D Simulation temperature kT = 0.1 MeV')
%     Frame = getframe(h1); 
%     writeVideo(v,Frame);
% close(v);