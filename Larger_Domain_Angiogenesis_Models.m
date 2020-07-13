% function [BCModel_StalkCellDensity, BCModel_TipCellDensity, PillayModel_StalkCellDensity, PillayModel_TipCellDensity, Omega, TimeMesh, Autonomous_BCModel_StalkCellDensity, Autonomous_BCModel_TipCellDensity, Autonomous_PillayModel_StalkCellDensity, Autonomous_PillayModel_TipCellDensity, BC_WaveFront, Pillay_WaveFront, Autonomous_BC_WaveFront, Autonomous_Pillay_WaveFront] = Larger_Domain_Angiogenesis_Models
function [BCModel_StalkCellDensity, BCModel_TipCellDensity, PillayModel_StalkCellDensity, PillayModel_TipCellDensity, Omega, TimeMesh] = Larger_Domain_Angiogenesis_Models
% load('9apr2020_Plots_LargerDomain_AutonomousAndIndependent_Models_BaselineParameters.mat', 'BCModel_StalkCellDensity', 'BCModel_TipCellDensity');
load('8jul2020_PillayCAModel_Fixed_P_p1e-1_P_m1_k100_1000Realizations.mat', 'TC_ColumnAverage', 'EC_ColumnAverage');
LW = 'linewidth'; IN = 'interpreter'; LA = 'latex'; FS = 'fontsize';
% TC_InitCond = 8e-1*[ones(201,1); zeros(1800,1)];
% TC_InitCond = 1e-3*[ones(201,1); zeros(1800,1)];
% TC_InitCond = 1e0*[TC_ColumnAverage(:,33);zeros(1800,1)];
% EC_InitCond = 1e0*[EC_ColumnAverage(:,33); zeros(1800,1)];

TC_InitCond = 1e0*TC_ColumnAverage(:,33);
EC_InitCond = 1e0*EC_ColumnAverage(:,33);

% TC_InitCond = 1e0*BCModel_TipCellDensity(:, 1);
% EC_InitCond = 1e0*BCModel_StalkCellDensity(:, 1);

% TC_InitCond = [BCModel_TipCellDensity(:, 1); zeros(2000,1)];
% EC_InitCond = [BCModel_StalkCellDensity(:, 1); zeros(2000,1)];


[BCModel_StalkCellDensity, BCModel_TipCellDensity, ~, ~] = SnailTrail_1D_PDE(TC_InitCond, EC_InitCond);
[PillayModel_StalkCellDensity, PillayModel_TipCellDensity, Omega, TimeMesh] = Pillay_1D_Model(TC_InitCond, EC_InitCond);
% [Autonomous_BCModel_StalkCellDensity, Autonomous_BCModel_TipCellDensity, ~, ~] = Autonomous_SnailTrail_1D_PDE(TC_InitCond, EC_InitCond);
% [Autonomous_PillayModel_StalkCellDensity, Autonomous_PillayModel_TipCellDensity, Omega, TimeMesh] = Autonomous_Pillay_1D_Model(TC_InitCond, EC_InitCond);

% BC_WaveFront = zeros(length(TimeMesh), 1);
% Pillay_WaveFront = BC_WaveFront;
% 
% for i = 1:length(TimeMesh)
% BC_WaveFront(i) = trapz(Omega, Omega.*BCModel_TipCellDensity(:, i))./trapz(Omega, BCModel_TipCellDensity(:, i));
% Pillay_WaveFront(i) = trapz(Omega, Omega.*PillayModel_TipCellDensity(:,i))./trapz(Omega, PillayModel_TipCellDensity(:, i));
% % Autonomous_BC_WaveFront(i) = trapz(Omega, Omega.*Autonomous_BCModel_TipCellDensity(:, i))./trapz(Omega, Autonomous_BCModel_TipCellDensity(:, i));
% % Autonomous_Pillay_WaveFront(i) = trapz(Omega, Omega.*Autonomous_PillayModel_TipCellDensity(:,i))./trapz(Omega, Autonomous_PillayModel_TipCellDensity(:, i));
% end

save('8jul2020_LargerDomain_OriginalModels_ae1e-3.mat');


plotting = false;


if plotting
figure;
p1 = plot(Omega, PillayModel_TipCellDensity(:, 1:32:end), 'b', LW, 1);
hold on;
% p4 = plot(Omega, Autonomous_PillayModel_TipCellDensity(:, 1:32:end), '--b', LW, 1);
p2 = plot(Omega, BCModel_TipCellDensity(:, 1:32:end), 'r', LW, 1);
% p3 = plot(Omega, Autonomous_BCModel_TipCellDensity(:, 1:32:end), '--r', LW, 1);
xlabel('$x$', IN, LA, FS, 16); ylabel('$N(x,t)$', IN, LA, FS, 16);

figure;
p5 = plot(Omega, PillayModel_StalkCellDensity(:, 1:32:end), 'b', LW, 1);
hold on;
% p8 = plot(Omega, Autonomous_PillayModel_StalkCellDensity(:, 1:32:end), '--b', LW, 1);
p6 = plot(Omega, BCModel_StalkCellDensity(:, 1:32:end), 'r', LW, 1);
% p7 = plot(Omega, Autonomous_BCModel_StalkCellDensity(:, 1:32:end), '--r', LW, 1);
xlabel('$x$', IN, LA, FS, 16); ylabel('$E(x,t)$', IN, LA, FS, 16);

% figure;
% p9 = plot(TimeMesh, max(abs(PillayModel_TipCellDensity-BCModel_TipCellDensity)./max(max([PillayModel_TipCellDensity,BCModel_TipCellDensity], [],1), [], 1), [], 1), 'k', LW, 1);
% xlabel('$t$', IN, LA, FS, 16); 
% 
% figure;
% p10 = plot(TimeMesh, max(abs(PillayModel_StalkCellDensity-BCModel_StalkCellDensity)./max(max([PillayModel_StalkCellDensity,BCModel_StalkCellDensity], [],1), [], 1), [], 1), 'k', LW, 1);
% xlabel('$t$', IN, LA, FS, 16); 

% figure;
% p11 = plot(TimeMesh, max(abs(Autonomous_PillayModel_TipCellDensity-Autonomous_BCModel_TipCellDensity)./max(Autonomous_BCModel_TipCellDensity, [], 1), [], 1), '--k', LW, 1);
% xlabel('$t$', IN, LA, FS, 16); 
% 
% figure;
% p12 = plot(TimeMesh, max(abs(Autonomous_PillayModel_StalkCellDensity-Autonomous_BCModel_StalkCellDensity)./max(Autonomous_BCModel_StalkCellDensity, [], 1), [], 1), '--k', LW, 1);
% xlabel('$t$', IN, LA, FS, 16); 
% 
% figure;
% p13 = plot(TimeMesh, 1./max(PillayModel_TipCellDensity, [], 1), 'b', LW, 1);
% hold on;
% p16 = plot(TimeMesh, 1./max(Autonomous_PillayModel_TipCellDensity, [], 1), '--b', LW, 1);
% p14 = plot(TimeMesh, 1./max(BCModel_TipCellDensity, [], 1), 'r', LW, 1);
% p15 = plot(TimeMesh, 1./max(Autonomous_BCModel_TipCellDensity, [], 1), '--r', LW, 1);
% legend([p13(1), p16(1), p14(1), p15(1)], 'Original Pillay Model', 'Autonomous Pillay Model', 'Original Snail-Trail Model', 'Autonomous Snail-Trail Model', IN, LA, FS, 14);
% ylabel('$\frac{1}{\max_x(N(x,t)}$', IN, LA, FS, 16);
% xlabel('$t$', IN, LA, FS, 16); xlim([TimeMesh(1) TimeMesh(end)]);
% 
% 
% figure;
% p17 = plot(TimeMesh, max(PillayModel_TipCellDensity, [], 1), 'b', LW, 1);
% hold on;
% p20 = plot(TimeMesh, max(Autonomous_PillayModel_TipCellDensity, [], 1), '--b', LW, 1);
% p18 = plot(TimeMesh, max(BCModel_TipCellDensity, [], 1), 'r', LW, 1);
% p19 = plot(TimeMesh, max(Autonomous_BCModel_TipCellDensity, [], 1), '--r', LW, 1);
% legend([p13(1), p16(1), p14(1), p15(1)], 'Original Pillay Model', 'Autonomous Pillay Model', 'Original Snail-Trail Model', 'Autonomous Snail-Trail Model', IN, LA, FS, 14);
% ylabel('$\max_x N(x,t)$', IN, LA, FS, 16);
% xlabel('$t$', IN, LA, FS, 16); xlim([TimeMesh(1) TimeMesh(end)]);
% 
% figure;
% plot(TimeMesh, BC_WaveFront, 'r', LW, 1);
% hold on;
% plot(TimeMesh, Pillay_WaveFront, 'b', LW, 1);
% title('Original Models'); xlabel('$t$', IN, LA, FS, 16); ylabel('$X(t)$', IN, LA, FS, 16);
% 
% figure;
% plot(TimeMesh, Autonomous_BC_WaveFront, 'r', LW, 1);
% hold on;
% plot(TimeMesh, Autonomous_Pillay_WaveFront, 'b', LW, 1);
% title('Autonomous Models'); xlabel('$t$', IN, LA, FS, 16); ylabel('$X(t)$', IN, LA, FS, 16);
end % if plotting
% save('24mar2020_LargerDomain_AutonomousAndIndependent_Models_BaselineParameters_ae1_1e-3TimesTipCellInitCond_StalkCellInitCond0.mat');
end