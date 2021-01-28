function [Rel_error, Rel_error2, Rel_error3] = Independent_Autonomous_Plots
LW = 'linewidth'; IN = 'interpreter'; LA = 'latex'; FS = 'fontsize'; PS = 'position'; position = [31, 60, 550, 432];
FN = 'fontname'; AR = 'arial'; 
load(['./',...
            '9apr2020_Plots_LargerDomain_AutonomousAndIndependent_Models_BaselineParameters.mat']);

%          '11apr2020_Plots_LargerDomain_OriginalModels_Chi4e-2_BaselineInitCond.mat']);

    %     '9apr2020_Plots_LargerDomain_AutonomousAndIndependent_Models_Lambda10_WideTildeLambda5.mat']);
%         '11apr2020_Plots_LargerDomain_OriginalModels_D1e-1_BaselineInitCond.mat']);

original_plots = false; save_yes_no = false; bvp = false;
% lambda = 0.16; mu = 160; D = 1e-3; chi = 0.4; eps = D*lambda/chi^2; a_e = 0.0391; beta = a_e*mu/lambda;
lambda = 0.16; mu = 160; D = 1e-3; chi = 0.4; eps = sqrt(D*lambda)/chi; a_e = 0.0391; beta = a_e*mu/lambda;
C0 = 0; dcdx = 1;

filename = 'R_28may2020_LargerDomain_OriginalModels_KazTWTest_DcDx_5e-1_C0_0_BaselineInitCond';

U_BC = mu/lambda*BCModel_TipCellDensity;
U_Pillay = mu/lambda*PillayModel_TipCellDensity;
W_BC = mu/lambda*a_e*BCModel_StalkCellDensity;
W_Pillay = mu/lambda*a_e*PillayModel_StalkCellDensity;
% X = lambda/chi*Omega;
X = sqrt(lambda/D)*Omega;
Tau = lambda*(TimeMesh);

% jj = find(Tau == 9.8*lambda);
jj = length(Tau);

% Abs_Linf_TC = abs(trapz(X, U_BC, 1)-trapz(X, U_Pillay, 1))./max([trapz(X, U_BC, 1); trapz(X,U_Pillay, 1)], [], 1);
% Abs_Linf_EC = abs(trapz(X, W_BC, 1)-trapz(X, W_Pillay, 1))./max([trapz(X, W_BC, 1); trapz(X,W_Pillay, 1)], [], 1);
% Abs_Linf_EC = max(abs(W_BC(127:end,:)-W_Pillay(127:end,:)), [], 1)./max([max(W_BC(127:end,:),[],1);max(W_Pillay(127:end,:),[],1)], [], 1);
% 
% figure;
% plot(Tau, Abs_Linf_TC);
% figure;
% plot(Tau, Abs_Linf_EC);

% Rel_Linf_TC2 = max(abs(U_BC-U_Pillay)./max([max(U_BC,[],1);max(U_Pillay,[],1)], [], 1), [], 1);
% Rel_Linf_EC2 = max(abs(W_BC-W_Pillay)./max([max(W_BC,[],1);max(W_Pillay,[],1)], [], 1), [], 1);
% 
% figure;
% plot(Tau, Rel_Linf_TC2);
% figure;
% plot(Tau, Rel_Linf_EC2);

Rel_Linf_TC = max(abs(U_BC-U_Pillay), [], 1)./max([max(U_BC,[],1);max(U_Pillay,[],1)], [], 1);
Rel_Linf_EC = max(abs(W_BC(11:end,:)-W_Pillay(11:end,:)), [], 1)./max([max(W_BC(11:end,:),[],1);max(W_Pillay(11:end,:),[],1)], [], 1);

figure;
plot(X, U_BC(:, 1:128:jj), '--r', LW, 1);
hold on;
plot(X, U_Pillay(:, 1:128:jj), 'b', LW, 1);
xlabel('$X$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$u(X,\tau)$', 'interpreter', 'latex', 'fontsize', 16);
xlim([X(1) X(end)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 16);
annotation('textarrow', [0.2 0.3], [0.5 0.5]);
text(0, 0, '$\tau$', 'interpreter', 'latex', 'fontsize', 16);
if save_yes_no
saveas(gcf, ['1_',filename,'_TipCells'], 'fig');
saveas(gcf, ['1_',filename,'_TipCells'], 'png');
saveas(gcf, ['1_',filename,'_TipCells'], 'eps');
end % if save

figure;
plot(X, W_BC(:, 1:128:jj), '--r', LW, 1);
hold on;
plot(X, W_Pillay(:, 1:128:jj), 'b', LW, 1);
xlabel('$X$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$w(X,\tau)$', 'interpreter', 'latex', 'fontsize', 16);
xlim([X(1) X(end)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 16);
annotation('textarrow', [0.2 0.3], [0.5 0.5]);
text(0, 0, '$\tau$', 'interpreter', 'latex', 'fontsize', 16);
if save_yes_no
saveas(gcf, ['2_',filename,'_StalkCells'], 'fig');
saveas(gcf, ['2_',filename,'_StalkCells'], 'png');
saveas(gcf, ['2_',filename,'_StalkCells'], 'eps');
end % if save

figure;
plot(Tau, Rel_Linf_TC, 'k', 'linewidth', 1);
xlabel('$\tau$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$\frac{||u_{ST}(X,\tau)-u_{P}(X,\tau)||_{\infty}}{\max\{||u_{ST}(X,\tau)||_\infty, \ ||u_{P}(X,\tau)||_\infty\}}$', 'interpreter', 'latex', 'fontsize', 16);
xlim([Tau(1) Tau(jj)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 16);
if save_yes_no
saveas(gcf, ['3_',filename,'_RelLInfError_TipCells'], 'fig');
saveas(gcf, ['3_',filename,'_RelLInfError_TipCells'], 'png');
saveas(gcf, ['3_',filename,'_RelLInfError_TipCells'], 'eps');
end % if save

figure;
plot(Tau, Rel_Linf_EC, 'k', 'linewidth', 1);
xlabel('$\tau$', 'fontsize', 16, 'interpreter', 'latex');
ylabel('$\frac{||w_{ST}(X,\tau)-w_{P}(X,\tau)||_{\infty}}{\max\{||w_{ST}(X,\tau)||_\infty, \ ||w_{P}(X,\tau)||_\infty\}}$', 'interpreter', 'latex', 'fontsize', 16);
xlim([Tau(1) Tau(jj)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 16);
if save_yes_no
saveas(gcf, ['4_',filename,'_RelLInfError_StalkCells'], 'fig');
saveas(gcf, ['4_',filename,'_RelLInfError_StalkCells'], 'png');
saveas(gcf, ['4_',filename,'_RelLInfError_StalkCells'], 'eps');
end % if save

figure; ii = 321;
loglog(Tau, max(U_BC, [], 1), '--r', LW, 1);
hold on;
loglog(Tau, max(U_Pillay, [], 1), 'b', LW, 1);
xlim([Tau(1) Tau(jj)]); xlabel('$\tau$', IN, LA, FS, 16);
ylabel('$\max_X u(X,\tau)$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 16);
line = polyfit(log(Tau(1:ii)), log(max(U_BC(:,1:ii), [],1)), 1);
% res = max(U_BC(:,1:ii),[],1) - exp(polyval(line,log(Tau(1:ii))));
% disp(norm(res));
fplot(@(tau) exp(polyval(line,log(tau))), [Tau(1), Tau(ii)], '-.k', LW, 1);
fplot(@(tau) max(U_BC(:,1),[],1).*(Tau(1)).^(1).*(tau).^(-1), [Tau(1), Tau(ii)], '.k', LW, 1);
text(1e-1, 1, ['Line of Best Fit: $\ln(\max_X u(X,\tau))= $', num2str(line(1), '%1.2f'),'$\ln(\tau) + $', num2str(line(2))], IN, LA, FS, 16);
if save_yes_no
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'fig');
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'png');
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'eps');

save('27may2020_Plots_LargerDomain_OriginalModels_KazTWTest_DcDx_1e-2_C0_9e-1_BaselineInitCond.mat');
end % if save

figure;
set(gcf, 'position', position);
BC_WaveFront = zeros(size(Tau));
Pillay_WaveFront = BC_WaveFront;
for i = 1:length(Tau)
BC_WaveFront(i) = trapz(X, X.*U_BC(:,i), 1)./trapz(X, U_BC(:,i), 1);
Pillay_WaveFront(i) = trapz(X, X.*U_Pillay(:,i), 1)./trapz(X, U_Pillay(:,i), 1);
end
plot(Tau, Pillay_WaveFront, 'b', 'linewidth', 1);
hold on;
plot(Tau, BC_WaveFront, '--r', 'linewidth', 1);
% ii = length(Tau);
ii = 65+64*11;
% ii = 65+64*3;
speed = polyfit(Tau(1:ii), BC_WaveFront(1:ii), 1);
fplot(@(tau) speed(1).*tau + speed(2), [.2*lambda 20*lambda], '-.k', 'linewidth', 1);
xlabel('$\tau$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$X(\tau)$', 'fontsize', 20, 'interpreter', 'latex');

kk = length(Tau);
ybar = mean(BC_WaveFront(1:kk));
SStot = sum((BC_WaveFront(1:kk)-ybar).^2);
SSres = sum((speed(1).*Tau(1:kk) + speed(2)-BC_WaveFront(1:kk)).^2);
R2 = 1-SSres/SStot;
fprintf(['R2 value for wavefront data on full time interval: ']);
disp(R2);

text(Tau(1), BC_WaveFront(1), ['$R^2 = $ ',num2str(R2)], 'interpreter', 'latex', 'fontsize', 16);
text(Tau(ii), BC_WaveFront(ii), ['$\widetilde{X}(\tau) = $',num2str(speed(1)),'$\tau + $',num2str(speed(2))],'interpreter', 'latex', 'fontsize', 16);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlim([Tau(1) Tau(end)]);

% Relative error with expected speed
mean_initial_x = trapz(Omega, Omega.*BCModel_TipCellDensity(:,1))./trapz(Omega, BCModel_TipCellDensity(:,1));
c = C0+dcdx*mean_initial_x;

fprintf(['c_0 = ',num2str(c),'\n']);
phi = (1+2*sqrt(c*eps));
fprintf(['Analytic Wavespeed is ',num2str(phi),'\n']);
fprintf('Relative error with the numerical result is ');
disp(abs(speed(1)-phi)./phi);

if save_yes_no
saveas(gcf, ['6_',filename,'_TipCellWaveFrontvsTime'], 'fig');
saveas(gcf, ['6_',filename,'_TipCellWaveFrontvsTime'], 'png');
saveas(gcf, ['6_',filename,'_TipCellWaveFrontvsTime'], 'eps');
end % if save

if bvp
z = zeros(length(Omega), length(Tau));
U_BC_Rescaled = z;
U_Pillay_Rescaled = z;
[max_u, ind] = max(U_BC(:,1));
s = speed(1);
% X_0 = X(ind)-s*Tau(1);
X_0 = speed(2);
for i = 1:length(Tau)
% z(:, i) = (X-s*(Tau(i))-X_0)./sqrt(eps*Tau(i));
z(:, i) = (X-s*(Tau(i))-X_0)./sqrt(Tau(i))./eps;
U_BC_Rescaled(:,i) = U_BC(:,i).*Tau(i).^(1);
U_Pillay_Rescaled(:,i) = U_Pillay(:,i).*Tau(i).^(1);
end
figure;
set(gcf, 'position', position);
colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; [0.4940, 0.1840, 0.5560]; [0.6350, 0.0780, 0.1840]];
for i = 1:5
plot(z(:,65+64*(i-1)), U_Pillay_Rescaled(:,65+64*(i-1)), 'linestyle', '-', 'color', colors(i,:)); hold on;
end
for i = 1:5
plot(z(:,65+64*(i-1)), U_BC_Rescaled(:,65+64*(i-1)), 'linestyle', '--', 'color', colors(i,:)); hold on;
end
legend('$\tau = 0.6\lambda$', '$\tau = \lambda$', '$\tau = 1.4\lambda$', '$\tau = 1.8\lambda$', '$\tau = 2.2\lambda$', 'interpreter', 'latex', 'fontsize', 14);
N = chebop(-100, 100);
N.init = chebfun('exp(-x.^2)', [-100, 100]);
N.op = @(x,y) diff(y,2)+x./2.*diff(y)+y-y.^2; N.lbc = 0; N.rbc = 0;
g = N\0;
delta = max_u*Tau(1)^(1)/max(g);
g = delta*g;
plot(g, '-.k', 'linewidth', 1, 'displayname', 'Self-Similar Solution');
xlim([-6 6]); ylim([0 0.45]);
xlabel('$z$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$\widetilde{U}(z)$', 'fontsize', 20, 'interpreter', 'latex');
set(gca, 'fontname', 'arial', 'fontsize', 16)

z = linspace(-100, 100, 20001)';
Rel_error = zeros(length(1:65+64*11), 1);
Rel_error2 = Rel_error;
Rel_error3 = Rel_error;
figure;
C1 = max(U_BC(:,1)).*Tau(1);
for i = 1:65+64*11
% gg = interp1(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), 1./Tau(i).*g(z), X);
gg = interp1(z.*eps.*sqrt(Tau(i))+X_0+s*Tau(i), 1./Tau(i).*g(z), X);
Rel_error(i) = max( [max( abs(gg-U_BC(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), gg], [], 1)) ), max( abs(gg-U_Pillay(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), gg], [], 1)) )] );
h = @(r) 1./(498*r.^2+1).*C1;
hh = 1./Tau(i).*h((X-s*Tau(i)-X_0)./sqrt(Tau(i)));
Rel_error2(i) = max( abs(hh-gg))./max(max([hh,gg], [], 1));
    if ~mod(i-1,64)
%         plot(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), 1./Tau(i).*g(z), '-.k', 'linewidth', 1); hold on;
        plot((z.*sqrt(Tau(i))+X_0+s*Tau(i)), 1./Tau(i).*g(z), '-.k', 'linewidth', 1); hold on;
    end % if
end
Rel_error = nonzeros(Rel_error);
plot(X, U_Pillay(:, 1:64:65+64*11), 'b');
plot(X, U_BC(:, 1:64:65+64*11), '--r');
xlim([0 1])
ylim([0 12])
set(gcf, 'position', position);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$u(X,\tau)$', 'fontsize', 20, 'interpreter', 'latex');

figure;
for i = 1:65+64*11
    h = @(r) 1./(498*r.^2+1).*max(U_BC(:,1)).*Tau(1);
    hh = 1./Tau(i).*h((X-s*Tau(i)-X_0)./sqrt(Tau(i)));
    if ~mod(i-1,64)
        plot(X, hh, '-.r', 'linewidth', 1); hold on;
        plot(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), 1./Tau(i).*g(z), '-k', 'linewidth', 1); hold on;
    end % if
end
xlim([0 1])
ylim([0 12])
set(gcf, 'position', position);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$u(X,\tau)$', 'fontsize', 20, 'interpreter', 'latex');
title('Full Self-Similar ODE solution vs. Approximate Asymptotic Solution to BVP');

figure;
l = @(z) Tau(1).*max(U_BC(:,1)).*exp(-z.^2./4);
for i = 1:65+64*11
ll = interp1(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), 1./Tau(i).*l(z), X);
Rel_error3(i) = max( [max( abs(ll-U_BC(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), ll], [], 1)) ), max( abs(ll-U_Pillay(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), ll], [], 1)) )] );
    if ~mod(i-1,64)
        plot(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), 1./Tau(i).*l(z), '-.k', 'linewidth', 1); hold on;
    end % if
end
plot(X, U_Pillay(:, 1:64:65+64*11), 'b');
plot(X, U_BC(:, 1:64:65+64*11), '--r');
xlim([0 1])
ylim([0 12])
set(gcf, 'position', position);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$u(X,\tau)$', 'fontsize', 20, 'interpreter', 'latex');
title('Exponential function');


[~, EC_BVP] = ode15s(@(t,e) ode(t,e,X,C1,s,X_0,eps,beta), Tau(1:65+64*11), W_BC(:,1));
EC_BVP = EC_BVP';
figure;
plot(X, EC_BVP(:, 1:32:65+64*11), '-.k', LW, 1);
hold on;
plot(X, W_Pillay(:, 1:32:65+64*11), 'b', LW, 1);
plot(X, W_BC(:, 1:32:65+64*11), '--r', LW,1);
xlim([0 1])
ylim([0 12])
set(gcf, 'position', position);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$w(X,\tau)$', 'fontsize', 20, 'interpreter', 'latex');
title('Exponential function');
end % if bvp-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if original_plots
figure;
p1 = plot(Omega, PillayModel_TipCellDensity(:, 1:64:end), 'b', LW, 1);
hold on;
p4 = plot(Omega, Autonomous_PillayModel_TipCellDensity(:, 1:64:end), '--b', LW, 1);
p2 = plot(Omega, BCModel_TipCellDensity(:, 1:64:end), 'r', LW, 1);
p3 = plot(Omega, Autonomous_BCModel_TipCellDensity(:, 1:64:end), '--r', LW, 1);
xlabel('$x$', IN, LA, FS, 16); ylabel('$N(x,t)$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['1_',filename,'_TipCells'], 'fig');
% saveas(gcf, ['1_',filename,'_TipCells'], 'png');
% saveas(gcf, ['1_',filename,'_TipCells'], 'eps');
% 
figure;
p5 = plot(Omega, PillayModel_StalkCellDensity(:, 1:64:end), 'b', LW, 1);
hold on;
p8 = plot(Omega, Autonomous_PillayModel_StalkCellDensity(:, 1:64:end), '--b', LW, 1);
p6 = plot(Omega, BCModel_StalkCellDensity(:, 1:64:end), 'r', LW, 1);
p7 = plot(Omega, Autonomous_BCModel_StalkCellDensity(:, 1:64:end), '--r', LW, 1);
xlabel('$x$', IN, LA, FS, 16); ylabel('$E(x,t)$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['2_',filename,'_StalkCells'], 'fig');
% saveas(gcf, ['2_',filename,'_StalkCells'], 'png');
% saveas(gcf, ['2_',filename,'_StalkCells'], 'eps');
% 
figure;
p9 = plot(TimeMesh, max(abs(PillayModel_TipCellDensity-BCModel_TipCellDensity)./max(BCModel_TipCellDensity, [], 1), [], 1), 'k', LW, 1);
xlabel('$t$', IN, LA, FS, 16); ylabel('$\frac{\max_x |N_P(x,t)-N_{ST}(x,t)|}{\max_x N_{ST}(x,t)}$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['3_',filename,'_RelDifference_OriginalModels_TipCells'], 'fig');
% saveas(gcf, ['3_',filename,'_RelDifference_OriginalModels_TipCells'], 'png');
% saveas(gcf, ['3_',filename,'_RelDifference_OriginalModels_TipCells'], 'eps');

figure;
p10 = plot(TimeMesh, max(abs(PillayModel_StalkCellDensity-BCModel_StalkCellDensity)./max(BCModel_StalkCellDensity, [], 1), [], 1), 'k', LW, 1);
xlabel('$t$', IN, LA, FS, 16); ylabel('$\frac{\max_x |E_P(x,t)-E_{ST}(x,t)|}{\max_x E_{ST}(x,t)}$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['4_',filename,'_RelDifference_OriginalModels_StalkCells'], 'fig');
% saveas(gcf, ['4_',filename,'_RelDifference_OriginalModels_StalkCells'], 'png');
% saveas(gcf, ['4_',filename,'_RelDifference_OriginalModels_StalkCells'], 'eps');

figure;
p11 = plot(TimeMesh, max(abs(Autonomous_PillayModel_TipCellDensity-Autonomous_BCModel_TipCellDensity)./max(Autonomous_BCModel_TipCellDensity, [], 1), [], 1), '--k', LW, 1);
xlabel('$t$', IN, LA, FS, 16); ylabel('$\frac{\max_x |N_{P}(x,t)-N_{ST}(x,t)|}{\max_x N_{ST}(x,t)}$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['5_',filename,'_RelDifference_AutonomousModels_TipCells'], 'fig');
% saveas(gcf, ['5_',filename,'_RelDifference_AutonomousModels_TipCells'], 'png');
% saveas(gcf, ['5_',filename,'_RelDifference_AutonomousModels_TipCells'], 'eps');

figure;
p12 = plot(TimeMesh, max(abs(Autonomous_PillayModel_StalkCellDensity-Autonomous_BCModel_StalkCellDensity)./max(Autonomous_BCModel_StalkCellDensity, [], 1), [], 1), '--k', LW, 1);
xlabel('$t$', IN, LA, FS, 16); ylabel('$\frac{\max_x |N_{P}(x,t)-N_{ST}(x,t)|}{\max_x N_{ST}(x,t)}$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['6_',filename,'_RelDifference_AutonomousModels_StalkCells'], 'fig');
% saveas(gcf, ['6_',filename,'_RelDifference_AutonomousModels_StalkCells'], 'png');
% saveas(gcf, ['6_',filename,'_RelDifference_AutonomousModels_StalkCells'], 'eps');

figure;
p13 = plot(TimeMesh, 1./max(PillayModel_TipCellDensity, [], 1), 'b', LW, 1);
hold on;
p16 = plot(TimeMesh, 1./max(Autonomous_PillayModel_TipCellDensity, [], 1), '--b', LW, 1);
p14 = plot(TimeMesh, 1./max(BCModel_TipCellDensity, [], 1), 'r', LW, 1);
p15 = plot(TimeMesh, 1./max(Autonomous_BCModel_TipCellDensity, [], 1), '--r', LW, 1);
legend([p13(1), p16(1), p14(1), p15(1)], 'Original Pillay Model', 'Autonomous Pillay Model', 'Original Snail-Trail Model', 'Autonomous Snail-Trail Model', IN, LA, FS, 14);
ylabel('$\frac{1}{\max_x(N(x,t)}$', IN, LA, FS, 16);
xlabel('$t$', IN, LA, FS, 16); xlim([TimeMesh(1) TimeMesh(end)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['7_',filename,'_InverseMax_TipCells'], 'fig');
% saveas(gcf, ['7_',filename,'_InverseMax_TipCells'], 'png');
% saveas(gcf, ['7_',filename,'_InverseMax_TipCells'], 'eps');


figure;
p13 = plot(TimeMesh, max(PillayModel_TipCellDensity, [], 1), 'b', LW, 1);
hold on;
p16 = plot(TimeMesh, max(Autonomous_PillayModel_TipCellDensity, [], 1), '--b', LW, 1);
p14 = plot(TimeMesh, max(BCModel_TipCellDensity, [], 1), 'r', LW, 1);
p15 = plot(TimeMesh, max(Autonomous_BCModel_TipCellDensity, [], 1), '--r', LW, 1);
legend([p13(1), p16(1), p14(1), p15(1)], 'Original Pillay Model', 'Autonomous Pillay Model', 'Original Snail-Trail Model', 'Autonomous Snail-Trail Model', IN, LA, FS, 14);
ylabel('$\max_x N(x,t)$', IN, LA, FS, 16);
xlabel('$t$', IN, LA, FS, 16); xlim([TimeMesh(1) TimeMesh(end)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['8_',filename,'_Max_TipCells'], 'fig');
% saveas(gcf, ['8_',filename,'_Max_TipCells'], 'png');
% saveas(gcf, ['8_',filename,'_Max_TipCells'], 'eps');

figure;
p21 = plot(TimeMesh, max(PillayModel_TipCellDensity, [], 1).^(-2/3), 'b', LW, 1);
hold on;
p24 = plot(TimeMesh, max(Autonomous_PillayModel_TipCellDensity, [], 1).^(-2/3), '--b', LW, 1);
p22 = plot(TimeMesh, max(BCModel_TipCellDensity, [], 1).^(-2/3), 'r', LW, 1);
p23 = plot(TimeMesh, max(Autonomous_BCModel_TipCellDensity, [], 1).^(-2/3), '--r', LW, 1);
legend([p21(1), p24(1), p22(1), p23(1)], 'Original Pillay Model', 'Autonomous Pillay Model', 'Original Snail-Trail Model', 'Autonomous Snail-Trail Model', IN, LA, FS, 14);
ylabel('$\max_x(N(x,t))^{-2/3}$', IN, LA, FS, 16);
xlabel('$t$', IN, LA, FS, 16); xlim([TimeMesh(1) TimeMesh(end)]);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['9_',filename,'_PowerofMax_TipCells'], 'fig');
% saveas(gcf, ['9_',filename,'_PowerofMax_TipCells'], 'png');
% saveas(gcf, ['9_',filename,'_PowerofMax_TipCells'], 'eps');

figure;
plot(TimeMesh, BC_WaveFront, 'r', LW, 1);
hold on;
plot(TimeMesh, Pillay_WaveFront, 'b', LW, 1);
title('Original Models'); xlabel('$t$', IN, LA, FS, 16); ylabel('$X(t)$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['10_',filename,'_WaveFront_OriginalModels'], 'fig');
% saveas(gcf, ['10_',filename,'_WaveFront_OriginalModels'], 'png');
% saveas(gcf, ['10_',filename,'_WaveFront_OriginalModels'], 'eps');

figure;
plot(TimeMesh, Autonomous_BC_WaveFront, 'r', LW, 1);
hold on;
plot(TimeMesh, Autonomous_Pillay_WaveFront, 'b', LW, 1);
title('Autonomous Models'); xlabel('$t$', IN, LA, FS, 16); ylabel('$X(t)$', IN, LA, FS, 16);
set(gcf, PS, position); set(gca, FN, AR, FS, 14);
% saveas(gcf, ['11_',filename,'_WaveFront_AutonomousModels'], 'fig');
% saveas(gcf, ['11_',filename,'_WaveFront_AutonomousModels'], 'png');
% saveas(gcf, ['11_',filename,'_WaveFront_AutonomousModels'], 'eps');
end % if original_plots
end

%% ODEs
function dedt = ode(t, ~, X, C1, s, X_0, eps, beta)
        dedt = beta*C1./t.*exp(-( (X-s*t-X_0)./sqrt(t)./sqrt(eps) ).^2./4);
end % ode