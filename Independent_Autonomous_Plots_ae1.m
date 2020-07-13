function Rel_error = Independent_Autonomous_Plots_ae1
LW = 'linewidth'; IN = 'interpreter'; LA = 'latex'; FS = 'fontsize'; PS = 'position'; position = [31, 60, 550, 432];
FN = 'fontname'; AR = 'arial'; 
load(['/Users/duncanmartinson/Documents/MATLAB/Oxford/Angiogenesis Models/',...
    '8jul2020_LargerDomain_OriginalModels_ae1_Lambda1e-2_1e-3_TCInitCond.mat']);
original_plots = false; save_yes_no = false; bvp = true;
lambda = 0.01; mu = 160; D = 1e-3; chi = 0.4; eps = D*lambda/chi^2; a_e = 1;

U_BC = mu/lambda*BCModel_TipCellDensity;
U_Pillay = mu/lambda*PillayModel_TipCellDensity;
W_BC = mu/lambda*a_e*BCModel_StalkCellDensity;
W_Pillay = mu/lambda*a_e*PillayModel_StalkCellDensity;
X = lambda/chi*Omega;
Tau = lambda*(TimeMesh);

jj = length(Tau);

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
fplot(@(tau) max(U_BC(:,1),[],1).*(Tau(1)).^(0.5).*(tau).^(-0.5), [Tau(1), Tau(ii)], '.k', LW, 1);
text(1e-1, 1, ['Line of Best Fit: $\ln(\max_X u(X,\tau))= $', num2str(line(1), '%1.2f'),'$\ln(\tau) + $', num2str(line(2))], IN, LA, FS, 16);
if save_yes_no
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'fig');
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'png');
saveas(gcf, ['5_',filename,'_LogLogMaxTipCellsVsTau'], 'eps');
end

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

z = zeros(length(Omega), length(Tau));
U_BC_Rescaled = z;
U_Pillay_Rescaled = z;
s = speed(1);
[~, ind] = max(U_BC(:,1));
% [~, ind] = max(U_Pillay(:,65));
% X_0 = X(ind)-s*Tau(1);
X_0 = speed(2);
for i = 1:length(Tau)
z(:, i) = (X-s*(Tau(i))-X_0)./sqrt(eps*Tau(i));
U_BC_Rescaled(:,i) = U_BC(:,i).*Tau(i).^(0.5);
U_Pillay_Rescaled(:,i) = U_Pillay(:,i).*Tau(i).^(0.5);
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
fplot(@(z) max([U_Pillay_Rescaled(:, 65+64), U_BC_Rescaled(:, 65+64)]).*exp(-z.^2./4), [-10, 10], '-.k', LW, 1);
legend('$\tau = 0.6\lambda$', '$\tau = \lambda$', '$\tau = 1.4\lambda$', '$\tau = 1.8\lambda$', '$\tau = 2.2\lambda$', 'Self-Similar Solution', 'interpreter', 'latex', 'fontsize', 14);
xlim([-10, 10]);

Rel_error = zeros(length(65+64:65+64*11), 1);
z = linspace(-10, 10, 2001)';
figure;
l = @(z) Tau(65+64)^(0.5).*max(U_BC(:,65+64)).*exp(-z.^2./4);
for i = 65+64:65+64*11
ll = interp1(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), Tau(i)^(-0.5).*l(z), X);
Rel_error(i) = max( [max( abs(ll-U_BC(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), ll], [], 1)) ), max( abs(ll-U_Pillay(:,i))./max(max([U_BC(:,i), U_Pillay(:,i), ll], [], 1)) )] );
    if ~mod(i-1,64)
        plot(z.*sqrt(eps*Tau(i))+X_0+s*Tau(i), Tau(i)^(-0.5).*l(z), '-.k', 'linewidth', 1); hold on;
    end % if
end
Rel_error = nonzeros(Rel_error);
plot(X, U_Pillay(:, 65+64:64:65+64*11), 'b');
plot(X, U_BC(:, 65+64:64:65+64*11), '--r');
xlim([0 .1])
ylim([0 0.05])
set(gcf, 'position', position);
set(gca, 'fontname', 'arial', 'fontsize', 16);
xlabel('$X$', 'fontsize', 20, 'interpreter', 'latex');
ylabel('$u(X,\tau)$', 'fontsize', 20, 'interpreter', 'latex');
title('Exponential function');

figure;
plot(Tau(65+64:65+64*11), Rel_error);
title('Relative error with self-similar solution');
end