% =====================================================================
% deformed_laplacian_perturbation.m
% ---------------------------------------------------------------------
% Numerical companion to the paper
%   "Perturbation bounds and error analysis for eigenvalues of the
%    deformed Laplacian matrix"
%
% It reproduces the numerical examples of Section 6 and verifies:
%   Example 1 (path P8)     : Thm 3 (single edge), Prop 10 (derivative),
%                            Prop 9 (Davis--Kahan bound on K8)
%   Sharpness and structure : Thm 5 (K_n - e -> K_n attains |s|+s^2),
%                            Thm 4 (monotonicity threshold |s| >= 1),
%                            Prop 2 (positivity thresholds)
%   Example 2 (star S8)     : Thm 6 (k edges, Hausdorff distance),
%                            equality in Hong's bound and Lemma 7
%   Example 3 (ER G(10,.5)) : Prop 12 (condition numbers),
%                            Thm 13 (residual-based backward error)
%
% Running this script prints all tables of the paper and saves the five
% figures as Figure_1..Figure_5 in both EPS (vector source format) and
% PNG in the current folder.  The LaTeX file includes Figure_N.pdf,
% obtained from Figure_N.eps via ghostscript (pdfwrite, -dEPSCrop).
%
% Run:  simply execute this script in MATLAB (R2016b or later).
% =====================================================================

clc; clear; close all;

% deformed Laplacian  M_G(s) = I - s A + s^2 (D - I)
MG = @(A,D,s,n) eye(n) - s*A + s^2*(D - eye(n));

fprintf('=====================================================\n');
fprintf(' Example 1: a tree (the path P8)\n');
fprintf('=====================================================\n');

%% -----------------------------------------------------------------
%% Example 1, part A: single-edge deletion (Theorem 3)
%% -----------------------------------------------------------------
n = 8;
A = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);   % path adjacency
D = diag(sum(A,1));
% delete the middle edge {4,5}
u = 4; v = 5;
A2 = A; A2(u,v) = 0; A2(v,u) = 0; D2 = diag(sum(A2,1));

s_grid = [-2 -1.2 -0.5 0.3 0.75 1.0 1.5 2.0];
fprintf('\nTable 1: single-edge deletion {4,5},  max|dlam| vs |s|+s^2\n');
fprintf('  %6s  %10s  %10s  %8s\n','s','actual','bound','ratio');
worst = 0;
for s = s_grid
    l1 = sort(eig(MG(A ,D ,s,n)));
    l2 = sort(eig(MG(A2,D2,s,n)));
    act = max(abs(l1-l2));  bnd = abs(s)+s^2;  r = act/bnd;
    worst = max(worst,r);
    fprintf('  %6.2f  %10.4f  %10.4f  %8.4f\n', s, act, bnd, r);
end
fprintf('  worst ratio = %.4f  (must be <= 1)\n', worst);

%% -----------------------------------------------------------------
%% Example 1, part B: derivative w.r.t. s (Proposition 10)
%% -----------------------------------------------------------------
s0 = 1.3;  h = 1e-6;
[V,w] = eig(MG(A,D,s0,n));  w = diag(w);  [w,ord] = sort(w);  V = V(:,ord);
gaps = diff(w);  i = find(gaps == max(gaps),1) + 1;      % simple eigenvalue
ui = V(:,i);
formDeriv = ui'*(-A + 2*s0*(D - eye(n)))*ui;
lamP = sort(eig(MG(A,D,s0+h,n)));
lamM = sort(eig(MG(A,D,s0-h,n)));
numDeriv = (lamP(i) - lamM(i))/(2*h);
fprintf('\nTable 2: derivative, i = %d, lambda_i = %.6f\n', i, w(i));
fprintf('  formula  u^T(-A+2s(D-I))u = %.10f\n', formDeriv);
fprintf('  finite difference         = %.10f\n', numDeriv);
fprintf('  absolute difference       = %.2e\n', abs(formDeriv-numDeriv));

%% -----------------------------------------------------------------
%% Example 1, part C: eigenvector bound on K8 (Proposition 9)
%% -----------------------------------------------------------------
fprintf('\nTable (Prop 9) on the complete graph K8:\n');
nK = 8;  AK = ones(nK) - eye(nK);  DK = diag(sum(AK,1));
s = 0.8;
wK = sort(eig(MG(AK,DK,s,nK)));  iK = 1;                 % isolated smallest
[Vt,wt] = eig(MG(AK,DK,s,nK));  wt = diag(wt);  [~,od] = sort(wt);  VK = Vt(:,od);
uiK = VK(:,iK);
delta = wK(iK+1) - wK(iK);                               % gap to the bulk
A4 = AK; A4(1,2) = 0; A4(2,1) = 0; D4 = diag(sum(A4,1));
w4 = sort(eig(MG(A4,D4,s,nK)));
j  = find(abs(w4 - wK(iK)) == min(abs(w4 - wK(iK))),1);
[V4t,w4t] = eig(MG(A4,D4,s,nK));  w4t = diag(w4t);  [~,od4] = sort(w4t);  V4 = V4t(:,od4);
utilde   = V4(:,j);
sinTheta = sqrt(max(0, 1 - (uiK'*utilde)^2));
normE    = norm(MG(A4,D4,s,nK) - MG(AK,DK,s,nK), 2);
bound4   = normE/delta;
fprintf('  s = %.2f, gap delta = %.4f, ||E|| = %.4f\n', s, delta, normE);
fprintf('  sin(theta) = %.4f,  bound = %.4f,  OK = %d\n', ...
        sinTheta, bound4, sinTheta <= bound4 + 1e-12);

%% -----------------------------------------------------------------
%% Example 1, part D: figure 1 (trajectories) and figure 2 (envelope)
%% -----------------------------------------------------------------
s_fine = linspace(-2,2,401);
LamP8  = zeros(n, numel(s_fine));
for t = 1:numel(s_fine)
    LamP8(:,t) = sort(eig(MG(A,D,s_fine(t),n)));
end

figure('Color','w','Position',[80 80 620 420]);
plot(s_fine, LamP8, 'LineWidth', 1.5); grid on;
xlabel('s'); ylabel('\lambda_i(s)');
title('Eigenvalues of M_G(s), path P_8');
legend(arrayfun(@(k)['\lambda_' num2str(k)],1:n,'UniformOutput',false), ...
       'Location','bestoutside');
set(gcf,'PaperUnits','points','PaperPosition',[0 0 620 420],'PaperSize',[620 420]);
print(gcf,'Figure_1.eps','-depsc');
print(gcf,'Figure_1.png','-dpng','-r150');

k   = 3;
lmax = LamP8(end,:);  lmin = LamP8(1,:);
env  = abs(s_fine)*min(k,sqrt(2*k)) + s_fine.^2*k;
figure('Color','w','Position',[80 80 620 420]);
% shade the band [lmin-b, lmax+b] itself: fills the interior and
% highlights the guaranteed containment region
patch([s_fine fliplr(s_fine)],[lmax+env fliplr(lmin-env)],[1 .93 .93],'EdgeColor','none');
hold on;
plot(s_fine, lmax, 'k', 'LineWidth', 2);
plot(s_fine, lmin, 'k-.', 'LineWidth', 1.2);
plot(s_fine, lmax+env, 'r--', 'LineWidth', 1.2);
plot(s_fine, lmin-env, 'r--', 'LineWidth', 1.2);
grid on; xlabel('s'); ylabel('eigenvalues / envelope');
ylim([1.05*min(lmin-env), 1.05*max(lmax+env)]);
title(['P_8: eigenvalue band [\lambda_{min}-b, \lambda_{max}+b], b = |s|min(k,sqrt(2k))+s^2k, k = ' num2str(k)]);
legend('\lambda_{max}', '\lambda_{min}', '\lambda_{max}+b', '\lambda_{min}-b', 'Location','best');
set(gcf,'PaperUnits','points','PaperPosition',[0 0 620 420],'PaperSize',[620 420]);
print(gcf,'Figure_2.eps','-depsc');
print(gcf,'Figure_2.png','-dpng','-r150');

fprintf('\n=====================================================\n');
fprintf(' Sharpness of the single-edge bound and structural thresholds\n');
fprintf('=====================================================\n');

fprintf('\nTable 3: sharpness of Thm 5 on K_n - e -> K_n, ratio to |s|+s^2\n');
fprintf('  %3s  %9s  %9s  %9s  %9s\n', 'n', 's=-1', 's=0.5', 's=1', 's=2');
for n = [4 6 10 60]
    Ak = ones(n) - eye(n);  Ag = Ak;  Ag(1,2) = 0;  Ag(2,1) = 0;
    Dk = diag(sum(Ak,1));    Dg = diag(sum(Ag,1));
    rr = zeros(1,4);  c = 1;
    for s = [-1 0.5 1 2]
        rr(c) = max(abs(sort(eig(MG(Ag,Dg,s,n))) - sort(eig(MG(Ak,Dk,s,n))))) / (abs(s)+s^2);
        c = c + 1;
    end
    fprintf('  %3d  %9.4f  %9.4f  %9.4f  %9.4f\n', n, rr(1), rr(2), rr(3), rr(4));
end
for n = [4 6 10 60]
    fprintf('  n=%2d  closed form at s=-1: %.4f\n', n, (sqrt(n^2+4*n-12)-n+2)/4);
end

fprintf('\nMonotonicity threshold (Thm 4), empty2 -> single edge:\n');
A0 = zeros(2); A1 = [0 1; 1 0];
for s = [0.5 1 2]
    D0 = diag(sum(A0,2)); D1 = diag(sum(A1,2));
    dd = sort(eig(MG(A1,D1,s,2))) - sort(eig(MG(A0,D0,s,2)));
    fprintf('  s=%.2f: increments = [%+.4f, %+.4f]\n', s, dd(1), dd(2));
end

fprintf('\nPositivity thresholds (Prop 2):\n');
AS = zeros(8); AS(1,2:8) = 1; AS = AS + AS';  DS = diag(sum(AS,2));
for s = [0.5 1 2]
    fprintf('  star S_8  s=%.1f: lambda_min = %+.4f\n', s, min(eig(MG(AS,DS,s,8))));
end
AK6 = ones(6) - eye(6); DK6 = diag(sum(AK6,1));
fprintf('  K_6 (5-regular): lambda_min(0.5) = %+.4f, lambda_min(2) = %+.4f\n', ...
    min(eig(MG(AK6,DK6,0.5,6))), min(eig(MG(AK6,DK6,2,6))));

fprintf('\n=====================================================\n');
fprintf(' Example 2: a highly non-regular tree (the star S8)\n');
fprintf('=====================================================\n');

%% -----------------------------------------------------------------
%% Example 2: k-edge deletion, Hausdorff distance (Theorem 6)
%% -----------------------------------------------------------------
n = 8;
A = zeros(n);  A(1,2:n) = 1;  A(2:n,1) = 1;  D = diag(sum(A,1));  % star
leaves = 2:n;
s_all  = [-2 -1 0.5 1 1.5 2];

fprintf('\nTable 4: delete k pendant edges, d_H vs |s|min{k,sqrt(2k-n''+1)}+s^2k\n');
fprintf('  %2s  %10s  %10s  %10s  %10s  %10s\n', ...
        'k','dH(s=1)','bnd','dH(s=2)','bnd','worst');
for kk = 1:4
    A3 = A;
    for t = 1:kk
        vv = leaves(t);  A3(1,vv) = 0;  A3(vv,1) = 0;
    end
    D3 = diag(sum(A3,1));
    dA = A3 - A;  npr = nnz(sum(abs(dA),2));   % n' = kk+1 vertices in the diff. graph
    fac = min(kk, sqrt(2*kk - npr + 1));       % sharpened adjacency factor = sqrt(kk)
    l1  = sort(eig(MG(A ,D ,1.0,n)));  l3  = sort(eig(MG(A3,D3,1.0,n)));
    l1b = sort(eig(MG(A ,D ,2.0,n)));  l3b = sort(eig(MG(A3,D3,2.0,n)));
    dH1 = max(abs(l1-l3));  dH2 = max(abs(l1b-l3b));
    worst = 0;
    for s = s_all
        la = sort(eig(MG(A ,D ,s,n)));
        lb = sort(eig(MG(A3,D3,s,n)));
        worst = max(worst, max(abs(la-lb))/(abs(s)*fac + s^2*kk));
    end
    fprintf('  %2d  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
        kk, dH1, fac+kk, dH2, 2*fac+4*kk, worst);
end
fprintf('  star family attains equality in Hong''s bound and in Lemma 7:\n');
for kk = 1:4
    A3 = A;
    for t = 1:kk, vv = leaves(t); A3(1,vv) = 0; A3(vv,1) = 0; end
    D3 = diag(sum(A3,1));
    dA = A3 - A;  dD = D3 - D;  npr = nnz(sum(abs(dA),2));
    fprintf('  k=%d: ||dA||_2^2+||dD||_2 = %.4f (2k=%d), sqrt(2k-np+1) = %.4f, ||dA||_2 = %.4f\n', ...
        kk, norm(dA,2)^2+norm(dD,2), 2*kk, sqrt(2*kk-npr+1), norm(dA,2));
end

%% -----------------------------------------------------------------
%% Example 2, figure 3: star spectra
%% -----------------------------------------------------------------
LamS = zeros(n, numel(s_fine));
for t = 1:numel(s_fine)
    LamS(:,t) = sort(eig(MG(A,D,s_fine(t),n)));
end
figure('Color','w','Position',[80 80 620 420]);
plot(s_fine, LamS, 'LineWidth', 1.5); grid on;
xlabel('s'); ylabel('\lambda_i(s)');
% tight y-limits (skip the tick-alignment padding) and the legend
% placed inside the empty middle region
ylim([min(LamS(:))-0.6, max(LamS(:))+1]);
title('Eigenvalues of M_G(s), star S_8');
legend(arrayfun(@(k)['\lambda_' num2str(k)],1:n,'UniformOutput',false), ...
       'Location','north');
set(gcf,'PaperUnits','points','PaperPosition',[0 0 620 420],'PaperSize',[620 420]);
print(gcf,'Figure_3.eps','-depsc');
print(gcf,'Figure_3.png','-dpng','-r150');

fprintf('\n=====================================================\n');
fprintf(' Example 3: ER G(10,0.5) - condition numbers & backward error\n');
fprintf('=====================================================\n');

%% -----------------------------------------------------------------
%% Example 3: condition numbers (Proposition 12)
%% -----------------------------------------------------------------
% Fixed realization (kept identical to the paper for reproducibility).
n = 10;
A = [0 0 0 1 1 0 1 0 0 1;
     0 0 1 1 0 0 0 0 0 0;
     0 1 0 1 1 0 1 0 0 0;
     1 1 1 0 0 1 1 1 0 1;
     1 0 1 0 0 0 1 0 0 0;
     0 0 0 1 0 0 0 1 0 0;
     1 0 1 1 1 0 0 1 0 1;
     0 0 0 1 0 1 1 0 0 1;
     0 0 0 0 0 0 0 0 0 1;
     1 0 0 1 0 0 1 1 1 0];
D     = diag(sum(A,1));
rhoA  = max(eig(A));
Delta = max(sum(A,1)) - 1;
fprintf('\n  degrees: %s\n', mat2str(sum(A,1)));
fprintf('  rho(A) = %.4f,  max_j(d_j-1) = %d\n', rhoA, Delta);

fprintf('\nTable 5: condition numbers of lambda_max\n');
fprintf('  %6s  %10s  %10s  %10s  %10s\n', 's','lam_max','kap_abs','bound','kap_rel');
for s = [-2.0 -1.5 -1.0 -0.5 0.5 1.0 1.5 2.0]
    [V,w] = eig(MG(A,D,s,n));  w = diag(w);  [w,od] = sort(w);  V = V(:,od);
    lam = w(n);  ue = V(:,n);
    kap_abs = abs(ue'*(-A + 2*s*(D-eye(n)))*ue);
    bnd     = rhoA + 2*abs(s)*Delta;
    kap_rel = kap_abs*abs(s)/abs(lam);
    fprintf('  %6.1f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            s, lam, kap_abs, bnd, kap_rel);
end

%% -----------------------------------------------------------------
%% Example 3, figure 4: condition numbers (two panels)
%% -----------------------------------------------------------------
BND1 = rhoA + 2*Delta;
% (a) all eigenvalues at s = 1
[V,w] = eig(MG(A,D,1.0,n));  w = diag(w);  [~,od] = sort(w);  V = V(:,od);
kapA = zeros(n,1);
for ii = 1:n
    kapA(ii) = abs(V(:,ii)'*(-A + 2*(D-eye(n)))*V(:,ii));
end
% (b) kappa_abs(s) for lambda_max (exclude the degenerate point s = 0)
ss = linspace(-2,2,401);  kapS = nan(size(ss));
for t = 1:numel(ss)
    if abs(ss(t)) < 1e-12, continue; end
    [Vv,wv] = eig(MG(A,D,ss(t),n));  wv = diag(wv);  [~,odv] = sort(wv);
    ue = Vv(:,odv(n));
    kapS(t) = abs(ue'*(-A + 2*ss(t)*(D-eye(n)))*ue);
end
bndS = rhoA + 2*abs(ss)*Delta;

figure('Color','w','Position',[80 80 1100 420]);
subplot(1,2,1);
bar(1:n, kapA, 0.65); hold on;
plot([0.5 n+0.5],[BND1 BND1],'r--','LineWidth',1.6);
grid on; xlabel('eigenvalue index i'); ylabel('\kappa_i^{abs}(1)');
title('(a) Condition numbers at s = 1 (Laplacian)');
legend('\kappa_i^{abs}(1)', 'bound \rho(A)+2\Delta_D', 'Location','northwest');
subplot(1,2,2);
plot(ss, kapS, 'b', 'LineWidth', 2); hold on;
plot(ss, bndS, 'r--', 'LineWidth', 1.5);
grid on; xlabel('s'); ylabel('absolute condition number');
title('(b) \kappa_{max}^{abs}(s) vs s');
legend('\kappa_{max}^{abs}(s)', '\rho(A)+2|s|\Delta_D', 'Location','northwest');
set(gcf,'PaperUnits','points','PaperPosition',[0 0 1100 420],'PaperSize',[1100 420]);
print(gcf,'Figure_4.eps','-depsc');
print(gcf,'Figure_4.png','-dpng','-r150');

%% -----------------------------------------------------------------
%% Example 3, part C: noise-perturbed eigenpair (Table 6) + figure 5
%% -----------------------------------------------------------------
fprintf('\nTable 6: noise-perturbed approximate eigenpair of lambda_max\n');
[V6,w6] = eig(MG(A,D,1.0,n));  w6 = diag(w6);  [~,o6] = sort(w6);
u6 = V6(:,o6(n));
v6 = ones(n,1) - (u6'*ones(n,1))*u6;  v6 = v6/norm(v6);  % fixed unit vector orthogonal to u6
epsgrid = [0.05 0.10 0.20 0.50 1.00];
resS1 = zeros(size(epsgrid));  resS2 = zeros(size(epsgrid));
kS1   = zeros(size(epsgrid));  kS2   = zeros(size(epsgrid));
fprintf('  %5s  %10s  %3s  %10s  %3s\n','eps','rho(s=1)','k','rho(s=2)','k');
for t = 1:numel(epsgrid)
    ut = u6 + epsgrid(t)*v6;  ut = ut/norm(ut);
    for sv = [1 2]
        Mv = MG(A,D,sv,n);
        lq = ut'*Mv*ut;
        rr = norm((Mv - lq*eye(n))*ut);
        tG = (sqrt(1+2*rr)-1)/sv;  kc = ceil(tG^2/2);
        if sv == 1, resS1(t) = rr;  kS1(t) = kc;
        else,       resS2(t) = rr;  kS2(t) = kc;  end
    end
    fprintf('  %5.2f  %10.4f  %3d  %10.4f  %3d\n', ...
            epsgrid(t), resS1(t), kS1(t), resS2(t), kS2(t));
end

figure('Color','w','Position',[80 80 1100 420]);
subplot(1,2,1);
plot(epsgrid, resS1, 'b-o', 'LineWidth', 1.6); hold on;
plot(epsgrid, resS2, 'r-s', 'LineWidth', 1.6);
grid on; xlabel('\epsilon'); ylabel('\rho_{res}');
title('(a) Residual \rho_{res} vs noise level');
legend('s=1','s=2','Location','northwest');
subplot(1,2,2);
plot(epsgrid, kS1, 'b-o', 'LineWidth', 1.6); hold on;
plot(epsgrid, kS2, 'r-s', 'LineWidth', 1.6);
grid on; xlabel('\epsilon'); ylabel('certified lower bound k');
title('(b) certified k vs noise \epsilon');
legend('s=1','s=2','Location','northwest');
set(gcf,'PaperUnits','points','PaperPosition',[0 0 1100 420],'PaperSize',[1100 420]);
print(gcf,'Figure_5.eps','-depsc');
print(gcf,'Figure_5.png','-dpng','-r150');

%% ------------------------------------------------------------------
%% Example 3, part D: residual from a genuine graph perturbation
%% ------------------------------------------------------------------
fprintf('\nGenuine 2-edge perturbation of G(10,0.5):\n');
A5 = A; A5(5,7)=0; A5(7,5)=0; A5(9,10)=0; A5(10,9)=0; D5 = diag(sum(A5,1));
for sv = [1.0 2.0]
    M0v = MG(A,D,sv,n);   Md5 = MG(A5,D5,sv,n);
    [Vg,wg] = eig(Md5); wg = diag(wg); [~,og] = sort(wg); ugg = Vg(:,og(n));
    lgg = wg(og(n));
    rG = norm((M0v - lgg*eye(n))*ugg);
    nE = norm(Md5 - M0v, 2);
    tG = (sqrt(1+2*rG)-1)/sv;  kc = ceil(tG^2/2);
    fprintf('  s=%.1f : rho_res = %.4f , ||E||_2 = %.4f , certified k >= %d , rho_res <= ||E||_2 : %d\n', sv, rG, nE, kc, rG<=nE);
end

fprintf('\nAll checks completed. Figures Figure_1..Figure_5 saved as .eps and .png.\n');