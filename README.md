deformed_laplacian_perturbation.m
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
