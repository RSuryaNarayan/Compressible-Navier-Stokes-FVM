clear all; close all; clc;

data = readmatrix("residuals_case3.txt");

iters = data(:, 1);
res_rho = data(:, 2);
res_rhou = data(:, 3);
res_rhov = data(:, 4);
res_Et = data(:, 5);

loglog(iters, res_rho,'LineWidth',2);
hold on;
loglog(iters, res_rhou,'--','LineWidth',2);
loglog(iters, res_rhov,'-.','LineWidth',2);
loglog(iters, res_Et,':','LineWidth',2);
legend('$\rho$','$\rho u$', '$\rho v$','$E_t$','interpreter','latex');

xlabel("Iterations","Interpreter","latex");
ylabel("$RMS(U^{n+1}-U^n)$","Interpreter","latex");
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gcf,'renderer','Painters')
set(gca,...
    "FontSize", 16, ...
    "FontName", "Computer Modern Roman");

title("Residuals vs. Iterations, Case 3 (Second order, with Flux Limiter)","Interpreter","latex", "FontSize",10);

exportgraphics(gca, "res_C3.pdf")