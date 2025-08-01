
data = readmatrix("plt_converged_case1.txt", Delimiter=',');

X = data(:,1);
Y = data(:,2);
F = data(:,end);

X = reshape(X,48,71);
Y = reshape(Y,48,71);
F = reshape(F,48,71);

contourf(X,Y,F,'ShowText','off','LineColor','none','LevelStep',0.001);
colormap(hot)
c=colorbar();
clim([1.6 2.4])
title("$M=\frac{\sqrt{(u^2+v^2)}}{c}$, Case 1 (First order)","Interpreter","latex")
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gcf,'renderer','Painters')
set(gca,...
    "FontSize", 16, ...
    "FontName", "Computer Modern Roman");
exportgraphics(gca, "M_C1.pdf")