clear all; close all; clc;

data = readmatrix("fine.txt");
N = 71*48;
data = data(1:N,:);
X = data(:,1);
Y = data(:,2);

X = reshape(X,48,71);
Y = reshape(Y,48,71);

for i=1:48
    plot(X(i,:),Y(i,:),'-b');
    hold on;
end

for i=1:71
    plot(X(:,i),Y(:,i),'-k');
    hold on;
end

xlim([-0.5,1.5]);
ylim([0,1.5]);
pbaspect([2 1.5 1]);

title("Fine grid, $71\times48$","Interpreter","latex")
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gcf,'renderer','Painters')
set(gca,...
    "FontSize", 16, ...
    "FontName", "Computer Modern Roman");
exportgraphics(gcf,"fine.pdf")
