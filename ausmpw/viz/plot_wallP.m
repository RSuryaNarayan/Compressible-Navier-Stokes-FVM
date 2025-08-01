clear all; close all; clc;

data1 = readmatrix("wallPressure_case1.txt");
x1 = data1(:, 1);
wallP1 = data1(:,2);

data2 = readmatrix("wallPressure_case2.txt");
x2 = data2(:, 1);
wallP2 = data2(:,2);

data3 = readmatrix("wallPressure_case3.txt");
x3 = data3(:, 1);
wallP3 = data3(:,2);

%compute analytical solution
gamma=1.4;
M=2;
p1=1/gamma;
theta = 10*pi/180;
Cp = 2*theta/sqrt(M^2-1);
pw1 = p1*(1 + 0.5 * gamma * M^2 * Cp);
pw2 = p1*(1 - 0.5 * gamma * M^2 * Cp);

x = x1;
wallPe = zeros(length(x1),1);
for i =1:length(x1)
    if (x(i)>=0.5)
        wallPe(i)=pw2;
    else
        wallPe(i)=pw1;
    end
end


plot(x1, wallP1,'k--','LineWidth',3);
hold on;
plot(x2, wallP2,'r-.','LineWidth',3);
plot(x3, wallP3,'g:','LineWidth',3);
plot(x, wallPe,'b-.','LineWidth',3)
hold off;

title("Wall Static Pressure","interpreter", "latex");
xlabel("$x$ [m]","Interpreter","latex");
ylabel("$P(x)|_{surface}$","Interpreter","latex");
legend("Case 1", "Case 2", "Case 3", "Exact","Interpreter", "latex")
xlim([0,1])

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gcf,'renderer','Painters')
set(gca,...
    "FontSize", 16, ...
    "FontName", "Computer Modern Roman");

exportgraphics(gca, "wallP.pdf")