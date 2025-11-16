%%Sabina Baken MATH519 Project
clear; close all; clc;
%%Ordinary model and its stability
syms x y real
a_val=1;
b_val=1;
c_val=1;
d_val=1;
a = a_val;
b = b_val;
c = c_val;
d = d_val;
prey = a*x-b*x*y;
predator = c*x*y-d*y;
M = [prey;predator];
J = jacobian(M, [x,y]);
disp('Jacobian matrix is:');
disp(J);
[x_eq,y_eq] = solve([prey==0, predator==0], [x,y]);
disp('Equilibrium points are:');
disp([x_eq,y_eq]);