%%Sabina Baken MATH519 Project
clear; close all; clc;
%%Ordinary model and it's stability
syms x y real
a_val=5;
b_val=2;
c_val=3;
d_val=5;
a = a_val;
b = b_val;
c = c_val;
d = d_val;
k_prey = 10;
prey = a*x*(1-x/k_prey)-b*x*y;
predator = c*x*y-d*y;
M = [prey;predator];
J = jacobian(M, [x,y]);
disp('Jacobian matrix is:');
disp(J);
[x_eq,y_eq] = solve([prey==0, predator==0], [x,y]);
disp('Equilibrium points are:');
disp([x_eq,y_eq]);
J_origin = subs(J, [x,y],[x_eq(1),y_eq(1)]);
disp('Jacobian matrix at origin is:');
disp(J_origin);

num_eq = length(x_eq);
for k = 2:num_eq
J_nontr = subs(J, [x,y],[x_eq(k),y_eq(k)]);
disp('Jacobian matrix at nontrivial equilibrium');
disp(J_nontr);
end
