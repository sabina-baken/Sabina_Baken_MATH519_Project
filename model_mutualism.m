%%Sabina Baken MATH519 Project
clear; close all; clc;
%%Mutualism model and it's stability
syms x y real
r_1_val=5;
r_2_val=2;
r_1 = r_1_val;
r_2 = r_2_val;
a_12 = 0.5;
a_21 = 3;
k_mut1 = 100;
k_mut2 = 100;
f1 = r_1*x*(1-(x-a_12*y)/k_mut1);
f2 = r_2*y*(1-(y-a_21*y)/k_mut2);
J = jacobian([f1;f2], [x,y]);
disp('Jacobian matrix is:');
disp(J);
[x_eq,y_eq] = solve([f1==0, f2==0], [x,y]);
disp('Equilibrium points are:');
disp([x_eq,y_eq]);
J_origin = subs(J, [x,y],[x_eq(1),y_eq(1)]);
disp('Jacobian matrix at origin is:');
disp(J_origin);
lambda_origin = eig(J_origin);
disp('Eigenvalues at origin equilibrium:');
disp(lambda_origin);
if all(real(lambda_origin) < 0)
    disp('Equilibrium is stable.');
elseif all(real(lambda_origin) > 0)
    disp('Equilibrium is unstable.');
elseif any(real(lambda_origin) < 0) && any(real(lambda_origin) > 0)
    disp('Equilibrium is an unstable saddle point.');
elseif all(real(lambda_origin) == 0)
    disp('Equilibrium is a center (might need another method).');
else
    disp('Mixed or degenerate case (might need another method).');
end

num_eq = length(x_eq);
for k = 2:num_eq
J_nontr = subs(J, [x,y],[x_eq(k),y_eq(k)]);
fprintf('Jacobian matrix at nontrivial equilibrium (x, y) = (%.2f, %.2f):\n', ...
        double(x_eq(k)), double(y_eq(k)));
lambda = eig(J_nontr);
disp('Eigenvalues at this equilibrium:');
disp(lambda);
if all(real(lambda) < 0)
    disp('Equilibrium is stable.');
elseif all(real(lambda) > 0)
    disp('Equilibrium is unstable.');
elseif any(real(lambda) < 0) && any(real(lambda) > 0)
    disp('Equilibrium is an unstable saddle point.');
elseif all(real(lambda) == 0)
    disp('Equilibrium is a center (might need another method).');
else
    disp('Mixed or degenerate case (might need another method).');
end
end

tspan = [0 50];
initial_cond = [50, 50];
f = @(t, X) [r_1*X(1)*(1-(X(1)-a_12*X(2))/k_mut1);
             r_2*X(2)*(1-(X(2)-a_21*X(1))/k_mut2);];
[t, sol] = ode45(f, tspan, initial_cond);

figure;
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 2);
xlabel('Mutual 1 (x)');
ylabel('Mutual 2 (y)');
title('Phase Portrait');
grid on;


[x_q, y_q] = meshgrid(0:1:110, 0:1:110);
u = r_1.*x_q.*(1-(x_q-a_12.*y_q)/k_mut1);
v = r_2.*y_q.*(1-(y_q-a_21.*y_q)/k_mut2);              
figure;
hold on;
quiver(x_q, y_q, u, v, 'Color', [0.5 0.5 0.5], 'AutoScaleFactor', 2);
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 2);
idx = 1:round(length(sol)/20):length(sol)-1;
quiver(sol(idx,1), sol(idx,2), ...
       sol(idx+1,1)-sol(idx,1), sol(idx+1,2)-sol(idx,2), ...
       0, 'Color', 'b', 'MaxHeadSize', 80);
xlabel('Mutual 1 (x)');
ylabel('Mutual 2 (y)');
title('Phase Portrait with Direction Field');
axis tight;
grid on;



figure;
plot(t, sol(:,1), 'b', t, sol(:,2), 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
legend('Mutual 1', 'Mutual 2');
title('Population vs Time');
grid on;
