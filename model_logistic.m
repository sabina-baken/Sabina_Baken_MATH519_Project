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
k_prey = 100;
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

tspan = [0 10];
initial_cond = [10, 10];
f = @(t, X) [a*X(1)*(1-X(1)/k_prey) - b*X(1)*X(2);
              c*X(1)*X(2) - d*X(2)];
[t, sol] = ode45(f, tspan, initial_cond);


figure;
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 2);
xlabel('Prey (x)');
ylabel('Predator (y)');
title('Phase Portrait');
grid on;


[x_q, y_q] = meshgrid(0:1:20, 0:1:20);
u = a*x_q.*(1 - x_q./k_prey) - b.*x_q.*y_q; 
v = c.*x_q.*y_q - d.*y_q;                  
figure;
hold on;
quiver(x_q, y_q, u, v, 'Color', [0.5 0.5 0.5], 'AutoScaleFactor', 2);
plot(sol(:,1), sol(:,2), 'b', 'LineWidth', 2);
idx = 1:round(length(sol)/20):length(sol)-1;
quiver(sol(idx,1), sol(idx,2), ...
       sol(idx+1,1)-sol(idx,1), sol(idx+1,2)-sol(idx,2), ...
       0, 'Color', 'b', 'MaxHeadSize', 80);
xlabel('Prey (x)');
ylabel('Predator (y)');
title('Phase Portrait with Direction Field');
axis tight;
grid on;

figure;
plot(t, sol(:,1), 'b', t, sol(:,2), 'r', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
legend('Prey', 'Predator');
title('Population vs Time');
grid on;



