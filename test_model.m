model_name = "";   
clearvars -except model_name;
close all; clc;
fprintf("Choose model to test:\n");
fprintf("  1) Logistic Predator–Prey\n");
fprintf("  2) Competition Model\n");
fprintf("  3) Mutualism Model\n");
choice = input("Enter 1 / 2 / 3: ");
switch choice
    case 1
        model_name = "model_logistic.m";
    case 2
        model_name = "model_competition.m";
    case 3
        model_name = "model_mutualism.m";
    otherwise
        error("Invalid choice. Write 1, 2, or 3.");
end
try
    run(model_name);
    Xe_all = double([x_eq, y_eq]);   
    numE = size(Xe_all,1);
    isStable = false(numE,1);
    for k = 1:numE
        Jk = double(subs(J, {x,y}, {Xe_all(k,1), Xe_all(k,2)}));
        lambda_k = eig(Jk);
        if all(real(lambda_k) < 0)
            isStable(k) = true;
        end
    end
    if any(isStable)
        idx = find(isStable,1,'first');
        Xeq = Xe_all(idx,:);
        fprintf("Selected STABLE equilibrium #%d: (%.6g, %.6g)\n", idx, Xeq(1), Xeq(2));
    elseif numE >= 2
        Xeq = Xe_all(2,:);
        fprintf("No stable equilibrium found. Using NONTRIVIAL equilibrium (index 2): (%.6g, %.6g)\n", Xeq(1), Xeq(2));
    else
        Xeq = Xe_all(1,:);
        fprintf("No stable equilibrium found. Using ONLY equilibrium: (%.6g, %.6g)\n", Xeq(1), Xeq(2));
    end
    f_sym = @(X) [ subs(f1, {x,y}, {X(1), X(2)});
                   subs(f2, {x,y}, {X(1), X(2)}) ];
    f_num = @(t,X) double(f_sym(X));    
    
    perturb_scale = 1e-1;           
    X0_near = Xeq .* (1 + perturb_scale); 
    X0_far  = Xeq + max(1, norm(Xeq))*[1, 1];

    if any(~isfinite(Xeq)) || any(imag(Xeq)~=0)
        error('Selected equilibrium is not a valid real finite point: %s', mat2str(Xeq));
    end
    tspan = [0 500];
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
    
    try
        [t_near, sol_near] = ode45(f_num, tspan, X0_near, opts);
        [t_far,  sol_far ] = ode45(f_num, tspan, X0_far,  opts);
    catch ME
        rethrow(ME);
    end
    
    % Compute distance-to-equilibrium error
    err_near = sqrt( (sol_near(:,1) - Xeq(1)).^2 + (sol_near(:,2) - Xeq(2)).^2 );
    err_far  = sqrt( (sol_far(:,1)  - Xeq(1)).^2 + (sol_far(:,2)  - Xeq(2)).^2 );
    tol = 1e-3; 
    
    isStableConverged = @(err) all(err(end) < tol);
    
    if isStableConverged(err_near)
        fprintf("✅ Test 1 PASSED: trajectory from near converged.\n");
    else
        fprintf("❌ Test 1 FAILED: trajectory from near did not converge.\n");
    end

    if isStableConverged(err_far)
        fprintf("✅ Test 2 PASSED: trajectory from far converged.\n");
    else
        fprintf("❌ Test 2 FAILED: trajectory from far did not converge.\n");
    end

    % PLOT 1: Distance vs Time (log scale) 
    figure('Name',sprintf('Log Distance to equilibrium'),'NumberTitle','off');
    semilogy(t_near, err_near, '-b', 'LineWidth', 1.6); hold on;
    semilogy(t_far,  err_far,  '--r', 'LineWidth', 1.2);
    xlabel('Time'); ylabel('Distance ||X(t)-X^*|| (log scale)');
    legend('near IC','far IC','Location','northeast');
    title(sprintf('Log-distance to equilibrium (%.3g, %.3g)', Xeq(1), Xeq(2)));
    grid on;
    
    % PLOT 2: Phase plane with trajectory and equilibrium 
    figure('Name',sprintf('Phase portrait'),'NumberTitle','off');

    xmax = max([max(sol_far(:,1)), max(sol_near(:,1)), Xeq(1)]) * 1.2;
    xmin = min([min(sol_far(:,1)), min(sol_near(:,1)), Xeq(1)]) * 0.8;
    ymax = max([max(sol_far(:,2)), max(sol_near(:,2)), Xeq(2)]) * 1.2;
    ymin = min([min(sol_far(:,2)), min(sol_near(:,2)), Xeq(2)]) * 0.8;
    
    if xmax==xmin; xmax = xmin+1; end
    if ymax==ymin; ymax = ymin+1; end
    
    [xg, yg] = meshgrid(linspace(xmin, xmax, 20), linspace(ymin, ymax, 20));
    ug = zeros(size(xg)); vg = zeros(size(yg));
    for ii = 1:numel(xg)
        V = double([subs(f1,{x,y},{xg(ii), yg(ii)}); subs(f2,{x,y},{xg(ii), yg(ii)})]);
        ug(ii) = V(1); vg(ii) = V(2);
    end

    quiver(xg, yg, ug, vg, 'AutoScaleFactor', 0.8, 'Color', [0.6 0.6 0.6]); hold on;
    plot(sol_near(:,1), sol_near(:,2), '-b', 'LineWidth', 1.6);
    plot(sol_far(:,1),  sol_far(:,2),  '--r', 'LineWidth', 1.2);
    plot(Xeq(1), Xeq(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize',8);
    text(Xeq(1), Xeq(2), sprintf('  (%.3g, %.3g)', Xeq(1), Xeq(2)), 'FontSize',10);
    xlabel('First species (x)'); ylabel('Second species (y)');
    title(sprintf('Phase portrait and trajectories'));
    legend('vector field','trajectory (near)','trajectory (far)','selected equilibrium','Location','best');
    axis tight; grid on;
    
catch ME
    rethrow(ME);
end

