function optimized_code()

    % Initial parameters and constraints
    D0 = 5e-6;   % Initial guess for diffusion coefficient
    lb = 1e-8;   % Lower bound for D
    ub = 1e-4;   % Upper bound for D
    options = optimoptions('fmincon', 'Display', 'iter');  % Optimization options
    
    % Perform optimization
    D_opt = fmincon(@costfun, D0, [], [], [], [], lb, ub, [], options);
    
    % Print optimized D value
    disp(['Optimized D value: ' num2str(D_opt)]);
    
    % Plotting
    global Time Diffusion_Fraction
    L = 0.04;  % Length parameter
    
    % Calculate M_1(t) with optimized D
    syms n t
    f_1(t) = (8 / ((2*n+1)^2 * pi^2)) * exp((-D_opt * (2*n+1)^2 * pi^2 * t) / (4 * L^2));
    M_1(t) = 1 - symsum(f_1(t), n, 0, Inf);
    
    % Define time vector T in seconds
    T = 0:0.2:7200;
    
    % Evaluate M_1(t) at time points T
    mt_minf_1 = double(M_1(T));
    
    % Plotting results
    figure();
    plot(T / 3600, 100 * mt_minf_1, 'r', 'LineWidth', 1);
    hold on;
    p = plot(Time, 100 * Diffusion_Fraction, 'ko ');
    p.MarkerFaceColor = [0 0 0];
    p.MarkerSize = 8;
    hold off;
    
    % Formatting plot
    xlabel('Time (hours)', 'Interpreter', 'latex');
    ylabel('$$M_{t}/M_{\infty} (\%)$$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 16);
    legend('Fitted Data', 'Location', 'best', 'Interpreter', 'latex');
    title(['Diffusion Fraction vs. Time ($$D = ' num2str(round(D_opt / 1e-6, 2)) ' \times 10^{-6} cm^{2}/s$$)'], 'Interpreter', 'latex');
    box on;
    
end

function cost = costfun(D)
    global Time Diffusion_Fraction
    
    % Parameters
    L = 0.04;
    syms n t
    
    % Define function f_1(t)
    f_1(t) = (8 / ((2*n+1)^2 * pi^2)) * exp((-D * (2*n+1)^2 * pi^2 * t) / (4 * L^2));
    
    % Define M_1(t) using symsum
    M_1(t) = 1 - symsum(f_1(t), n, 0, Inf);
    
    % Evaluate M_1(t) at Time points
    T = 3600 * Time;
    mt_minf_1 = double(M_1(T));
    
    % Calculate cost (objective function)
    output = mt_minf_1';
    cost = 1e4 * sum((Diffusion_Fraction' - output).^2);
end
