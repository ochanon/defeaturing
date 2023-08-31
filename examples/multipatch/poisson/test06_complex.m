clear
close all
clc

% WARNING: This simulation may be long and may require a large amount of memory.
% Indeed, it requires a high number of degrees of freedom (in particular
% around and in the feature) in order to correcly capture the behavior of 
% the defeaturing error. Otherwise, the numerical approximation error
% influences the defeaturing error and cannot be neglected. 

%% Set problem data
% -- To reproduce [BCV2022], Section 6.2.1, geometry \Omega_\varepsilon^6 --
% eps_values = 1e-2 ./ 2.^(0:6);
% shifts = eps_values ./ 2;
% problem_data.c_diff = @(x, y) ones(size(x));
% problem_data.g = @(x, y, ind) zeros(size(x));
% problem_data.f = @(x, y) 10 * cos(3 * pi * x) .* sin(5 * pi * y); 
% problem_data.h = @(x, y, ind) zeros(size(x));
% --

% -- To reproduce [AC2023+], Section 6.1.2, geometry \Omega_\delta, \delta < 0 --
% -- For delta > 0, see test08_neg_pos.m ; for delta = 0, see test07_complex.m --
% -- Uncomment below the value of delta you want to test --
delta = -1e-3;
% delta = -9.9e-2;
% -- Uncomment the following lines of code --
eps_values = 0.1;
shifts = -delta; 
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 40 * cos(pi * x) + 10 * cos(5 * pi * x);
% --

filename = 'results/test06_complex';
saveIt = false;
plotIt = true;

[problem_data, problem_data_0, problem_data_F] = determineBC(problem_data);

method_data.degree = [4 4];
method_data.regularity = [3 3];
method_data.nsub = [16 16]; 
method_data.nquad = [6 6];


%% Main
number_of_epsilons = numel(eps_values);
error_h1s = zeros(1, number_of_epsilons);
error_h1s_Omega_star = zeros(1, number_of_epsilons);
error_h1s_F = zeros(1, number_of_epsilons);
error_h1s_from_boundary = zeros(1, number_of_epsilons);
error_h1s_from_positive_boundary = zeros(1, number_of_epsilons);
error_h1s_from_negative_boundary = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
estimator_positive_feature = zeros(1, number_of_epsilons);
estimator_negative_feature = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_h1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = eps_values(iter);
    shift = shifts(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_Fpos] = buildGeometry(epsilon, shift); 
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_F.geo_name = srf_Fpos;

    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);

    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    problem_data_F.h = buildDirichletConditionExtPb(msh_0, problem_data_0, space_0, u_0);

    % 3b) SOLVE THE EXTENSION PROBLEM
    [F, msh_F, space_F, u_0tilde] = my_mp_solve_laplace (problem_data_F, method_data);

    % 4a) COMPUTE ERROR
    [error_h1s(iter), error_h1s_Omega_star(iter), error_h1s_F(iter)] = errh1s_positive(msh, space, u, msh_0, ...
        space_0, u_0, msh_F, space_F, u_0tilde, problem_data.omega0_patches, ...
        problem_data_F.F_patches, problem_data_0.omega_patches);

    % 4b) COMPUTE ESTIMATOR AND ERROR FROM BOUNDARY
    [estimator_positive_feature(iter), ~, ~, error_h1s_from_positive_boundary(iter)] = ...
        est_positive(msh_F, space_F, u_0tilde, problem_data_F.gamma0_sides, problem_data_F.gammae_sides,...
            problem_data_0.g, problem_data.g, problem_data_F.F_patches, ...
            problem_data.omega0_patches, space, u);
    [estimator_negative_feature(iter), ~, error_h1s_from_negative_boundary(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);

    estimator(iter) = sqrt(estimator_positive_feature(iter)^2 + estimator_negative_feature(iter)^2); 
    error_h1s_from_boundary(iter) = sqrt(error_h1s_from_positive_boundary(iter)^2 + error_h1s_from_negative_boundary(iter)^2);
    
    norm_of_u(iter) = errh1s_positive(msh, space, u, ...
        msh_0, space_0, zeros(size(u_0)), msh_F, space_F, zeros(size(u_0tilde)), problem_data.omega0_patches, ...
        problem_data_F.F_patches, problem_data_0.omega_patches);
    relative_error_h1s(iter) = error_h1s(iter)/norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'eps_values', 'error_h1s', 'error_h1s_Omega_star', 'error_h1s_F', ...
        'error_h1s_from_boundary', 'error_h1s_from_positive_boundary', 'error_h1s_from_negative_boundary',...
        'estimator', 'estimator_positive_feature', 'estimator_negative_feature', 'norm_of_u', 'relative_error_h1s')
end
if plotIt
    fig = figure;
    loglog(eps_values, error_h1s, '+-r', eps_values, estimator, '+-b', ...
           eps_values, eps_values*400, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf("For epsilon (feature's side length) = %e, \n", epsilon)
fprintf('    * Estimator                           = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_h1s(end))
fprintf('    * Effectivity index                   = %f \n', estimator(end)/error_h1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_Fpos] = buildGeometry(epsilon, shift)
    srf_tot(1) = nrbsquare([0, 0], 0.5 - epsilon + shift / 2, 1 - epsilon);
    srf_tot(2) = nrbsquare([0.5 - epsilon + shift / 2, 0], epsilon - shift, 1 - epsilon);
    srf_tot(3) = nrbsquare([0.5 - shift / 2, 0], shift, 1 - epsilon);
    srf_tot(4) = nrbsquare([0.5 + shift / 2, 0], epsilon - shift, 1 - epsilon);
    srf_tot(5) = nrbsquare([0.5 + epsilon - shift / 2, 0], 0.5 - epsilon + shift / 2, 1 - epsilon);
    
    srf_tot(6) = nrbsquare([0, 1 - epsilon], 0.5 - epsilon + shift / 2, epsilon);
    srf_tot(7) = nrbsquare([0.5 - epsilon + shift / 2, 1 - epsilon], epsilon - shift, epsilon);
    srf_tot(8) = nrbsquare([0.5 - shift / 2, 1 - epsilon], shift, epsilon);
    srf_tot(9) = nrbsquare([0.5 + shift / 2, 1 - epsilon], epsilon - shift, epsilon);
    srf_tot(10) = nrbsquare([0.5 + epsilon - shift / 2, 1 - epsilon], 0.5 - epsilon + shift / 2, epsilon);
    
    srf_tot(11) = nrbsquare([0.5 - epsilon + shift / 2, 1], epsilon - shift, epsilon);
    srf_tot(12) = nrbsquare([0.5 - shift / 2, 1], shift, epsilon);
    
    knewknots = 8;
    srf_tot(1) = nrbkntins(srf_tot(1), {1 - epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (0.5 - epsilon + shift / 2), ...
                                        1 - epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (1 - epsilon)}); 
    for iptc = 2:4
        srf_tot(iptc) = nrbkntins(srf_tot(iptc), {[2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)], ...
                                                   1 - epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (1 - epsilon)});
    end
    srf_tot(5) = nrbkntins(srf_tot(5), {epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (0.5 - epsilon + shift / 2), ...
                                        1 - epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (1 - epsilon)}); 
    srf_tot(6) = nrbkntins(srf_tot(6), {1 - epsilon .* 2.^(-(knewknots - 1) *2:0) ./ (0.5 - epsilon + shift / 2), ...
                                        [2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)]}); 
    for iptc = 7:9
        srf_tot(iptc) = nrbkntins(srf_tot(iptc), {[2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)], ...
                                                  [2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)]});
    end
    srf_tot(10) = nrbkntins(srf_tot(10), {epsilon .* 2.^(-(knewknots - 1) * 2:0) ./ (0.5 - epsilon + shift / 2), ...
                                          [2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)]});  
    srf_tot(11) = nrbkntins(srf_tot(11), {[2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)], 2.^(-(knewknots) * 2 + 1:-1)});
    srf_tot(12) = nrbkntins(srf_tot(12), {[2.^(-knewknots:-1), 1 - 2.^-(2:knewknots)], 2.^(-(knewknots) * 2 + 1:-1)});
    
    srf_0 = srf_tot(1:10);
    srf_Fpos = srf_tot([11 12]);
    srf = srf_tot([1:7 10:12]);
end

function [problem_data, problem_data_0, problem_data_F] = determineBC(problem_data)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F = rmfield(problem_data_F, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [5 7 11:13 15:20]; 
    problem_data.drchlt_sides = [1:4 6 8:10 14];
    problem_data.gamma_sides = [5 7 12 13]; % for the negative feature component
    problem_data.omega0_patches = 1:8;
    
    % Simplified problem
    problem_data_0.nmnn_sides = [9:12 14];
    problem_data_0.drchlt_sides = [1:8 13];
    problem_data_0.gamma0_sides = 10; % for the Dirichlet extension
    problem_data_0.omega_patches = [1:7 10];
    problem_data_0.gamma_sides = cell(10, 1); % relative to each patch
    problem_data_0.gamma_sides([3 4 7 10]) = {4, 4, 2, 1};
    
    % Extension problem
    problem_data_F.nmnn_sides = [1 3:6]; 
    problem_data_F.drchlt_sides = 2; 
    problem_data_F.gamma0_sides = 2; 
    problem_data_F.gammae_sides = cell(2, 1); % relative to each patch
    problem_data_F.F_patches = 1:2; 
end
