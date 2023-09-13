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
epsilon_values = 0.1;
shifts = -delta; 
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 40 * cos(pi * x) + 10 * cos(5 * pi * x);
% --

filename = 'results/test06_complex';
saveIt = false;
plotIt = true;

[problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data);

method_data.degree = [4 4];
method_data.regularity = [3 3];
method_data.nsub = [16 16]; 
method_data.nquad = [6 6];


%% Main
number_of_epsilons = numel(epsilon_values);
error_H1s = zeros(1, number_of_epsilons);
error_H1s_Omega_star = zeros(1, number_of_epsilons);
error_H1s_Fp = zeros(1, number_of_epsilons);
error_H1s_boundary_representation = zeros(1, number_of_epsilons);
error_H1s_boundary_representation_Fp = zeros(1, number_of_epsilons);
error_H1s_boundary_representation_Fn = zeros(1, number_of_epsilons);
estimator = zeros(1, number_of_epsilons);
estimator_Fp = zeros(1, number_of_epsilons);
estimator_Fn = zeros(1, number_of_epsilons);
norm_of_u = zeros(1, number_of_epsilons);
relative_error_H1s = zeros(1, number_of_epsilons);

for iter = 1:number_of_epsilons
    epsilon = epsilon_values(iter);
    shift = shifts(iter);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_Fp] = build_geometry(epsilon, shift); 
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_Fp.geo_name = srf_Fp;

    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_generalized(problem_data, method_data);

    % 3a) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_generalized(problem_data_0, method_data);
    problem_data_Fp.h = extract_boundary_dofs(msh_0, space_0, u_0, problem_data_0.gamma0_sides);

    % 3b) SOLVE THE EXTENSION PROBLEM
    [Fp, msh_Fp, space_Fp, u_0tilde] = mp_solve_laplace_generalized(problem_data_Fp, method_data);

    % 4a) COMPUTE THE DEFEATURING ESTIMATOR
    [estimator_Fp(iter), ~, ~, error_H1s_boundary_representation_Fp(iter)] = ...
        est_positive(msh_Fp, space_Fp, u_0tilde, problem_data_Fp.gamma0_sides, problem_data_Fp.gammae_sides,...
            problem_data_0.g, problem_data.g, problem_data_Fp.F_patches, ...
            problem_data.omega0_patches, space, u);
    [estimator_Fn(iter), ~, error_H1s_boundary_representation_Fn(iter)] = ...
        est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);

    estimator(iter) = sqrt(estimator_Fp(iter)^2 + estimator_Fn(iter)^2); 

    error_H1s_boundary_representation(iter) = sqrt(error_H1s_boundary_representation_Fp(iter)^2 ...
                                                   + error_H1s_boundary_representation_Fn(iter)^2);
    
    % 4b) COMPUTE THE DEFEATURING ERROR
    [error_H1s(iter), error_H1s_Omega_star(iter), error_H1s_Fp(iter)] = ...
        defeaturing_error_H1s(msh_0, space_0, u_0, problem_data_0.omega_patches, ...
                              msh, space, u, problem_data.omega0_patches, ...
                              msh_Fp, space_Fp, u_0tilde, problem_data_Fp.F_patches);

    norm_of_u(iter) = error_H1s_in_patches(msh, space, u, 1:msh.npatch, ...
                                           msh, space, zeros(size(u)));
    relative_error_H1s(iter) = error_H1s(iter) / norm_of_u(iter);
end


%% Display and save the results
if saveIt
    save(filename, 'epsilon_values', 'error_H1s', 'error_H1s_Omega_star', 'error_H1s_Fp', ...
        'error_H1s_boundary_representation', 'error_H1s_boundary_representation_Fp', ...
        'error_H1s_boundary_representation_Fn', 'estimator', 'estimator_Fp', ...
        'estimator_Fn', 'norm_of_u', 'relative_error_H1s')
end
if plotIt
    fig = figure;
    loglog(epsilon_values, error_H1s, '+-r', epsilon_values, estimator, '+-b', ...
           epsilon_values, epsilon_values*400, 'k:');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon', 'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

fprintf("For epsilon (feature's side length) = %e, \n", epsilon)
fprintf('    * Estimator                           = %e \n', estimator(end))
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_H1s(end))
fprintf('    * Effectivity index                   = %f \n', estimator(end) / error_H1s(end))


%% Auxiliary functions
function [srf_0, srf, srf_Fpos] = build_geometry(epsilon, shift)
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

function [problem_data, problem_data_0, problem_data_Fp] = set_boundary_conditions(problem_data)
    problem_data_0 = problem_data;
    problem_data_Fp = problem_data;
    problem_data_Fp = rmfield(problem_data_Fp, 'h');
    
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
    problem_data_Fp.nmnn_sides = [1 3:6]; 
    problem_data_Fp.drchlt_sides = 2; 
    problem_data_Fp.gamma0_sides = 2; 
    problem_data_Fp.gammae_sides = cell(2, 1); % relative to each patch
    problem_data_Fp.F_patches = 1:2; 
end
