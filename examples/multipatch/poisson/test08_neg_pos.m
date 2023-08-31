clear
close all
clc

% -- Test reproducing [AC2023+], Section 6.1.2, geometry \Omega_\delta, \delta > 0 --
% -- For delta = 0, see test07_complex.m ; for delta < 0, see test06_complex.m --

%% Set problem data
% -- Choose below the desired delta > 0 to test --
% delta = 2e-1;
delta = 2e-4; 
% --

feature_side_length = 0.1;

filename = 'results/test08_neg_pos';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 40 * cos(pi * x) + 10 * cos(5 * pi * x);

[problem_data, problem_data_0, problem_data_F] = determineBC(problem_data);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32]; 
method_data.nquad = [5 5];


%% Main

% 1) BUILD GEOMETRY
[srf_0, srf, srf_Fpos] = buildGeometry(feature_side_length, delta); 
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
[error_h1s, error_h1s_from_boundary, error_h1s_F] = errh1s_positive(msh, space, u, msh_0, ...
    space_0, u_0, msh_F, space_F, u_0tilde, problem_data.omega0_patches, ...
    problem_data_F.F_patches, problem_data_0.omega_patches);

% 4b) COMPUTE ESTIMATOR AND ERROR FROM BOUNDARY
[estimator_positive_part, measure_of_gamma0, ~, error_h1s_from_positive_boundary] = ...
    est_positive(msh_F, space_F, u_0tilde, problem_data_F.gamma0_sides, problem_data_F.gammae_sides,...
        problem_data_0.g, problem_data.g, problem_data_F.F_patches, ...
        problem_data.omega0_patches, space, u);
[estimator_negative_part, measure_of_gamma, error_h1s_from_negative_boundary] = ...
    est_negative(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
        problem_data_0.omega_patches, problem_data.gamma_sides, ...
        problem_data.omega0_patches, msh, space, u);
        
estimator = sqrt(estimator_positive_part^2 + estimator_negative_part^2); 
errh1s_interface = sqrt(error_h1s_from_positive_boundary^2 + error_h1s_from_negative_boundary^2);


%% Display and save the results
if saveIt
    save(filename, 'feature_side_length', 'delta', 'error_h1s', 'error_h1s_from_boundary', ...
        'error_h1s_F', 'errh1s_interface', 'measure_of_gamma0', 'measure_of_gamma', ...
        'error_h1s_from_positive_boundary', 'error_h1s_from_negative_boundary', ...
        'estimator', 'estimator_positive_part', 'estimator_negative_part')
end

fprintf('For delta = %e, \n', delta)
fprintf('    * Estimator                           = %e \n', estimator)
fprintf('    * Defeaturing error |u-u_0|_{1,Omega} = %e \n', error_h1s)
fprintf('    * Effectivity index                   = %f \n', estimator/error_h1s)


%% Auxiliary functions
function [srf_0, srf, srf_Fpos] = buildGeometry(epsilon, delta)
    extension_factor = 4; 

    srf_tot(1) = nrb4surf([0, 0], [0.5 - epsilon - delta / 2, 0], ...
                          [0, 1 - epsilon], [0.5 - epsilon - delta / 2, 1 - epsilon]);
    srf_tot(1) = nrbkntins(srf_tot(1), {(1 - extension_factor * epsilon) / (1 - epsilon), ...
                                         1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(2) = nrb4surf([0.5 - epsilon - delta / 2, 0], [0.5 - delta / 2, 0], ...
                          [0.5 - epsilon - delta / 2, 1 - epsilon], [0.5 - delta / 2, 1 - epsilon]);
    srf_tot(2) = nrbkntins(srf_tot(2), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(3) = nrb4surf([0.5 - delta / 2, 0], [0.5 + delta / 2, 0], ...
                          [0.5 - delta / 2, 1 - epsilon], [0.5 + delta / 2, 1 - epsilon]);
    srf_tot(3) = nrbkntins(srf_tot(3), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(4) = nrb4surf([0.5 + delta / 2, 0], [0.5 + epsilon + delta / 2, 0], ...
                          [0.5 + delta / 2, 1 - epsilon], [0.5 + epsilon + delta / 2, 1 - epsilon]);
    srf_tot(4) = nrbkntins(srf_tot(4), {min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(5) = nrb4surf([0.5 + epsilon + delta / 2, 0], [1, 0], ...
                          [0.5 + epsilon + delta / 2, 1 - epsilon], [1, 1 - epsilon]);
    srf_tot(5) = nrbkntins(srf_tot(5), {1 - (1 - extension_factor * epsilon) / (1 - epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    
    srf_tot(6) = nrb4surf([0, 1 - epsilon], [0.5 - epsilon - delta / 2, 1 - epsilon], ...
                          [0, 1], [0.5 - epsilon - delta / 2, 1]);
    srf_tot(6) = nrbkntins(srf_tot(6), {(1 - extension_factor * epsilon) / (1 - epsilon), ...
                                         1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(7) = nrb4surf([0.5 - epsilon - delta / 2, 1 - epsilon], [0.5 - delta / 2, 1 - epsilon], ...
                          [0.5 - epsilon - delta / 2, 1], [0.5 - delta / 2, 1]);
    srf_tot(7) = nrbkntins(srf_tot(7), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(8) = nrb4surf([0.5 - delta / 2, 1 - epsilon], [0.5 + delta / 2, 1 - epsilon], ...
                          [0.5 - delta / 2, 1], [0.5 + delta / 2, 1]);
    srf_tot(8) = nrbkntins(srf_tot(8), {0.5, 1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(9) = nrb4surf([0.5 + delta / 2, 1 - epsilon], [0.5 + epsilon + delta / 2, 1 - epsilon], ...
                          [0.5 + delta / 2, 1], [0.5 + epsilon + delta / 2, 1]);
    srf_tot(9) = nrbkntins(srf_tot(9), {min(0.5, (delta / 2 * extension_factor) / epsilon), ...
                                        1 - (extension_factor - 1) * epsilon / 2});
    
    srf_tot(10) = nrb4surf([0.5 + epsilon + delta / 2, 1 - epsilon], ...
                           [1, 1 - epsilon], [0.5 + epsilon + delta / 2, 1], [1, 1]);
    srf_tot(10) = nrbkntins(srf_tot(10), {1-(1-extension_factor*epsilon)/(1-epsilon), ...
                                          1-(extension_factor-1)*epsilon/2});
    
    srf_tot(11) = nrb4surf([0.5 - epsilon - delta / 2, 1], [0.5 - delta / 2, 1], ...
                           [0.5 - epsilon - delta / 2, 1 + epsilon], [0.5 - delta / 2, 1 + epsilon]);
    srf_tot(11) = nrbkntins(srf_tot(11), {1 - min(0.5, (delta / 2 * extension_factor) / epsilon), 0.5});
    

    srf_0 = srf_tot(1:10);
    srf_Fpos = srf_tot(11);
    srf = srf_tot([1:8 10:11]);
end

function [problem_data, problem_data_0, problem_data_F] = determineBC(problem_data)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F = rmfield(problem_data_F, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1 6 7 9:18]; 
    problem_data.drchlt_sides = [2 3 4 5 8];
    problem_data.gamma_sides = [6 11 13]; % for negative features
    problem_data.omega0_patches = 1:9;

    % Simplified problem
    problem_data_0.nmnn_sides = [1 6 8:14];
    problem_data_0.drchlt_sides = [2 3 4 5 7];
    problem_data_0.gamma0_sides = 10; % for the Dirichlet extension
    problem_data_0.omega_patches = [1:8 10];
    problem_data_0.gamma_sides = cell(10, 1); % relative to each patch
    problem_data_0.gamma_sides([4 8 10]) = {4, 2, 1};
    
    % Extension problem
    problem_data_F.nmnn_sides = [1 2 4]; 
    problem_data_F.drchlt_sides = 3; 
    problem_data_F.gamma0_sides = 3; 
    problem_data_F.gammae_sides = {[]}; % relative to each patch
    problem_data_F.F_patches = 1; 
end
