clear
close all
clc

% -- Test reproducing [BCV2022], Section 6.3.2 --

%% Set problem data
% -- Choose below the desired feature extension, following the notation in [BCV2022] --
% feature_extension = 'F';
% feature_extension = 'Ftilde1';
feature_extension = 'Ftilde2';
% --

filename = 'results/test03_pos';
saveIt = false;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) cos(pi * x) + 10 * cos(5 * pi * x);

[problem_data, problem_data_0, problem_data_F] = determineBC(problem_data, feature_extension);

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [2 2].^5; 
method_data.nquad = [5 5];


%% Main
% 1) BUILD GEOMETRY
[srf_0, srf, srf_F, srf_Ftilde] = buildGeometry(feature_extension);
problem_data.geo_name = srf;
problem_data_0.geo_name = srf_0;
if strcmp(feature_extension, 'F')
    problem_data_F.geo_name = srf_F;
else
    problem_data_F.geo_name = srf_Ftilde;
end

% 2) SOLVE THE EXACT PROBLEM
[omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);

% 3a) SOLVE THE DEFEATURED PROBLEM
[omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
problem_data_F.h = buildDirichletConditionExtPb(msh_0, problem_data_0, space_0, u_0);

% 3b) SOLVE THE (EXTENDED) EXTENSION PROBLEM
[Ftilde, msh_Ftilde, space_Ftilde, u_0tilde] = my_mp_solve_laplace (problem_data_F, method_data);

% 4) COMPUTE ERROR AND ESTIMATOR
[error_h1s, error_h1s_0, error_h1s_F] = errh1s_positive(msh, space, u, ...
    msh_0, space_0, u_0, msh_Ftilde, space_Ftilde, u_0tilde, ...
    problem_data.omega0_patches, problem_data_F.F_patches);

[estimator, measure_of_gamma0, measure_of_gammae, error_h1s_from_boundary] = ...
    est_positive(msh_Ftilde, space_Ftilde, u_0tilde, problem_data_F.gamma0_sides,...
        problem_data_F.gammae_sides, problem_data_0.g, problem_data.g,...
        problem_data_F.F_patches, problem_data.omega0_patches, space, u);


%% Display and save the results
if saveIt
    save(filename, 'feature_extension', 'error_h1s', 'error_h1s_0', 'error_h1s_F', ...
        'error_h1s_from_boundary', 'estimator', 'measure_of_gamma0', 'measure_of_gammae')
end

fprintf('For feature extension = %s, \n', feature_extension)
fprintf('    * Estimator                                          = %e \n', estimator)
fprintf('    * Defeaturing error |u-u_0|_{1,Omega}                = %e \n', error_h1s)
fprintf('    * Defeaturing error in Omega_0, |u-u_0|_{1,Omega_0}  = %e \n', error_h1s_0)
fprintf('    * Defeaturing error in F, |u-u_0|_{1,F}              = %e \n', error_h1s_F)
fprintf('    * Effectivity index                                  = %f \n', estimator/error_h1s)


%% Auxiliary functions
function [problem_data, problem_data_0, problem_data_F] = determineBC(problem_data, feature_extension)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F = rmfield(problem_data_F, 'h');
    
    % Exact problem
    problem_data.nmnn_sides = [1:3 5 7 8]; 
    problem_data.drchlt_sides = [4 6];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [1:5 7 8 10];
    problem_data_0.drchlt_sides = [6 9];
    problem_data_0.gamma0_sides = [2 5 7 10]; 
    
    if strcmp(feature_extension, 'F')
        % Extension problem
        problem_data_F.nmnn_sides = [2 5]; 
        problem_data_F.drchlt_sides = [1 3 4 6];
        problem_data_F.gamma0_sides = [1 3 4 6];
        problem_data_F.F_patches = 1:3;
        problem_data_F.gammae_sides = cell(3,1);
    else
        % Extended extension problem
        problem_data_F.nmnn_sides = 5:8; 
        problem_data_F.drchlt_sides = 1:4;
        problem_data_F.gamma0_sides = 1:4;
        problem_data_F.F_patches = 1:3;
        problem_data_F.gammae_sides = {4, [], 2}; % side local to patch, for each patch
    end
end

function [srf_0, srf, srf_F, srf_Ftilde] = buildGeometry(feature_extension)
    epsilon = 0.5;
    srf(6) = nrbsquare([0.5, 0.5], epsilon / 5, epsilon / 5);
    srf(6) = nrbdegelev(srf(6), [1, 1]);
    srf(6) = nrbkntins(srf(6), {0.5, 0.5});

    circle = nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], pi, 5 * pi / 4);
    line = nrbline([0.5, 0.5 + epsilon / 5], [0.5 + epsilon / 5, 0.5 + epsilon / 5]);
    srf(5) = nrbruled(line, circle);
    srf(5) = nrbdegelev(srf(5), [0, 1]);
    srf(5) = nrbkntins(srf(5), {0.5, 0.5});

    circle = nrbreverse(nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], 5 * pi / 4, 3 * pi / 2));
    line = nrbline([0.5 + epsilon / 5, 0.5], [0.5 + epsilon / 5, 0.5 + epsilon / 5]);
    srf(7) = nrbruled(line, circle);
    srf(7) = nrbdegelev(nrbtransp(srf(7)), [1, 0]);
    srf(7) = nrbkntins(srf(7), {0.5, 0.5});

    srf(1) = nrbsquare([0, 0.5 + epsilon / 5], 0.5, 4 * epsilon / 5);

    line1 = nrbline([0, 0], [0.5, 0.5]);
    line2 = nrbline([0, 0.5 + epsilon / 5], [0.5, 0.5 + epsilon / 5]);
    srf(2) = nrbruled(line1, line2);

    line1 = nrbline([0, 0], [0.5 + epsilon / 5, 0]);
    line2 = nrbline([0.5, 0.5], [0.5 + epsilon / 5, 0.5]);
    srf(3) = nrbruled(line1, line2);

    srf(4) = nrbsquare([0.5 + epsilon / 5, 0], 4 * epsilon / 5, 0.5);

    for ii = 1:4
        srf(ii) = nrbdegelev(srf(ii), [1, 1]);
        srf(ii) = nrbkntins(srf(ii), {0.5, 0.5});
    end

    if strcmp(feature_extension, 'Ftilde1')
        circle = nrbreverse(nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], pi, 5 * pi / 4));
        line = nrbline([1 - 3 * epsilon / 5, 1 - 3 * epsilon / 5], [1 - 3 * epsilon / 5, 1]);
        srf(8) = nrbtransp(nrbruled(circle, line));
        srf(8) = nrbkntins(nrbdegelev(srf(8), [1, 0]), {0.5, 0.5});

        circle = nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], 5 * pi / 4, 3 * pi / 2);
        line = nrbline([1 - 3 * epsilon / 5, 1 - 3 * epsilon / 5], [1, 1 - 3 * epsilon / 5]);
        srf(9) = nrbtransp(nrbruled(circle, line));
        srf(9) = nrbkntins(nrbdegelev(srf(9), [1, 0]), {0.5, 0.5});

        srf(10) = nrbsquare([1 - 3 * epsilon / 5, 1 - 3 * epsilon / 5], ...
                            3 * epsilon / 5, 3 * epsilon / 5);
        srf(10) = nrbkntins(nrbdegelev(srf(10), [1,1]), {0.5, 0.5});

    elseif strcmp(feature_extension, 'Ftilde2')
        internal_radius = 0.25; 
        cir1 = nrbreverse(nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], pi, 5*pi/4));
        cir2 = nrbreverse(nrbcirc(internal_radius, [0.5 + epsilon, 0.5 + epsilon], pi, 5 * pi / 4));
        srf(8) = nrbtransp(nrbruled(cir1, cir2));
        srf(8) = nrbkntins(srf(8), {0.5, 0.5});

        cir1 = nrbreverse(nrbcirc(epsilon, [0.5 + epsilon, 0.5 + epsilon], 5 * pi / 4, 3 * pi / 2));
        cir2 = nrbreverse(nrbcirc(internal_radius, [0.5 + epsilon, 0.5 + epsilon], ...
                                  5 * pi / 4, 3 * pi / 2));
        srf(9) = nrbtransp(nrbruled(cir1, cir2));
        srf(9) = nrbkntins(srf(9), {0.5, 0.5});
    end

    srf_0 = srf(1:4);
    srf_F = srf(5:7);
    if strcmp(feature_extension, 'Ftilde1')
        srf_Ftilde = srf(5:10);
    elseif strcmp(feature_extension, 'Ftilde2')
        srf_Ftilde = srf(5:9);
    elseif strcmp(feature_extension, 'F')
        srf_Ftilde = srf_F;
    else
        error('Unknown feature extension.')
    end
    srf = srf(1:7);
end
