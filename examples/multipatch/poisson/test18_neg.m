clear
close all
clc

eps_values = 1e-2 ./ 2.^(0:6); % sqrt(1e-2^2/pi); % sqrt(0.1^2/pi); % 2*0.1/pi;
filename = 'results/test18_drchlt0';
saveIt = false;
plotIt = true;

problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) ones(size(x)); 
problem_data.h = @(x, y, ind) zeros(size(x));
[problem_data, problem_data_0] = determineBC(problem_data);
ext = 4;

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [32 32];
method_data.nquad = [5 5];

%%
nn = numel(eps_values);
errh1s = zeros(1,nn);
errh1s_interface = zeros(1,nn);
est = zeros(1,nn);
meas_gamma = zeros(1,nn);
normu = zeros(1,nn);
errh1s_rel = zeros(1,nn);

for ii = 1:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace (problem_data, method_data);
    
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace (problem_data_0, method_data);
    
    % 4) COMPUTE ERROR AND ESTIMATOR
    errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    normu(ii) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    errh1s_rel(ii) = errh1s(ii)/normu(ii);
%     [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
%         est_negative(space_0, u_0, problem_data.gamma_sides, problem_data_0.omega_patches, problem_data.g, msh, space, u);

    [est(ii), meas_gamma(ii), errh1s_interface(ii)] = ...
        est_negative_v2(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
            problem_data_0.omega_patches, problem_data.gamma_sides, ...
            problem_data.omega0_patches, msh, space, u);
end

%%
if saveIt
    save(filename, 'eps_values', 'errh1s', 'errh1s_interface', 'est', 'meas_gamma', 'normu', 'errh1s_rel')
end
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r', eps_values, est, '+-b', eps_values, eps_values.^2.*abs(log(eps_values))*0.8, 'k:', eps_values, eps_values.^2*5, 'k.-');
    grid on
    legend('|u-u_0|_{1,\Omega}', 'Estimator', '\epsilon^2 log(\epsilon)', '\epsilon^2',  'Location', 'northwest')
    if saveIt
        saveas(fig, filename, 'epsc');
    end
end

%%
function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    L = sqrt(2)*epsilon/4;
    srf_0(1) = nrbruled (nrbcirc(epsilon, [0 0], pi/4, 3*pi/4), nrbcirc(1, [0 0], pi/4, 3*pi/4));
    srf_0(2) = nrbtform(srf_0(1), vecrotz(pi/2));
    srf_0(3) = nrbtform(srf_0(2), vecrotz(pi/2));
    srf_0(4) = nrbtform(srf_0(3), vecrotz(pi/2));
    for ii = 1:4
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (ext-1)*epsilon/(1-epsilon)});
    end
    srf_0(5) = nrbruled (nrbline([L L], [-L L]), nrbcirc(epsilon, [0 0], pi/4, 3*pi/4));
    srf_0(6) = nrbtform(srf_0(5), vecrotz(pi/2));
    srf_0(7) = nrbtform(srf_0(6), vecrotz(pi/2));
    srf_0(8) = nrbtform(srf_0(7), vecrotz(pi/2));
    for ii = 5:8
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L -L], [L, -L], [-L L], [L L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:4);
    srf_F = srf_0(5:9);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [1 3 5 7];
    problem_data.drchlt_sides = [2 4 6 8];
    problem_data.gamma_sides = [1 3 5 7];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
    problem_data_0.omega_patches = 1:4;
    problem_data_0.gamma_sides = cell(9,1); 
    problem_data_0.gamma_sides(1:4) = {3, 3, 3, 3}; % relative to each patch
end