% MARK_FEATURES: mark features for geometric refinement according to 
%  the marking strategy in adaptivity_data.mark_strategy, taking into account 
%  the computed error estimators in the variable estimators.
%
%   [marked, nmarked] = mark_features (estimators, adaptivity_data)
%
% INPUT:
%
%   estimators:      result obtained from the estimate step (see estimate_defeaturing_error_H1s)
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it must contain the following fields:
%      - mark_strategy: the possible marking strategies are:
%           GERS: guaranteed error reduction strategy (Dörfler's),
%           MS:   maximum strategy,
%           GR:   global (uniform) refinement.
%      - mark_param:    parameter for marking, 0 < mark_param < 1.
%
% OUTPUT:
%
%    marked:  indices of the marked features
%    nmarked: total number of marked features
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
%               2024 Ondine Chanon
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [marked, nmarked] = mark_features(estimators, adaptivity_data)
    aux_marked = zeros(size (estimators));
    
    switch adaptivity_data.mark_strategy
        case 'GR'
            aux_marked = ones(size (estimators));
        case 'MS'
            max_estimators = max(estimators);
            aux_marked(estimators > adaptivity_data.mark_param * max_estimators) = 1;
        case 'GERS'
            est_sum2 = sum(estimators.^2);
            [est2_ordered, perm] = sort(estimators.^2, 'descend');
            index = find(cumsum (est2_ordered) > (1 - adaptivity_data.mark_param)^2 * est_sum2, 1, 'first');
            index = find(est2_ordered > 0.999 * est2_ordered(index), 1, 'last');
            aux_marked(perm(1:index)) = 1;
    end
    
    marked = find(aux_marked);
    nmarked = numel(marked);
end
