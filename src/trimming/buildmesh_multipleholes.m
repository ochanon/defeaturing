function geo_name = buildmesh_multipleholes(nel, centers, radiuses, degree)

nel_in_extension = 2 * nel;
nel_in_feature = 2 * nel;
extension_factor = 4; 

if length(nel) == 1
    nel = nel*ones(1, size(centers,2));
end
if length(nel_in_extension) == 1
    nel_in_extension = nel_in_extension*ones(1, size(centers,2));
end
if length(nel_in_feature) == 1
    nel_in_feature = nel_in_feature*ones(1, size(centers,2));
end

knots = cell(1, size(centers, 2));

for idim = 1:size(centers, 2)
    [centers_x, Isort] = unique(centers(:, idim)); % centers sorted
    radiuses_x = radiuses(Isort);
    h_out = 1/nel(idim);

    h_ext = (extension_factor-1)*radiuses_x./nel_in_extension(idim); 
    h_bd = radiuses_x./nel_in_feature(idim); 
    cpr = [centers_x-radiuses_x; centers_x+radiuses_x];
    hh = [h_bd/2; h_bd/2];
    hh_ext = [h_ext; h_ext];
    indicators = [ones(length(centers_x), 1); -ones(length(centers_x), 1)];
    
    [cpr, sortid] = sort(cpr);
    indicators = indicators(sortid);
    hh = hh(sortid);
    hh_ext = hh_ext(sortid);
    
    id1 = find(cpr>=0, 1, 'first');
    idend = find(cpr<=1, 1, 'last');
    cpr = cpr(id1:idend);
    indicators = indicators(id1:idend);
    hh = hh(id1:idend);
    hh_ext = hh_ext(id1:idend);
    
    indicators = cumsum(indicators);
    left = id1 - 1;
    
    % get rid of multiple centers
    [cpr, idunique] = unique(cpr);
    hhunique = zeros(length(idunique), 1);
    hh_extunique = zeros(length(idunique), 1);
    indicatorsunique = zeros(length(idunique), 1);
    idunique = [idunique; length(hh)];
    for ii = 1:length(idunique)-1
        if ii ~= length(idunique)-1
            ids = idunique(ii):idunique(ii+1)-1;
        else 
            ids = idunique(ii):idunique(ii+1);
        end
        hhunique(ii) = max(hh(ids));
        hh_extunique(ii) = max(hh_ext(ids));
        indicatorsunique(ii) = indicators(ids(end));
    end
    indicators = indicatorsunique;
    hh = hhunique;
    hh_ext = hh_extunique;
    
    % take care of boundaries
    if cpr(1)~=0
        indicators = [left; indicators];
        cpr = [0; cpr];
        hh = [0; hh];
        hh_ext = [0; hh_ext];
        flag0 = true;
    else
        flag0 = false;
    end
    if indicators(end)==0 && cpr(end)~=1
        indicators = [indicators; 0];
        cpr = [cpr; 1];
        hh = [hh; 0];
        hh_ext = [hh_ext; 0];
        flag1 = true;
    elseif cpr(end)~=1
        indicators = [indicators; 1];
        cpr = [cpr; 1];
        hh = [hh; 0];
        hh_ext = [hh_ext; 0];
        flag1 = true;
    else
        flag1 = false;
    end
    
    % find ids of boundaries that begin a "non hole" zone
    extid = find(~indicators(1:end-1));
    
    % do the same for the extension
    cprext = [centers_x-extension_factor*radiuses_x; centers_x+extension_factor*radiuses_x]; 
    indicatorsext = [ones(length(centers_x), 1); -ones(length(centers_x), 1)];
    [cprext, sortidext] = sort(cprext);
    indicatorsext = indicatorsext(sortidext);
    
    id1 = find(cprext>=0, 1, 'first');
    idend = find(cprext<=1, 1, 'last');
    cprext = cprext(id1:idend);
    indicatorsext = indicatorsext(id1:idend);
    
    indicatorsext = cumsum(indicatorsext)+id1-1;
    left = id1 - 1;
    if min(indicatorsext)<0
        indicatorsext = indicatorsext-min(indicatorsext);
    end
    
    % get rid of multiple centers
    [cprext, iduniqueext] = unique(cprext);
    indicatorsuniqueext = zeros(length(iduniqueext), 1);
    if ~isempty(iduniqueext)
        iduniqueext = [iduniqueext; iduniqueext(end)];
        for ii = 1:length(iduniqueext)-1
            if ii ~= length(iduniqueext)-1
                idsext = iduniqueext(ii):iduniqueext(ii+1)-1;
            else
                idsext = iduniqueext(ii):iduniqueext(ii+1);
            end
            indicatorsuniqueext(ii) = indicatorsext(idsext(end));
        end
    end
    indicatorsext = indicatorsuniqueext;
    
    % take care of boundaries
    if isempty(cprext) || cprext(1)~=0
        indicatorsext = [left; indicatorsext];
        cprext = [0; cprext];
    end
    if ~isempty(indicatorsext) && indicatorsext(end)==0 && cprext(end)~=1
        indicatorsext = [indicatorsext; 0];
        cprext = [cprext; 1];
    elseif cprext(end)~=1
        indicatorsext = [indicatorsext; 1];
        cprext = [cprext; 1];
    end
    
    % find ids of boundaries that begin a "non hole" zone
    extidext = find(~indicatorsext(1:end-1));
    extid_real = [];
    for index = 1:length(extidext)
        extid_temp = find(cpr <= cprext(extidext(index)));
        extid_real = [extid_real, extid_temp(end)];
    end

    knots_x = [];
    hh_icpr = hh;
    for icpr = 1:length(cpr)-1
        if (cpr(icpr)+hh(icpr) > cpr(icpr+1)-hh(icpr+1))
            hh([icpr, icpr+1]) = (cpr(icpr+1)-cpr(icpr))/2;
        end
    end
    
    for icpr = 1:length(cpr)-1
        if icpr == 1 && ~flag0
            knots_x = [knots_x, 0];
        end
        
        if ~ismember(icpr, extid) % inside features
            h_icpr = 2*min(hh_icpr([icpr, icpr+1])); 
            if h_icpr == 0 && icpr ~= length(hh_icpr)
                h_icpr = 2*max(hh_icpr([icpr, icpr+1])); 
            end
            knots_x = [knots_x, linspace(cpr(icpr)+hh(icpr), cpr(icpr+1)-hh(icpr+1), ...
                ceil( (cpr(icpr+1)-hh(icpr+1)-cpr(icpr)-hh(icpr))/h_icpr )+1)]; % h_feat
            
        elseif ~ismember(icpr, extid_real) % in the extension zone between features or at boundaries
            h_icpr = min(hh_ext([icpr, icpr+1])); 
            if h_icpr == 0
                h_icpr = max(hh_ext([icpr, icpr+1])); 
            end
            knots_x = [knots_x, linspace(cpr(icpr)+hh(icpr), cpr(icpr+1)-hh(icpr+1), ...
                ceil( (cpr(icpr+1)-hh(icpr+1)-cpr(icpr)-hh(icpr))/h_icpr )+1)]; % h_ext
            
        else % extension + outside zones between features or at boundaries
            ids_cprext = find(cprext>cpr(icpr)+hh(icpr) & cprext<cpr(icpr+1)-hh(icpr+1));
            if icpr == 1 && extid_real(1) == 1
                knots_x = [knots_x, linspace(cpr(icpr)+hh(icpr), cprext(ids_cprext(1)), ...
                    ceil( (cprext(ids_cprext(1))-cpr(icpr)-hh(icpr))/h_out )+1)];
            else
                h_icpr = hh_ext(icpr); 
                if h_icpr == 0 && icpr ~= length(hh_ext)
                    h_icpr = hh_ext(icpr+1); 
                elseif h_icpr == 0 && icpr ~= 1
                    h_icpr = hh_ext(icpr-1);
                end
                knots_x = [knots_x, linspace(cpr(icpr)+hh(icpr), cprext(ids_cprext(1)), ...
                    ceil( (cprext(ids_cprext(1))-cpr(icpr)-hh(icpr))/h_icpr )+1)]; % h_ext
            end
            knots_x = knots_x(1:end-1);
            
            for icprext = 1:length(ids_cprext)-1
                if ismember(ids_cprext(icprext), extidext) % outside zone
                    knots_x = [knots_x, linspace(cprext(ids_cprext(icprext)), cprext(ids_cprext(icprext)+1),...
                        ceil( (cprext(ids_cprext(icprext)+1)-cprext(ids_cprext(icprext)))/h_out )+1)]; 
                    knots_x = knots_x(1:end-1);
                    
                else % extension zone
                    h_icpr = hh_ext(icpr); 
                    if h_icpr == 0 && icpr ~= length(hh_ext)
                        h_icpr = hh_ext(icpr+1); 
                    elseif h_icpr == 0 && icpr ~= 1
                        h_icpr = hh_ext(icpr-1);
                    end
                    knots_x = [knots_x, linspace(cprext(ids_cprext(icprext)), cprext(ids_cprext(icprext)+1),...
                        ceil( (cprext(ids_cprext(icprext)+1)-cprext(ids_cprext(icprext)))/h_icpr )+1)]; % h_ext
                    knots_x = knots_x(1:end-1);
                    
                end
            end
            if ismember(ids_cprext(end), extidext)
                knots_x = [knots_x, linspace(cprext(ids_cprext(end)), cpr(icpr+1)-hh(icpr+1), ...
                    ceil( (cpr(icpr+1)-hh(icpr+1)-cprext(ids_cprext(end)))/h_out )+1)];
            else
                h_icpr = hh_ext(icpr+1); 
                if h_icpr == 0
                    h_icpr = hh_ext(icpr);
                end
                knots_x = [knots_x, linspace(cprext(ids_cprext(end)), cpr(icpr+1)-hh(icpr+1), ...
                    ceil( (cpr(icpr+1)-hh(icpr+1)-cprext(ids_cprext(end)))/h_icpr )+1)]; % h_ext
            end
        end
        
        if icpr+1 == length(cpr) && ~flag1
            knots_x = [knots_x, 1];
        end
    end
    
    knots{idim} = knots_x;
end

knots{1} = knots{1}(2:end-1);
knots{2} = knots{2}(2:end-1);

square = nrbsquare([0,0], 1, 1);
square = nrbdegelev(square, degree-1);
square = nrbkntins(square, knots);

geo_name = square; 

end
