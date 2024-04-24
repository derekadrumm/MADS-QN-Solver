%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Non-Selfintersecting NURBS --------------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines a constraint function to check that a nurbs curve with 
% given cpts is non-self-intersecting, by checking if the control polygon
% has a self intersection. If the control polygon has a self intersection,
% then it is likely that the spline will have a self intersection as well.
% (This is not necessary or sufficient! Just a good rule-of-thumb.)
%
% We check each mincomb-tuple of control points, and see if the corresponding
% control polygon has a self intersection. This idea comes from
% [3] "Zhaou and Zhu. Injectivity of NURBS curves. 2015."
%
% Output 1 if there is an intersection, or -1 if no intersection found.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = poly_selfintersection(cpts,m,mincomb)
% INPUT:
% cpts:    control points of spline as a nx1 collumn vector [x;y]
% m:       how many dimensions of cpts we have
% mincomb: minimum vertex combinations to check (between 4 and N)
%
% OUTPUT:
% out:     1 if intersection found, -1 otherwise

N = length(cpts)/m;
cpts_x = cpts(1:N);
cpts_y = cpts(N+1:2*N);
cpts = [cpts_x cpts_y]';

% search from N-1 to 4 vertex combinations to find intersections
i = N+1;
while i > mincomb
    i = i-1;
    % generate combinations of control points to check
    combos = nchoosek([1:N],i);

    maxcount = size(combos,1);
    count = 0;
    INTERSECT_FOUND = false;
    while (~INTERSECT_FOUND) && (count < maxcount)
        count = count + 1;
        row_i = combos(count,:);
        vx = cpts_x(row_i)';
        vy = cpts_y(row_i)';
        % search for intersections of sub-polygons
        % find all polygon selfintersections (includes vertices)
        [x,y] = polyxpoly(vx,vy,vx([2:end 1]),vy([2:end 1]));

        % find points which dont occur in vx,vy
        XY = [x y]';
        pts = [vx;vy];
        check = ismember(XY',pts','rows')';
        XYint = unique(XY(:,~check)','rows')';

        INTERSECT_FOUND = ~isempty(XYint);
    end
    
    % stop if find an intersection
    if INTERSECT_FOUND
        break
    end

end
% output 1 if intersection found, or -1 otherwise
out = 2*sum(INTERSECT_FOUND) - 1;


end

