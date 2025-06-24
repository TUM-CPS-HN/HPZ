% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a hybrid zonotope in 1D, 2D, or 3D by
%       treating hybrid zonotope at the union of constrained zonotopes
%   Syntax:
%       [v,f] = plotAsConZono(Z,optSolver)
%   Inputs:
%       Z - 1D, 2D, or 3D hybrid zonotope in HCG-Rep (hybZono object)
%       optSolver - solver options needed for linear and mixed-integer linear propgrams
%   Outputs:
%       v - nV x 3 matrix, each row denoting the x (first column), y (second column), 
%                          and z (third column) positions of the nV vertices
%       f - nF x nMax matrix, each row denoting the vertices (up to nMax) contained
%                          in the nF faces (padded with NaN if face
%                          contains less than nMax vertices)
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotAsConPolyZono(obj,optSolver)
% optSolver='gurobi';
% figure;
hold on;
% Determine number of non-empty constrained zonotopes
[leaves] = getPolyLeaves(obj,optSolver);
if isempty(leaves)
    warning('zonoLAB:EmptyZonotope','Hybrid zonotope is empty and cannot be plotted.')
    v = []; f = [];
    return
end
nLeaves = size(leaves,2);

optPlot = plotOptions('Display','off','SolverOpts',optSolver);
v = [];
f = [];
nVerts = zeros(nLeaves,1);
waitbarHandle = waitbar(0,['Plotting hybrid polynomial zonotope with ',num2str(nLeaves),' leaves.']);
%  If \xib denotes the i^th column of leaves, then the corresponding
%  non-empty constrained zonotope is of the form:
%  Z = { (c + Gb \xib) + Gc \xic | ||\xic||_inf <= 1, Ac \xi = b - Ab \xib }
emptyLeaves = false;
for i = 1:nLeaves

    isIdE = isequal( obj.E, eye(size(obj.E)) );
    isIdEC = isequal( obj.EC, eye(size(obj.EC)) );
    if isIdE && isIdEC          % && 
    % conZonotope(c, Gc, Ac, b)
    Zi = conZonotope( ...
             obj.c  + obj.Gb * leaves(:, i), ...
             obj.Gc, ...
             obj.Ac, ...
             obj.b  - obj.Ab * leaves(:, i) );
    % plot(Zi);
    else
    % conPolyZono(c, Gc, E, Ac, b, EC)
    Zi = conPolyZono( ...
             obj.c  + obj.Gb * leaves(:, i), ...
             obj.Gc, ...
             obj.E, ...
             obj.Ac, ...
             obj.b  - obj.Ab * leaves(:, i), ...
             obj.EC );
    % plot(Zi,[1,2],'Splits',16);
    end
    % if isIdE==1 & isIdEC==1
    %     Zi = conZonotope(obj.c+obj.Gb*leaves(:,i),obj.Gc,obj.Ac,obj.b-obj.Ab*leaves(:,i));  
    %     else
    % % Zi = conZono(obj.Gc,obj.c+obj.Gb*leaves(:,i),obj.Ac,obj.b-obj.Ab*leaves(:,i));
    %     Zi = conPolyZono(obj.c+obj.Gb*leaves(:,i),obj.Gc,obj.E,obj.Ac,obj.b-obj.Ab*leaves(:,i),obj.EC);
    %     end
    % % [vi,fi] = plot(Zi,optPlot);
    % % [vi,fi] = plotHPZ(Zi,optPlot);
    % end
    % redOp=40;
    % Zi=reduce(Zi,'girard',redOp);

    color={'r','b','g','c'};
    try
        % plot(Zi,[1,2],'Splits',8,'FaceColor',color{i}, 'EdgeColor','k', 'FaceAlpha',0.3);
        plot(Zi,[1,2],'Splits',8,'FaceColor','b', 'EdgeColor','k', 'FaceAlpha',0.3);
    catch
        
    end
        % plot(Zi);
    hold on;
    % if size(vi,1) == 0
    %     emptyLeaves = true;
    % end
    % nVerts(i) = size(vi,1);
    % v = [v;vi];
    % if size(fi,2) > size(f,2)
    %     f = [f nan(size(f,1),size(fi,2)-size(f,2))]; 
    % end
    % if size(fi,2) < size(f,2)
    %     fi = [fi nan(size(fi,1),size(f,2)-size(fi,2))];
    % end
    % f = [f;fi+sum(nVerts(1:i-1))];
    waitbar(i/nLeaves,waitbarHandle)
end
close(waitbarHandle)

% Check for unplotted conZono leaves
if emptyLeaves
    warning('zonoLAB:Tolerance','Some leaves of the hybrid zonotope did not plot and may be caused by constraints that are nearly redundant (close to the plotting/optimization tolerances). Check the validity of other leaves.')
end
% res=0;
v=0;
f=0;
end