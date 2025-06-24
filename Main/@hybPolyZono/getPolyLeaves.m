% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return the binary factors \xib \in {-1,1}^nGb that result in 
%       non-empty constrained zonotopes 
%   Syntax:
%       [leaves] = getLeaves(Z,optSolver)
%   Inputs:
%       Z - 1D, 2D, or 3D hybrid zonotope in HCG-Rep (hybZono object)
%       optSolver - solver options needed for mixed-integer linear propgrams
%   Outputs:
%       leaves - nGb x nLeaves matrix, each column denoting the \xi vector
%                corresponding to one of the nLeaves non-empty constrained
%                zonotopes
%   Notes:
%       If \xib denotes the i^th column of leaves, then the corresponding
%       non-empty constrained zonotope is of the form:
%       Z = { (c + Gb \xib) + Gc \xic | ||\xic||_inf <= 1, Ac \xi = b - Ab \xib }
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [leaves] = getPolyLeaves(obj,optSolver)

% Full tree when there are no constrains (nLeaves = 2^nGb)
if isempty(obj.b)
	combos = ff2n(obj.nGb);
	combos(combos==0) = -1;
	leaves = combos';
	return
end

% Continue if set has equality constraints
if nargin < 2 || isempty(optSolver)
    optSolver = solverOptions;
end

if ~strcmp(optSolver.milpSolver,'gurobi')
    error('Gurobi is currently required for this function.')
end

% Get up to 2^nGb solutions
optSolver.nSolutions = 2^obj.nGb;
optSolver.MIPFocus = 0; % Find feasible solutions

leaves = allSignedComb(obj.nGb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem data for mixed-integer linear program (MILP)
% Aeq = [obj.Ac 2*obj.Ab];
% beq = obj.b+obj.Ab*ones(obj.nGb,1);
% lb = -ones(obj.nGc+obj.nGb,1);
% ub =  ones(obj.nGc+obj.nGb,1);
% lb((obj.nGc+1):end) = 0; % Set lower bound on binary factors to zero
% vType(1:obj.nGc) = 'C';
% vType(obj.nGc+1:(obj.nGc+obj.nGb)) = 'B';
% [x,~,~] = solveMILP([],[],[],Aeq,beq,lb,ub,vType,optSolver);
% 
% % convert binaries back to {-1 1}
% leaves = (x((obj.nGc+1):end,:)-0.5)/0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%% MINLP %%%%%%%%%%%%%%%%%%%%%%%%
% beq = obj.b+obj.Ab*ones(obj.nGb,1);
% lb = -ones(obj.nGc+obj.nGb,1);
% ub =  ones(obj.nGc+obj.nGb,1);
% lb((obj.nGc+1):end) = 0; % Set lower bound on binary factors to zero
% vtype = repmat('C',1,obj.nGc+obj.nGb);   % 先全是 'C'
% vtype(obj.nGc+1:end) = 'B';    
% model.lb    = lb;
% model.ub    = ub;
% model.vtype = vtype;
% model.modelsense = 'max';      % 这里只是找可行解
% model.A     = sparse(0,obj.nGc+obj.nGb);     % 没有线性约束
% model.rhs   = [];
% model.f   = [];
% model.sense = '';
% colIdx = cell(1,size(obj.EC,2));
% for i = 1:size(obj.EC,2)
%     colIdx{i} = find(obj.EC(:,i) ~= 0);   % 行号列表
% end
% powM = obj.EC;                       %  n_var × n_mono
% 
% nonlinFun = @(vars) obj.Ac * prod( (vars(1:obj.nGc)'.^powM) , 1 ).' + 2* obj.Ab * vars(obj.nGc + (1:obj.nGb));
% gen.type  = 'func';
% gen.fun   = nonlinFun;   % 就是新 g(x)
% gen.fvars = 1:(obj.nGc+obj.nGb);  % 依赖的决策变量编号
% gen.name  = 'Aeq_zero';
% gen.sense = '=';         % 等式
% gen.rhs   = beq;           % g(x) = 0
% model.gencon = gen;
% 
% % m = numel(beq);		
% % gen(m) = struct; % 预分配		
% % for k = 1:m		
% %     gen(k).type  = 'func';
% %     gen(k).fun   = @(x) obj.Ac(k,:) * prod((x(1:obj.nGc)'.^powM), 1).' ...
% %                        + 2*obj.Ab(k,:) * x(obj.nGc+(1:obj.nGb));
% %     gen(k).fvars = 1:(obj.nGc+obj.nGb);
% %     gen(k).sense = '=';
% %     gen(k).rhs   = beq(k);
% %     gen(k).name  = sprintf('eq_%d',k);
% % end
% % model.gencon = gen;
% 
% params.OutputFlag    = 1;
% params.FuncNonlinear = 1;   % 必开
% params.NonConvex     = 2;   % 非凸全局
% params.Threads       = 1;
% if optSolver.nSolutions > 1
%             params.PoolSearchMode = 2;
%             params.PoolSolutions = min(optSolver.nSolutions,2e9);
% end
% % params.PoolSolutions = 4;
% % params.PoolGap        = 0;   % 0 = 丢掉所有劣解
% % params.PoolGapAbs     = 0;   % （两行都设最严格）
% % params.MIPFocus = optSolver.MIPFocus;
% % params.FeasibilityTol = 1e-9;
% % params.IntFeasTol     = 1e-9;
% % params.FuncTol        = 1e-9;   % ≥ Gurobi 11.0
% result = gurobi(model, params);
% x = [result.pool.xn];
% leaves = (x((obj.nGc+1):end,:)-0.5)/0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %   取出的 b = x(idxB,:)  ,   ξ_b = 2*b-1
% xi_b  = 2*x(idxB,:) - 1;
% xi_c  =   x(1:nGc,:);

end