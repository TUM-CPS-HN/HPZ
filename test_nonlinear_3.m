%% 清理环境
clear; close all; clc;
%% 1) 定义系统与初始集
% 二次系统参数 Q
Q{1} = 0.1*[-0.17 0.028;
    0 -0.17];
Q{2} = [-0.1 0;
    0 -0.1];
% 线性系统参数 M
M{1} = -Q{1};
M{2} = Q{2};
M1 = [0.75,0.25;-0.25,0.75];
Bd1 = [-0.25;-0.5];
M2 = [0.75,-0.25;0.25,0.75];
Bd2 = [0.25;-0.5];
% 初始 CPZ P0
c0 = [ -0.201 ; 1.395 ]; % Crosses the guard once
G0 = [0.2 0 0;
    0 0.2 0];
E0 = [1 0 1;
    0 1 1];
P0 = conPolyZono(c0, G0, E0);

steps = 10;
colors = interp1([1;steps+1],[0 0 1;1 0 0],1:1:steps+1);
cx = [ -0.201 ; 1.395  ];	% Crosses the guard once
Gx = 0.2*eye(2);
Z0 = zonotope(cx,Gx);	
% Propagate backwards two st
figure;
plot(Z0, [1,2],'FaceColor', colors(1,:), ...      % 填充色
     'FaceAlpha', 0.3, ...      
     'EdgeColor', 'k', ...     
     'EdgeAlpha', 0.5) % Plot initial condition set

%% 假设 P 是你的 conPolyZono 对象
% —— 3) 常数偏移 —— 把常数当成中心
c_const = [0; 0]; % 若 f 有常数项，可写成 zonotope(c_const,[])
%% 2) 定义半空间 S: x1 < h 及其补集
h = 10e-6;
H = [1 0];
S_lo = struct('A', H, 'b', h); % x1 <= h 片段
S_hi = struct('A', -H, 'b', -h); % -x1 <= -h 片段 (x1 >= h)
%% 3) 三步迭代：半空间交集 + 对应映射

Reach = cell(steps+1,1);
numreduce = 200;
 P0 = reduce(P0,'girard',numreduce);
Reach{1} = { P0 };
% colors = lines(steps+1);


   Z_quad = quadMap(P0, Q); % 用二次系统
            Z_lin = M1 * P0;
            Z_2 = Z_quad + Z_lin + conZonotope(Bd1, [], []);
            Z_2 = reduce(Z_2,'girard',numreduce);
            Reach{2} ={Z_2};

P_lo = my_and_(Z_2, S_lo);
P_hi = my_and_(Z_2, S_hi);


  Z_quad1 = quadMap(P_lo, Q); % 用二次系统
            Z_lin1 = M1 * P_lo;
            Z_lo1 = Z_quad1 + Z_lin1 + conZonotope(Bd1, [], []);
             Z_lo1 = reduce(Z_lo1,'girard',numreduce);
            nextList{1} = Z_lo1;
             
            
             Z_quad2 = quadMap(P_hi, Q); % 用二次系统
            Z_lin2 = M2 * P_hi;
            Z_hi2 = Z_quad2 + Z_lin2 + conZonotope(Bd2, [], []);
             Z_hi2 = reduce(Z_hi2,'girard',numreduce);
            nextList{end+1} = Z_hi2;

            Reach{3} = nextList;

for k = 3:steps
    prevList = Reach{k};
    nextList = {};
        % ——— 下半空间：x1 < h ———
            P_lo = Reach{k}{1};
            Z_quad = quadMap(P_lo, Q); % 用二次系统
            Z_lin = M1 * P_lo;
            Z_lo = Z_quad + Z_lin + conZonotope(Bd1, [], []);
            Z_lo = reduce(Z_lo,'girard',numreduce);
            nextList{end+1} = Z_lo;
    
        % ——— 上半空间：x1 >= h ———
            P_hi = Reach{k}{2};
            Z_quad = quadMap(P_hi, Q); % 用二次系统
            Z_lin = M2 * P_hi;
            Z_hi = Z_quad + Z_lin + conZonotope(Bd2, [], []);
            Z_hi = reduce(Z_hi,'girard',numreduce);
            nextList{end+1} = Z_hi;
    Reach{k+1} = nextList;
end
%% 4) 分步绘图（带空集过滤）


for k = 2:(steps+1)
    for i = 1:numel(Reach{k})
        P = Reach{k}{i};
        % 如果 P 是空集，跳过
        if isemptyobject(P)
            continue;
        end
        hold on
        % 否则正常绘制
     try   plot(P, [1,2], ...
                'FaceColor', colors(k,:), ...      % 填充色
     'FaceAlpha', 0.5, ...      
     'EdgeColor', 'k', ...     
     'EdgeAlpha', 0.5, ...      
     'Splits', 16);
            % plot(P, [1,2], ...
            % 'Color', colors(k,:));

     end
    end
end
p = plot([ 0 0 ] , [ -2 2 ], 'g--','linewidth',2); % Guard
legend(p,{'Guard'},'interpreter','latex')
legend(arrayfun(@(k) sprintf('Step %d',k), 0:steps, 'Uni',false), ...
    'Location','bestoutside');
xlabel('$x_{(1)}$','Interpreter','latex');
ylabel('$x_{(2)}$','Interpreter','latex');
title('Reachable Sets over 3 Steps with Half-space Switching');
enlargeAxis(1.05);
grid on;