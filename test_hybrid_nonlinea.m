%% 清理环境
clear; close all; clc;




%% 1) 定义系统与初始集
% 二次系统参数 Q
Q{1} = [0.1 -0.12;
    0.5 0];
Q{2} = [-0.1 0;
    0.5 0];
% 线性系统参数 M
M{1} = [1.2 -1;
    -1 0.1];
M{2} = [1 0;
    0 2];


M1 = [0.75,0.25;-0.25,0.75];
Bd1 = [-0.25;-0.25];
M2 = [0.75,-0.25;0.25,0.75];
Bd2 = [0.25;-0.25];

% 初始 CPZ P0
c0 = [-1;5];
G0 = [-0.75 -0.25 0;
    0.75 -0.25 0];
E0 = [1 0 1;
    0 1 1];
P0 = conPolyZono(c0, G0, E0);



%% 假设 P 是你的 conPolyZono 对象

% —— 3) 常数偏移 —— 把常数当成中心
c_const = [0; 0];       % 若 f 有常数项，可写成 zonotope(c_const,[])





%% 2) 定义半空间 S: x1 < h 及其补集
h = 10e-6;
H = [1 0];
S_lo = struct('A', H, 'b', h); % x1 <= h 片段
S_hi = struct('A', -H, 'b', -h); % -x1 <= -h 片段 (x1 >= h)
%% 3) 三步迭代：半空间交集 + 对应映射
steps = 3;
Reach = cell(steps+1,1);
Reach{1} = { P0 };
colors = lines(steps+1);
for k = 1:steps
    prevList = Reach{k};
    nextList = {};
    for i = 1:numel(prevList)
        P = prevList{i};
        % ——— 下半空间：x1 < h ———
        P_lo = my_and_(P, S_lo);
        % 如果 P_lo 为空，跳过
        if ~isemptyobject(P_lo)
            Z_quad = quadMap(P_lo, Q); % 用二次系统
            Z_lin = M1 * P_lo;
            Z_lo = Z_quad + Z_lin + conZonotope(c_const, [], []);
        
            nextList{end+1} = Z_lo;
        end
        % ——— 上半空间：x1 >= h ———
        P_hi = my_and_(P, S_hi);
        if ~isemptyobject(P_hi)
             Z_quad = quadMap(P_hi, Q); % 用二次系统
            Z_lin = M2 * P_hi;
            Z_hi = Z_quad + Z_lin + conZonotope(c_const, [], []);
            nextList{end+1} = Z_hi;
        end
    end
    Reach{k+1} = nextList;
end
%% 4) 分步绘图（带空集过滤）
colors = lines(steps+1);
for k = 1:(steps+1)
    figure;
    for i = 1:numel(Reach{k})
        P = Reach{k}{i};
        % 如果 P 是空集，跳过
        if isemptyobject(P)
            continue;
        end
        hold on
        % 否则正常绘制
        plot(P, [1,2], ...
            'Color', colors(k,:), ...
            'Splits', 6);
    end
end
legend(arrayfun(@(k) sprintf('Step %d',k-1), 0:steps, 'Uni',false), ...
    'Location','bestoutside');
xlabel('$x_{(1)}$','Interpreter','latex');
ylabel('$x_{(2)}$','Interpreter','latex');
title('Reachable Sets over 3 Steps with Half-space Switching');
enlargeAxis(1.05);
grid on;