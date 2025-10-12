function [id, E1, E2, R1, R2] = mergeID_HPZ(id1, id2, E1, E2, R1, R2)
% mergeID_CPZ - Merge the ID-vectors of two CPZs and adapt exponent/residual matrices
%   Implements the mergeID operator from the proposition:
%     id = [id1; H], where H are ids in id2 but not in id1 (stable, unique)
%     E1 <- [E1; zeros(|H|, h1)], R1 <- [R1; zeros(|H|, q1)]
%     E2 <- aligned to the merged 'id' (first occurrence rule; zero if absent)
%     R2 <- aligned to the merged 'id' (first occurrence rule; zero if absent)
%
% Inputs:
%   id1, id2 : ID vectors (row/column) for the two CPZs
%   E1, E2   : exponent matrices, sizes [length(id1) x h1], [length(id2) x h2]
%   R1, R2   : (optional) residual matrices, sizes [length(id1) x q1], [length(id2) x q2]
%
% Outputs:
%   id       : merged ID vector
%   E1, E2   : adapted exponent matrices aligned to 'id'
%   R1, R2   : adapted residual matrices aligned to 'id' (if not provided, returned as [])
%
% Notes:
%   - Duplicates inside one CPZ are first removed via removeRedundantIds(E, id).
%     为保持一致，R 会在去冗余后按新 id 自动对齐（取首次出现行，其他丢弃）。
%   - 若 id1 与 id2 完全相同（同序），直接返回；矩阵不改。
%
% Author:    zhen zhang & assistant
% Date:      2025-10-12

% ------------------------------ BEGIN CODE -------------------------------

% ---- Normalize shapes ----
if isrow(id1), id1 = id1.'; end
if isrow(id2), id2 = id2.'; end
if nargin < 5 || isempty(R1), R1 = []; end
if nargin < 6 || isempty(R2), R2 = []; end

% ---- Basic size checks ----
L1_in = numel(id1);
L2_in = numel(id2);
assert(size(E1,1) == L1_in, 'E1 rows (%d) must match length(id1) (%d).', size(E1,1), L1_in);
assert(size(E2,1) == L2_in, 'E2 rows (%d) must match length(id2) (%d).', size(E2,1), L2_in);
if ~isempty(R1), assert(size(R1,1) == L1_in, 'R1 rows (%d) must match length(id1) (%d).', size(R1,1), L1_in); end
if ~isempty(R2), assert(size(R2,1) == L2_in, 'R2 rows (%d) must match length(id2) (%d).', size(R2,1), L2_in); end

% ---- Keep copies for later re-alignment of R ----
id1_orig = id1;
id2_orig = id2;

% ---- Per-object uniqueness (remove redundant ids inside each CPZ) ----
%     用户环境已有 removeRedundantIds(E, id)
[E1, id1] = removeRedundantIds(E1, id1);
[E2, id2] = removeRedundantIds(E2, id2);

% After id1/id2 changed, align R1/R2 accordingly (if provided)
if ~isempty(R1)
    R1 = alignByMergedId(id1, id1_orig, R1);  % 取首次出现的行；无则零
end
if ~isempty(R2)
    R2 = alignByMergedId(id2, id2_orig, R2);
end

L1 = numel(id1);
L2 = numel(id2);

% ---- Trivial case: identical ID order ----
if L1 == L2 && all(id1 == id2)
    id = id1;
    return
end

% ---- Build H = ids in id2 but not in id1; stable order, unique ----
maskNew = ~ismember(id2, id1);
H = unique(id2(maskNew), 'stable');

% ---- Merged id ----
id = [id1; H];
L  = numel(id);

% ---- Adapt E1 / R1: append |H| zero rows ----
if size(E1,1) ~= L1
    error('Internal: E1 rows changed unexpectedly.');
end
E1 = [E1; zeros(L - L1, size(E1,2))];

if ~isempty(R1)
    assert(size(R1,1) == L1, 'Internal: R1 rows changed unexpectedly.');
    R1 = [R1; zeros(L - L1, size(R1,2))];
end

% ---- Adapt E2 / R2: align to merged id using first-occurrence rule ----
E2 = alignByMergedId(id, id2, E2);   % 合并对齐（出现则取首行；否则零）
if ~isempty(R2)
    R2 = alignByMergedId(id, id2, R2);
else
    R2 = [];
end

% ------------------------------ END OF CODE ------------------------------
end


% ======================== Helper: alignByMergedId =========================
function M_out = alignByMergedId(id_merged, id_part, M_part)
% Align a per-id matrix M_part (rows correspond to id_part) to id_merged.
% Rule:
%   - 对于 id_merged 中的每个 κ：
%       若 κ 在 id_part 中出现，取其在 id_part 的第一次出现行复制到 M_out 对应行；
%       否则该行置零。
%   - 若 id_part 内部有重复，我们自然采用“首出现”行（与命题允许的 任一 j 一致）。
%
% 输入:
%   id_merged : 合并后的 id（列向量）
%   id_part   : 原矩阵对应的 id（列向量）
%   M_part    : 对应矩阵（行数 = numel(id_part)）
%
% 输出:
%   M_out     : 行数 = numel(id_merged) 的对齐矩阵
%
if isrow(id_part), id_part = id_part.'; end
L = numel(id_merged);
M_out = zeros(L, size(M_part, 2));

if isempty(id_part) || isempty(M_part)
    return;
end

% 为每个唯一 id（按出现顺序）选择其首次出现行
u = unique(id_part, 'stable');
for k = 1:numel(u)
    iMerged = find(id_merged == u(k), 1, 'first');
    if ~isempty(iMerged)
        iPart = find(id_part == u(k), 1, 'first');
        M_out(iMerged, :) = M_part(iPart, :);
    end
end
end
