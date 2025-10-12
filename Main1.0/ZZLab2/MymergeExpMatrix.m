function [id,E1,E2] = MymergeExpMatrix(id1,id2,E1,E2)
% mergeExpMatrix - Merge the ID-vectors of two polyZonotope objects
%                  and adapte the exponent matrices accordingly
%
% Syntax:
%    [id,E1,E2] = mergeExpMatrix(id1,id2,E1,E2)
%
% Inputs:
%    id1 - ID-vector of the first polynomial zonotope
%    id2 - ID-vector of the second polynomial zonotope
%    E1 - exponent matrix of the first polynomial zonotope
%    E2 - exponent matrix of the second polynomial zonotope
%
% Outputs:
%    id - merged ID-vector
%    E1 - adapted exponent matrix of the first polynomial zonotope
%    E2 - adapted exponent matrix of the second polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       25-June-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure uniqueness
[E1,id1] = removeRedundantIds(E1,id1);
[E2,id2] = removeRedundantIds(E2,id2);

L1 = length(id1);
L2 = length(id2);

% ID vectors are identical
if L1 == L2 && all(id1 == id2)

    id = id1;

% ID vectors not identical -> MERGE
else

    % merge the two sets
    id = id1;
    ind2 = zeros(size(id2));
    for i = 1:length(id2)
       ind = find(id == id2(i));
       if isempty(ind)
          id = [id;id2(i)];
          ind2(i) = length(id);
       else
          ind2(i) = ind;
       end
    end

    % construct the new exponent matrices
    L = length(id);

    E1 = [E1;zeros(L-L1,size(E1,2))];

    temp = zeros(L,size(E2,2));
    temp(ind2,:) = E2;
    E2 = temp;
end

% ------------------------------ END OF CODE ------------------------------



% function [id,E1,E2] = MymergeExpMatrix(id1,id2,E1,E2)
% % mergeExpMatrix - Merge the ID-vectors of two polyZonotope objects
% %                  and adapt the exponent matrices accordingly
% %
% %   Implements the mergeID operator consistent with the theory:
% %     id = [id1; H],  H = ids in id2 but not in id1 (stable order, unique)
% %     E1_new = [E1; zeros(|H|, size(E1,2))]
% %     E2_new maps rows by id into the merged id; rows for ids not present are zeros
% %     If id2 contains duplicates, their exponent rows must be identical;
% %     otherwise an error is thrown (theory requires equality).
% %
% % Inputs:
% %   id1,id2 : column or row vectors of IDs (integers)
% %   E1,E2   : exponent matrices of the two polynomial zonotopes
% %
% % Outputs:
% %   id      : merged ID vector
% %   E1,E2   : adapted exponent matrices aligned to 'id'
% %
% % Authors:   @zhenzhang32768 (refined by assistant)
% % Date:      2025-10-12
% 
% % ------------------------------ BEGIN CODE -------------------------------
% 
% % Normalize id vectors to column
% if isrow(id1), id1 = id1.'; end
% if isrow(id2), id2 = id2.'; end
% 
% % Ensure uniqueness per object (user-provided helper)
% [E1,id1] = removeRedundantIds(E1,id1);
% [E2,id2] = removeRedundantIds(E2,id2);
% 
% L1 = numel(id1);
% L2 = numel(id2);
% 
% % Trivial case: identical IDs (same order & values)
% if L1 == L2 && all(id1 == id2)
%     id = id1;
%     return
% end
% 
% % --- Build H = ids in id2 but not in id1; keep stable order & unique -----
% in1 = ismember(id2, id1);
% H_candidates = id2(~in1);
% % unique with 'stable' to keep the first appearance order in id2
% H = unique(H_candidates,'stable');
% 
% % Merged id as in the theory: [id1; H]
% id = [id1; H];
% L  = numel(id);
% 
% % --------- Adapt E1: append zeros for the new |H| rows at the bottom -----
% if size(E1,1) ~= L1
%     error('E1 has %d rows but id1 has %d elements.', size(E1,1), L1);
% end
% E1 = [E1; zeros(L - L1, size(E1,2))];
% 
% % ---------------- Adapt E2: map rows into merged id order -----------------
% if size(E2,1) ~= L2
%     error('E2 has %d rows but id2 has %d elements.', size(E2,1), L2);
% end
% 
% % Check duplicates in id2: rows for the same id must be identical
% if ~isempty(id2)
%     [uIds, ~, grp] = unique(id2,'stable');
%     for g = 1:numel(uIds)
%         rows = find(grp == g);
%         if numel(rows) > 1
%             block = E2(rows, :);
%             % All equal? (allow exact equality; adapt here if you need tol.)
%             if any(any(block ~= block(1,:)))
%                 error(['mergeID: inconsistent exponent rows for duplicated id %d in id2. ', ...
%                        'Rows associated with the same dependent factor must be identical to satisfy the theory.'], uIds(g));
%             end
%         end
%     end
% end
% 
% % Build new E2 with L rows; place each (unique) id2 row at its merged index
% E2_new = zeros(L, size(E2,2));
% if ~isempty(id2)
%     [lia, locb] = ismember(id2, id);   % locb(k) is the row in merged 'id' for id2(k)
%     if ~all(lia)
%         error('Internal error: some id2 elements not found in merged id.');
%     end
%     % For duplicates, we already validated equality; keep first occurrence.
%     [~, firstIdx] = unique(id2,'stable');
%     rows_keep = firstIdx;
%     E2_new(locb(rows_keep), :) = E2(rows_keep, :);
% end
% E2 = E2_new;
% 
% % ------------------------------- END CODE --------------------------------
% end
