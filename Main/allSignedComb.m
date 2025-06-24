function M = allSignedComb(nGb)
% allSignedComb  Return an nGb × 2^nGb matrix of ±1 combinations
%
%   M = allSignedComb(nGb)  gives a matrix M whose columns enumerate every
%   length-nGb vector whose entries are –1 or 1.
%
%   Example
%   -------
%   >> allSignedComb(3)
%   ans =
%       -1 -1 -1 -1  1  1  1  1
%       -1 -1  1  1 -1 -1  1  1
%       -1  1 -1  1 -1  1 -1  1

    % --- 1)  0/1  -------------------------
    bits = dec2bin(0 : 2^nGb - 1, nGb); %  = 0/1 
    
    % --- 2)  ±1 ------------------------------
    M = 2 * (bits.' - '0') - 1;         % 
    
    % --- 3)  int8  -----------------------------
    % M = int8(M);   %  double 2 int8
end
