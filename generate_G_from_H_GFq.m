function [G, Hsys, col_perm] = generate_G_from_H_GFq(H, add_mat, mul_mat)
% ================================================================
%  generate_G_from_H_GFq
% ================================================================
%  Purpose:
%   Generate generator matrix G from a parity-check matrix H
%   over GF(2^p), using lookup tables for addition/multiplication.
%
%  Inputs:
%   H         : MxN parity-check matrix over GF(2^p) (values 0..q-1)
%   add_mat   : qxq lookup table for GF addition
%   mul_mat   : qxq lookup table for GF multiplication
%
%  Outputs:
%   G         : KxN generator matrix (H * G' = 0)
%   Hsys      : systematic version of H  [I_M | P]
%   col_perm  : column permutation applied to reach systematic form
%
%  Example:
%   load('arith_4.mat');
%   [G, Hsys, perm] = generate_G_from_H_GFq(h_GF4, add_mat, mul_mat);

% ================================================================

H = full(H);
[M, N] = size(H);
q = size(add_mat,1);
K = N - M;

% Finite field helper functions 
gf_add = @(a,b) add_mat(a+1,b+1);   % 0-indexed field
gf_mul = @(a,b) mul_mat(a+1,b+1);

% Precompute multiplicative inverse lookup
inv_vec = zeros(1,q);
inv_vec(1) = 0;
for a = 2:q
    row = mul_mat(a,:);
    inv_idx = find(row == 1, 1);
    inv_vec(a) = inv_idx - 1;
end
gf_inv = @(a) inv_vec(a+1);

% Gaussian elimination over GF(q) 
Hwork = H;
col_perm = 1:N;
row = 1;
for col = 1:N
    if row > M
        break;
    end
    % Find a pivot (nonzero) in current column from row..M
    nz = find(Hwork(row:end,col) ~= 0, 1);
    if isempty(nz)
        continue; % move to next column
    end
    piv_row = row + nz - 1;

    % Swap rows if needed
    if piv_row ~= row
        tmp = Hwork(row,:);
        Hwork(row,:) = Hwork(piv_row,:);
        Hwork(piv_row,:) = tmp;
    end

    % Normalize pivot to 1
    piv = Hwork(row,col);
    if piv ~= 1
        invp = gf_inv(piv);
        for j = 1:N
            Hwork(row,j) = gf_mul(Hwork(row,j), invp);
        end
    end

    % Eliminate all other rows in this column
    for r2 = 1:M
        if r2 == row, continue; end
        a = Hwork(r2,col);
        if a ~= 0
            for j = 1:N
                Hwork(r2,j) = gf_add(Hwork(r2,j), gf_mul(a, Hwork(row,j)));
            end
        end
    end

    row = row + 1;
    if row > M, break; end
end

% Attempt to find pivot columns (identity part) 
pivot_cols = [];
for i = 1:M
    col_idx = find(Hwork(i,:) == 1);
    if numel(col_idx) == 1
        pivot_cols(end+1) = col_idx;
    end
end

non_piv = setdiff(1:N, pivot_cols, 'stable');
new_order = [pivot_cols non_piv];
Hsys = Hwork(:, new_order);
col_perm = col_perm(new_order);

% Extract [I_M | P] and build G = [P' | I_K] 
P = Hsys(:, M+1:end);
Gp = zeros(K, N);
Gp(:, 1:M) = P.';
for i = 1:K
    Gp(i, M+i) = 1;
end

% Undo the column permutation
inv_perm = zeros(1, N);
inv_perm(col_perm) = 1:N;
G = Gp(:, inv_perm);

end
