function y = gf_mat_mul(a,b,add_mat,mul_mat)
[ma, na] = size(a); % a=info_seq
[~, nb] = size(b);  % b=G, generator matrix
y = zeros(ma, nb);
for i = 1 : ma    % Loop over rows of info_seq
  for j = 1 : nb  % Loop over columns of generator matrix
  rw = a(i,:);   % Extract the i-th row of a
  cl = b(:,j);   % Extract the j-th col of b 
    for k=1 : na  % Loop over cols of info_seq
        % Perform GF multiplication and addition using lookup tables
    y(i,j) = add_mat(1+y(i,j), 1+mul_mat(1+rw(k),1+cl(k)));
    end
  end
end
