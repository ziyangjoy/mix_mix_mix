function [L,U,p] = lu3(A)
[n,n] = size(A);
p = (1:n)';

for k = 1:n-1

   % Find index of largest element below diagonal in k-th column
   [r,m] = max(abs(double(A(k:n,k))));
   m = m+k-1;

   % Skip elimination if column is zero
   if (A(m,k) ~= 0)
   
      % Swap pivot row
      if (m ~= k)
         A([k m],:) = A([m k],:);
         p([k m]) = p([m k]);
      end

      % Compute multipliers
      i = k+1:n;
      A(i,k) = A(i,k)/A(k,k);

      % Update the remainder of the matrix
      j = k+1:n;
      A(i,j) = A(i,j) - A(i,k)*A(k,j); 
   end
end
A = double(A);
% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);
