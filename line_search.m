function alpha = line_search(U,G,X,p)

alpha = 0.1/0.8;
c = 0.8; 


normG = sum(cellfun(@(x,y)norm(x(:)/y),G,p).^2);
T = double(tensor(ktensor(U)));
f = norm(T(:)-X(:))^2;
f_new = 1e+8;

while alpha>=0.001&& f-f_new <= 0.8*alpha*normG
   alpha = alpha*c;
   U_tmp = cellfun(@(x,y,z)x-alpha*y/z,U,G,p,'UniformOutput',false);
   T_tmp = double(tensor(ktensor(U_tmp)));
   f_new = norm(T_tmp(:)-X(:))^2;
end