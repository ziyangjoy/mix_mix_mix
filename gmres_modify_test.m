N = ndims(X);
Xnlist = {};
Xcontent = X.data;
for n = 1:N
    Xnlist{end+1}=n_mode_unfold(Xcontent,n);
end
Uinit = cell(N,1);
for n = 1:N
    Uinit{n} = rand(size(X,n),20);
end
U = Uinit;
n = 1;
uk = 0;
uf = 2;
uw = 2;
ur = 2;
giters = 3;
gtol = 0.0001;
Kr_product = kr3_for_cp(U,n,uk);
Xn = Xnlist{n};
Xnt = Xn';
smallXn=Xnt./10;
Ant = getx(Kr_product,Xnt,uf,uw,ur,giters,gtol);
A = rand(300,20); 
B = rand(300,10);
x = getx(A,B,uf,uw,ur,giters,gtol);
error = norm(A*x-B)