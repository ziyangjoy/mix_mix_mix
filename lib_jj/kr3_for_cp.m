function pro = kr3_for_cp(Uinit,mode,prec)
U = {Uinit{[1:mode-1,mode+1:end]}};
N = length(U);
Uinv = {};
for i = 1:N
    Uinv{end+1} = U{end-i+1};
end
pro = Uinv{1};
for i = 2:N
    pro = kr3(pro,Uinv{i},prec);
end
