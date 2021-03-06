function p = kr3_all_except_mode(list,mode,prec)
N = length(list);
if mode == 0
    p = list{1};
    for i = 2:N
        p = kr3(p,list{i},prec);
    end
end
list_except_mode = {list{[1:mode-1,mode+1:N]}};
p = list_except_mode{1};
for i = 2:N-1
    p = kr3(p,list_except_mode{i},prec);
end


