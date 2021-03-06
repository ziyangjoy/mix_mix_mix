condnums_single = [1e7,1e8,1e9,1e10,1e18];
condnum_string = [string('1e7'),string('1e8'),string('1e9'),string('1e10'),string('1e18')];
sizes = [10,20,30,40,50,100,200];
colors = ['y','r','g','b','k'];


for i = 1:numel(condnums_single)
    error_list = [];
    for n = sizes
        A = gallery('randsvd',n,condnums_single(i),3);
        B = gallery('randsvd',n,condnums_single(i),2);
        real_value = kr3(A,B);
        halfA = half(A);
        halfB = half(B);
        half_value = kr3(halfA,halfB);
        error = norm(real_value - double(half_value),'fro');
        error_list(end+1) = error;
    end
    figure('Name','condition number = '+ string(condnum_string(i)) );
    plot(sizes,error_list,'Color',colors(i));
    hold on;
    %disp(error_list);
end
