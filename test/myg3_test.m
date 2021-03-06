condnums_single = [1e1,1e2,1e3,1e4,1e5,1e6,1e7];
condnum_string = [string('condition number = 1e1'),string('condition number = 1e2'),string('condition number = 1e3'),string('condition number = 1e4'),string('condition number = 1e5'),string('condition number = 1e6'),string('condition number = 1e7')];
iterations = [1,2,3]; % usually converge in three iterations 
colors = ['y','r','g','b','k']; 
n = 20; % size of the matrix

% each plot is number of iterations vs error under different condition
% numbers.
for i = 1:numel(condnums_single)
    error_list_double = [];
    error_list_half = [];
    for j = iterations
        A = gallery('randsvd',n,condnums_single(i),3);
        b = rand(n,1);
        real_value = myg3(A,b,2,2,2,j,0.0001);
        half_value = myg3(A,b,0,2,2,j,0.0001);
        error_double = norm(A*real_value-b,'fro');
        error_half = norm(A*half_value-b,'fro');
        error_list_double(end+1) = error_double;
        error_list_half(end+1)=error_half;
    end
    figure('Name','condition number = '+ string(condnum_string(i)) );
    plot(iterations,error_list_double,'Color','g');
    hold on;
    plot(iterations,error_list_half,'Color','b');
end
% green line: error of double precision;
% blue line: error of using half precision as factorization precision. 
