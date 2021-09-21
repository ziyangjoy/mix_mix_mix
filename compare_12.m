%load aminoacids.mat
N = 3;
for size = [50]
%     X = randi(size,size/2,size*2,size);
%     X = randn(size/2,size*2,size);
%     X = tensor(X);

    s = [size/2,size*2,size,size/2,size/2];
    r = 10;
    A = cell(N,1);
    for i = 1:N
        A{i} = randn(s(i),r);
    end
    X = ktensor(A);
    
    
    iterCG = 10;
    iterSG = 50;
    step = 0.01;
    initial = cell(N,1);
    
    for i = 1:N
        initial{i} = randn(s(i),r);
    end
    
    [nA,e1,~] = SGD_outer2(0,initial,X,iterCG,iterSG,0.0001,r,step);
    [nA2,e2,~] = SGD_outer2(2,initial,X,iterCG,iterSG,0.0001,r,step);
%     [~,e1,~] = SGD_outer2(0,X,20,50,0.0001,20,1);
%     [~,e2,~] = SGD_outer2(2,X,20,50,0.0001,20,1);

%original setting
%     [~,e1,~] = SGD_outer2(0,X,20,50,0.0001,20,1);
%     [~,e2,~] = SGD_outer2(2,X,30,50,0.0001,20,1);

    normX = norm(X);
    
    figure
    plot(e1/normX)
%     semilogy(e1/normX)
    hold on
    plot(e2/normX)
%     semilogy(e2/normX)
    legend('half precision','double precision')
    
    xlabel('number of iterations')
    ylabel('error')
    title('number of iteration vs error, number of elements = ',size^3)
end
% load aminoacids.mat
% [~,e1,~] = SGD_outer2(0,X,20,50,0.0001,20,1);
% [~,e2,~] = SGD_outer2(2,X,30,50,0.0001,20,1);
% 
% figure
% plot(e1)
% hold on
% plot(e2)
% legend('half precision','double precision')
% 
% xlabel('number of iterations')
% ylabel('error')
% title('number of iteration vs error, number of elements = 61305')
% % % % 
% figure
% plot(T1,e1)
% hold on
% plot(T2,e2)
% legend('single precision','double precision')
% xlabel('time')
% ylabel('error')
% title('time vs error')