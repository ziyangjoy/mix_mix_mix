%load aminoacids.mat
for size = [10]
    X = randi(size,size/2,size*2,size);
    X = tensor(X);
    [~,e1,~] = SGD_outer2(0,X,20,50,0.0001,20,1);
    [~,e2,~] = SGD_outer2(2,X,30,50,0.0001,20,1);
    
    figure
    plot(e1)
    hold on
    plot(e2)
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