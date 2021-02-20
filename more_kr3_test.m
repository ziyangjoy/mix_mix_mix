R = 10;
sizes = [10,20,30,40,50,60,70,80,90,100];
%matrix = zeros(1,10);
%matrix(1)=gallery('randsvd',n,100,1);    
errors1 = [];
errors2 = [];
for i = sizes
    A = rand(i,R);
    B = rand(10,R);
    hA = half(A);
    hB = half(B);
    realans = kr3(A,B);
    halfans = kr3(hA,hB);
    error = (norm(realans-double(halfans),'fro'));
    errors1(end+1) = error;
end
figure

plot(sizes,errors1)
title('Size vs error')
xlabel('row number of A')
ylabel('fro error')
disp(errors1)