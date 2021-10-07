% normX = norm(X);
% maxiter = 150;
% figure
% semilogy(error_SVRG_half(1:maxiter)/normX)
% 
% hold on
% semilogy(error_SVRG_half_old(1:maxiter)/normX)
% 
% hold on
% semilogy(error_SVRG_half_pure(1:maxiter)/normX)
% 
% hold on
% semilogy(error_SVRG_double(1:maxiter)/normX)
% 
% hold on
% semilogy(error_SVRG_double_old(1:maxiter)/normX)
% 
% title('relative errors for [100,50,25]')
% legend('SVRG mixed','SVRG mixed old','SVRG half','SVRG double','SVRG double old')

normX = norm(X);
figure
semilogy(error_SGD/normX)

hold on
semilogy(error_half/normX)

hold on
semilogy(error_ADAM/normX)

hold on
semilogy(error_half/normX)
legend('SGD newG','SGD newG half','ADAM oldG')
title('relative errors for [140,70,35]')


