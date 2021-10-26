normX = norm(X(:));
figure
semilogy(error_double/normX);
hold on 
semilogy(min(error_w_all,error_all)/normX);
% legend('without SWA','with SWA')
legend('double','half')
title('relative errors for [50,50,50,50]')

xlabel('number of epoches')
ylabel('error')
% figure
% semilogy(error_SGD/normX)
% 
% hold on
% semilogy(error_half/normX)
% 
% hold on
% semilogy(error_ADAM/normX)
% 
% hold on
% semilogy(error_half/normX)
% legend('SGD newG','SGD newG half','ADAM oldG')
% title('relative errors for [140,70,35]')


