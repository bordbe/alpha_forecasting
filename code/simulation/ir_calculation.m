function [ir] = ir_calculation(alpha_realized,a_simul_ls, theta,gap)

fin = size(alpha_realized,2);
ir = zeros(1,fin);

for t=1+gap:fin
    ir(1,t) = alpha_realized(1,t)/sqrt((a_simul_ls(:,t-gap))'*theta*(a_simul_ls(:,t-gap)));
end



end

