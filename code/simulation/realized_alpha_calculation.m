function [alpha_realized] = realized_alpha_calculation(alpha_simul,a_simul_ls,gap)

fin = size(alpha_simul,1);
alpha_realized = zeros(1,fin);

for t=1+gap:fin
    alpha_realized(1,t) = sum(alpha_simul(t,:).*a_simul_ls(:,t-gap)','omitnan')';
end


end

