function [alpha_simul_gk,IC_estimate] = estimate_calculation(rdt_simul,rdt_bench, z_simul, rolling)

taille = size(z_simul,1);
fin = size(z_simul,3);
factor=size(z_simul,2);

k_simul = zeros(fin,factor);
C_simul = zeros(fin,factor,factor);
residuals_return = zeros(fin,taille);

%Calcul des vecteurs k et matrices C pour chaque periode
for t=2:fin
    for i=1:taille
        X = rdt_bench(1,max(1,t-rolling+1):t)';
        Y = rdt_simul(max(1,t-rolling+1):t,i);
        B = (X'*X)\(X'*Y);
        residuals_return(t,i) = rdt_simul(t,i) - (rdt_bench(1,t).*B)';
    end
    k_simul(t,:) = corr(residuals_return(t,:)',squeeze(z_simul(:,:,t-1)),'rows','complete');
    C_simul(t,:,:) = corr(z_simul(:,:,t),'rows','complete');


end

ICm = movavg(k_simul,'simple',rolling);
Cm = movavg(C_simul,'simple',rolling);

IC_estimate = zeros(fin-rolling, factor);
alpha_simul_gk = zeros(fin-rolling,taille);

for t=(rolling+2):fin
    IC_estimate(t-rolling,:)=ICm(t,:);
    alpha_simul_gk(t-rolling-1,:)=ICm(t,:)*inv(squeeze(Cm(t,:,:)))*z_simul(:,:,t-1)';
end

end

