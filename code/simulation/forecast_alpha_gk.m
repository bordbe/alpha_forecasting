function [alpha_gk,IC,C] = forecast_alpha_gk(rdt,signals,z_scores, rates)

fin = size(z_scores,3);
taille = size(z_scores,1);
factor = size(z_scores,2);

alpha_gk = NaN(fin,taille);
C = zeros(factor,factor,fin);
IC = zeros(fin,factor);

%for each month :
for t=2:fin
    
    %C matrix computation
    C(:,:,t) = corr(signals(:,:,t),'rows','complete');
    
    %IC vector computation
    for i=1:factor
        IC(t,i) = corr(rdt(t,:)'-rates(t,1),signals(:,i,t-1),'rows','complete');
    end
    
    %Alphas computation
    alpha_gk(t,:) = std(rdt(t,:)-rates(t,1),'omitnan')*IC(t,:)*inv(C(:,:,t))*z_scores(:,:,t-1)';

end


end

