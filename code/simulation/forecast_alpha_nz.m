function [alpha_nz] = forecast_alpha_nz(z_scores)

fin = size(z_scores,3);
taille = size(z_scores,1);

alpha_nz = NaN(fin,taille);
%for each month :
for t=2:fin
    %Alphas computation
    alpha_nz(t,:)= nanmean(z_scores(:,:,t-1),2)';
end

end

