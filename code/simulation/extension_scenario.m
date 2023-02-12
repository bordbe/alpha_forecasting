function [output] = extension_scenario(rdt_snp,rdt,high_IC_ext,high_C_ext,rolling,gap)

%In this function we compute the realized alpha with the 3 approaches
%(TD,GK BU and NZ BU)

fin = size(rdt,1);
taille = size(rdt,2);

%Simulation
[rdt_ext,alpha_ext,~,rdt_bench_ext,~,z_simul_ext,~] = simulation(rdt_snp,rdt,(high_IC_ext'*ones(1,fin))',high_C_ext + zeros(2,2,fin),taille,2,fin);

% GK
%alpha GK estimation
[alpha_simul_gk_ext,~] = estimate_calculation(rdt_ext,rdt_bench_ext, z_simul_ext, rolling);
%estimated weight computation
[~,a_simul_ls_gk_ext] = weight_calculation(alpha_simul_gk_ext);
%realized alpha GK computation
alpha_realized_gk_ext = mean(realized_alpha_calculation(alpha_ext(rolling+1:end,:),a_simul_ls_gk_ext,gap)*12*100);

%NZ 
%alpha NZ estimation
alpha_simul_nz_ext = squeeze(sum(z_simul_ext(:,:,rolling+1:end),2))';
%estimated weight computation
[~,a_simul_ls_nz_ext] = weight_calculation(alpha_simul_nz_ext);
%realized alpha NZ computation
alpha_realized_nz_ext = mean(realized_alpha_calculation(alpha_ext(rolling+1:end,:),a_simul_ls_nz_ext,gap)*12*100);

%TD : we compute separatly our 2-single factor performance
%single factor 1
alpha_simul_td_ext1 = squeeze(sum(z_simul_ext(:,1,rolling+1:end),2))';
[~,a_simul_ls_td_ext1] = weight_calculation(alpha_simul_td_ext1);
alpha_realized_td_ext1 = mean(realized_alpha_calculation(alpha_ext(rolling+1:end,:),a_simul_ls_td_ext1,gap)*12*100);

%single factor 2
alpha_simul_td_ext2 = squeeze(sum(z_simul_ext(:,2,rolling+1:end),2))';
[~,a_simul_ls_td_ext2] = weight_calculation(alpha_simul_td_ext2);
alpha_realized_td_ext2 = mean(realized_alpha_calculation(alpha_ext(rolling+1:end,:),a_simul_ls_td_ext2,gap)*12*100);

%Creation of our TD multifactor portfolio
alpha_realized_td_ext = (alpha_realized_td_ext1 + alpha_realized_td_ext2)/2;

%output
output = {alpha_realized_td_ext,alpha_realized_nz_ext,alpha_realized_gk_ext};

end

