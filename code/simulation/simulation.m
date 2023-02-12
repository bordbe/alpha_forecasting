function [rdt_simul,alpha_simul,beta,rdt_bench,epsilon,z_simul,theta] = simulation(rdt_snp,rdt,IC,C,taille,factor,fin)

%Calcul par regression des returns residuels
residuals = NaN(fin,taille);
for i=1:taille
    X = rdt_snp';
    Y = rdt(:,i);
    B = (X'*X)\(X'*Y);
    residuals(:,i) = rdt(:,i) - (rdt_snp.*B)';
end

%Création de la matrice de covariance theta (le trident du papier)
theta = cov(residuals);
theta(isnan(theta))=0;

%Génération des epsilons
epsilon = mvnrnd(zeros(1,taille),theta,fin);

%Génération du benchmark par Maximum Likelyhood
max_likly = mle(rdt_snp);
rdt_bench = normrnd(max_likly(1), max_likly(2),1,fin);
beta = linspace(0.5,1.5,taille);

Ccs = zeros(factor,factor,fin);
z_simul = zeros(taille,factor,fin);
%Simulation des z_scores
for t=1:fin
    Ccs(:,:,t) = C(:,:,randi([2 fin],1));
    z_simul(:,:,t) = mvnrnd(zeros(1,factor),Ccs(:,:,t),taille);
end

alpha_simul = zeros(fin,taille);
rdt_simul = zeros(fin,taille);
%Calcul du alpha et des rendements simulés
for t=2:fin
    alpha_simul(t,:) = std(epsilon(t,:),0,2,'omitnan')*IC(randi(fin,1),:)*inv(Ccs(:,:,t))*z_simul(:,:,t-1)';
    rdt_simul(t,:)=alpha_simul(t,:)+beta*rdt_bench(1,t)+epsilon(t,:);
end


end