function [a_long,a_long_short] = weight_calculation(alpha)

fin = size(alpha,1);
taille = size(alpha,2);

a_long = max(alpha,zeros(fin,taille))'./sum(max(alpha,zeros(fin,taille))', 'omitnan');
a_short = max(-1.*alpha,zeros(fin,taille))'./sum(max(-1.*alpha,zeros(fin,taille))', 'omitnan');
a_long_short = a_long - a_short;

end

