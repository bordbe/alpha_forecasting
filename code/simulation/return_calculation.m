function [rdt_w] = return_calculation(rdt,w,gap)

fin = size(rdt,1);
rdt_w = zeros(1,fin);

for t=1+gap:fin
    rdt_w(1,t) = sum(rdt(t,:).*w(:,t-gap)','omitnan');
end

end

