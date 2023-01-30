function [AUC,AUC_PDtau,AUC_PFtau,Pf,Pd]=clcAUCv1(ano_map,Sinput)



[Pf,Pd,Tau] = perfcurve(ano_map,Sinput,'1') ;

AUC         = -sum((Pf(1:end-1)-Pf(2:end)).*(Pd(2:end) + Pd(1:end-1))/2)
AUC_PDtau   = sum((Tau(1:end-1)-Tau(2:end)).*(Pd(2:end)+Pd(1:end-1))/2)
AUC_PFtau   = sum((Tau(1:end-1)-Tau(2:end)).*(Pf(2:end)+Pf(1:end-1))/2)


tau2 = sort(Sinput,'descend');
map1 =Sinput;
for i = 1: length(tau2)
    map1(Sinput>=tau2(i))=1;
    map1(Sinput<tau2(i))=0;
    [Pd(i),Pf(i)] = cal_pdpf(map1,ano_map);
end

end

function [PD,PF] = cal_pdpf(result,GT)

BKG = find(GT == 0);
g_result = find(GT == 1);
d_result = find(result(:) == 1);
PD(:) = length(intersect(d_result,g_result))/length(g_result);
PF(:) = length(intersect(d_result,BKG))/length(BKG);
end