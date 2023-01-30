close all
clear
time=[];
dc=1;
for ii=1:length(dc)
    switch (num2str(dc(ii)))
        case '1';load sandiego;dctype='sandiego';
    end
    
    [Y,data,label,nr,nc,nb,L,N,anomaly_map,normal_map,f_show]=ppData(data,map);
   

    
    tic
    [Smafra]=runSaFraV1(Y,data,label,dctype);
    time=toc;
    [AUCpdpf,AUCpdtau,AUCpftau,Pf,Pd]=clcAUCv1(label,Smafra);
    
    %% -------Map

    Mmafra=reshape(Smafra',nr,nc,1);
    figure
    subplot_tight(1, 3, 1,[.003 .003]); imagesc(f_show,[0,1]);  axis image;axis off;
    subplot_tight(1, 3, 2,[.003 .003]); imagesc(map,[0,1]);  axis image;axis off;
    subplot_tight(1, 3, 3,[.003 .003]); imagesc(Mmafra,[0,1]);  axis image;axis off;

    %% -------Display
    
    figure
    semilogx(Pf, Pd, '-', 'LineWidth', 2);

    xlabel('False Alarm Rate'); ylabel('Probability of Detection');
    legend('SaFra','location','northwest')
    set(gca,'Fontname','arial')
    %axis([0,0.07,0.6,1])
    set(gca,'FontSize',12);
    grid on

end