
function [Y,data,label,nr,nc,nb,L,N,anomaly_map,normal_map,f_show]=ppData(data,map)

[nr,nc,nb]=size(data);f_show=data(:,:,[37,18,8]);
for i=1:3;max_f=max(max(f_show(:,:,i)));min_f=min(min(f_show(:,:,i)));f_show(:,:,i)=(f_show(:,:,i)-min_f)/(max_f-min_f);end
figure;subplot_tight(1, 2, 1,[.08 .08]); imagesc(f_show);  title('RGM Map','fontsize',8);
subplot_tight(1, 2, 2,[.08 .08]); imagesc(map); title('Label','fontsize',8);
for i=1:nb;data(:,:,i) = (data(:,:,i)-min(min(data(:,:,i)))) / (max(max(data(:,:,i))-min(min(data(:,:,i)))));end
Y=reshape(data, nr*nc, nb)';[L,N]=size(Y);
label=reshape(map,1,nr*nc);anomaly_map = logical(double(label)>0);normal_map = logical(double(label)==0);


end