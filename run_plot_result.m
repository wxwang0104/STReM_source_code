%% Overlay the recovery with the raw iamge

rec_pos = sortrows(rec_pos,3);
pos=rec_pos;
tt=rec_pos(:,3);
t1=tt/210*100;
pos(:,3)=t1;
pos=sortrows(pos,1);
flag=1; 
figure;plot(1);colorbar;
ab=[];

t11=1-pos(:,3)/max(pos(:,3));
for i=1:1:numel(t11)
    f=t11(i);
    cm = colormap; % returns the current color map
    colorID = max(1, sum(f > [0:1/length(cm(:,1)):1])); 
    myColor = cm(colorID, :); % returns color
    ab=[ab;myColor];%store the color information of the emitters
end
close;

figure(101);imagesc(im_raw);colormap gray;hold on;axis square;

params2=pos;
params2(:,1:2) = params2(:,1:2)/3;
for i=1:1:size(ab,1)
    plot(params2(i,1),params2(i,2),'+','Markersize',7,'linewidth',3,'Color',[ab(i,1)*1,ab(i,2)*1,ab(i,3)]*1);
end
title('Overlaid of recovery on raw image')

