%% simulation of the continuous movement
fprintf('Simulating a cosin wave trajectory...\n')
close all
[X, Y] = meshgrid(1:200, 1:200);
rm = sqrt((X-100.5).^2+(Y-100.5).^2)<=100;
kdph = 1-0.5*((X-100.5).^2+(Y-100.5).^2)/100^2*(2.7/2/100)^2;
m_z=2;
np=210;
i = 1:200;
[ii,jj] = meshgrid(i,i);
cir = sqrt((ii-100.5).^2+(jj-100.5).^2)<100;

ii = poissrnd(1500*1, 1, np);
x_m=0.05*sin((1:np)*2*pi/np)*1;
y_m = 0.1*0.5*ones(1,np);
z_m=(rand(1, np)-0.5)*10;

x=[];y=[];z=[];
x(1)=0;
y(1)=-5.0;
z(1)=105;

for i=2:1:np
    x(i)=x(i-1)+x_m(i);
    y(i)=y(i-1)+y_m(i);
    z(i)=z(i-1)+z_m(i);
end
z=1:210;

x_record=x;
y_record=y;
z_record=z;
pos1=[x;y;z/max(abs(z))*18]';

im = zeros(N, N);
imm = zeros(N, N);
immm = zeros(N,N);

noise_posi = zeros(rsz,rsz);
for n = 1:numel(z)
    ph_xy = abs(X-100.5+1i*(Y-100.5)).*abs(x(n)+1i*y(n))*7/100.*cos(angle(X-100.5+1i*(Y-100.5))-angle(x(n)+1i*y(n)));
    t2 = fftshift(ifft2(cir.*mask_record(:,:,ceil(mod(z(n),211))).*exp(1i*ph_xy), N,N));
    I = abs(t2).^2/10;
    I = I/sqrt(sum(I(:).^2))*sqrt(ii(n))/1*sqrt(98);
    im = im+I;   
    noise_posi=noise_posi+poissrnd(50, rsz,rsz)*1;

end

im2 = im(N/2-rsz/2+1:N/2+rsz/2,N/2-rsz/2+1:N/2+rsz/2);
im2 = im2+noise_posi;
zn = 21;
pos1_real = [-x(1:1:end)'*6.6545*31/32*(N/200)+68, -y(1:1:end)'*6.6545*31/32*(N/200)+68, -z(1:1:end)'*(zn-1)/m_z+ceil(zn/2)];% 74 for rsz = 56 
close all;
figure;imagesc(im2);axis square;colormap gray;
title('Raw image of trajectory')
pos1_real(:,3)=0+1*pos1(:,3)/max(abs(pos1(:,3)))*21;


