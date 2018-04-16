%% Generate PSF matrix and 1st order derivative
function [A, dx, dy, dz] = gene_PSF_matrix(mask_record)

fprintf('Generating Matrices...\n')

[X, Y] = meshgrid(1:200, 1:200);
kdph = 1-0.5*((X-100.5).^2+(Y-100.5).^2)/100^2*(2.7/2/100)^2;% 2.7 mm beam size
rsz = 44;

i = 1:200;
[ii,jj] = meshgrid(i,i);
cir = sqrt((ii-100.5).^2+(jj-100.5).^2)<100;%Circular pupil aparture
N = 320;

zn = 21;
up = 3;
A = zeros(up*rsz, up*rsz, zn);

for n = 1:zn
    for j=1:1:1
        u2=fftshift(ifft2(cir.*mask_record(:,:,(n-1)*10+j),N*up,N*up));
        u2=u2(round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,:);
        temp=abs(u2).^2/10;
        A(:,:,n)=A(:,:,n)+temp;%A matrix contains 21 different layers
    end
end

phase_max = 0;
ri = 1/max(A(:));
A = A*ri;

delta = 0.5;
sft = delta/6.6545;

PSF_temp1 = psf3d_kernel_07302015(mask_record, kdph, zn, up, -sft, 0, 0, N, ri, phase_max);
PSF_temp2 = psf3d_kernel_07302015(mask_record, kdph, zn, up, sft, 0, 0, N, ri, phase_max);


PSF_dev_x = (PSF_temp2-PSF_temp1)/delta/2;
dx = PSF_dev_x(round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,:);

clear PSF_dev_x

PSF_temp1 = psf3d_kernel_07302015(mask_record, kdph, zn, up, 0, -sft, 0, N, ri, phase_max);
PSF_temp2 = psf3d_kernel_07302015(mask_record, kdph, zn, up, 0, sft, 0, N, ri, phase_max);


PSF_dev_y = (PSF_temp2-PSF_temp1)/delta/2;
dy = PSF_dev_y(round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,:);

clear PFS_dev_y


PSF_temp1 = psf3d_kernel_07302015(mask_record, kdph, zn, up, 0, 0, -5, N, ri, phase_max);
PSF_temp2 = psf3d_kernel_07302015(mask_record, kdph, zn, up, 0, 0, +5, N, ri, phase_max);


PSF_dev_z = (PSF_temp2-PSF_temp1)/delta/2;
dz = PSF_dev_z(round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,round(N/2*up)-rsz/2*up+1:round(N/2*up)+rsz/2*up,:);
%fdz = fftn(dz);
clear PSF_dev_z
clear PSF_temp*

end



