function [paratrj] = step3_solver_ww(u1, A, dx, dy, dz, im2, max_it,lambda,ratio)
[center, elx, ely, elz] = local_max(u1);
[a,b,c] = size(center);
elx = elx(center>0);
ely = ely(center>0);
elz = elz(center>0);
I = center(center>0);
q = find(center>0);
[iy,ix,iz] = ind2sub(size(center), q);
figure(101);imagesc(zeros(size(u1,1)/3,size(u1,2)/3));colormap gray;
hold on;
plot(ix/3,iy/3,'ro');
n = numel(I);
paratrj = zeros(n, 4*max_it+4);
paratrj(:,1:4) = [I, elx+ix, ely+iy, elz+iz];
m = numel(im2(:));
up = 3;
img = zeros(size(center(:,:,1)));
img(ceil(up/2):up:end, ceil(up/2):up:end) = im2;
X = zeros(m*up*up,4*n);
for i = 1:numel(I)
    ta = zeros(size(A,1)+30,size(A,2)+30);
    ta(16:end-15,16:end-15)=A(:,:,end+1-iz(i));
    

    point = zeros(size(ta)+0);
    point(iy(i)+15,ix(i)+15) = 1;
    ta = fftshift(ifft2(fft2(ta).*fft2(point)));
    ta = ta(16:end-15,16:end-15);
    X(:,i) = ta(:);

    ta_x = zeros(size(A,1)+30,size(A,2)+30);
    ta_x(16:end-15,16:end-15)=-dx(:,:,end+1-iz(i));
    ta_x = fftshift(ifft2(fft2(ta_x).*fft2(point)));
    ta_x = ta_x(16:end-15,16:end-15);
    X(:,n+i) = ta_x(:);

    ta_y = zeros(size(A,1)+30,size(A,2)+30);
    ta_y(16:end-15,16:end-15)=-dy(:,:,end+1-iz(i));
    ta_y = fftshift(ifft2(fft2(ta_y).*fft2(point)));
    ta_y = ta_y(16:end-15,16:end-15);
    X(:,2*n+i) = ta_y(:);

    ta_z = zeros(size(A,1)+30,size(A,2)+30);
    ta_z(16:end-15,16:end-15)=-dz(:,:,end+1-iz(i));
    %ta_z = -dz(:,:,end+1-iz(i));
    ta_z = fftshift(ifft2(fft2(ta_z).*fft2(point)));
    ta_z = ta_z(16:end-15,16:end-15);
    X(:,3*n+i) = ta_z(:);

end

paras = (X'*X)\(X'*img(:));
rank(X'*X);
paras = reshape(paras, n, 4);
paras(:,2) = paras(:,2)./paras(:,1);
paras(:,3) = paras(:,3)./paras(:,1);
paras(:,4) = paras(:,4)./paras(:,1);

paratrj(:,5:8) = [paras(:,1), lambda*paras(:,2)+ix, lambda*paras(:,3)+iy, lambda*paras(:,4)+iz];

paras = paratrj(:,5:8);

    ix = paras(:,2);
    iy = paras(:,3);
    iz = paras(:,4);

%
%% usually not necessary
all_n = n;
all_inx = 1:all_n;
I = paras(:,1);
%ratio = 0.1;
thread = max(I)*ratio;
%max_it = 10;
for s = 2:max_it
    inx1 = all_inx(iz>c|iz<1|ix>b|ix<1|iy>a|iy<1 | paras(:,1)<=thread);
    all_inx(inx1) = 0;
    paras(inx1, :) = [];

    
    ix = paras(:,2);
    iy = paras(:,3);
    iz = paras(:,4);
    
    n = numel(round(iy));
    X = zeros(m*up*up,4*n);
    for i = 1:n
        if abs(round(iz))<100
        ta = A(:,:,end+1-round(iz(i)));
        point = zeros(size(ta));
        point(round(iy(i)),round(ix(i))) = 1;
        ta = fftshift(ifft2(fft2(ta).*fft2(point)));
        X(:,i) = ta(:);
        %
        ta_x = -dx(:,:,end+1-round(iz(i)));
        ta_x = fftshift(ifft2(fft2(ta_x).*fft2(point)));
        X(:,n+i) = ta_x(:);

        ta_y = -dy(:,:,end+1-round(iz(i)));
        ta_y = fftshift(ifft2(fft2(ta_y).*fft2(point)));
        X(:,2*n+i) = ta_y(:);

        ta_z = -dz(:,:,end+1-round(iz(i)));
        ta_z = fftshift(ifft2(fft2(ta_z).*fft2(point)));
        X(:,3*n+i) = ta_z(:);
        end
    end
    paras = (X'*X)\(X'*img(:));
    paras = reshape(paras, n, 4);
    paras(:,2) = paras(:,2)./paras(:,1);
    paras(:,3) = paras(:,3)./paras(:,1);
    paras(:,4) = paras(:,4)./paras(:,1);

    
    paratrj(all_inx(all_inx>0), s*4+1:s*4+4)=[paras(:,1),lambda*paras(:,2)+ix,lambda*paras(:,3)+iy,lambda*paras(:,4)+iz];    
   
    paras = paratrj(:,s*4+1:s*4+4);
    I = paras(:,1);

    thread = max(I)*ratio;

    ix = paras(:,2);
    iy = paras(:,3);
    iz = paras(:,4);
    all_inx = 1:all_n;
    figure(101);imagesc(zeros(size(u1,1)/3,size(u1,2)/3));colormap gray;
    hold on;
    plot(ix/3,iy/3,'ro');pause(0.0001);
end
