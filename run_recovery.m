
%% Recovery using ADMM & L1 norm constraint

function [im2, recovery_pos] = run_recovery(A, dx, dy, dz, mask_record)

np = 40;
itr_max = 1;
N = 320;
rsz = 44;
zn = 21;
up = 3;
for j=1:itr_max
    
    cos_movement; %Simulation of the cos_trajectory
    im2=im2-min(im2(:));
    fprintf('Running L_1 norm constrained minimization...\n');
    [u1,p1,pp1,img_raw_new] = ADMM_lessfft_learning(im2, A,2000);% ADMM L1 norm constrained minimization
    lambda=0.1; %step size
    ratio=0.04;
    max_loop = 1;
    paratrj = step3_solver_ww(u1, A, dx,dy,dz,img_raw_new,max_loop,lambda,ratio);% L2 norm constrained confinement
    I = paratrj(:,end-3);
    recovery_pos = paratrj(I>0,end-2:end);%Pos2: Recovered position
end
end