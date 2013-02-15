%%
% CS 6640 : Image Processing Final Project
%
% Author : Arthur COSTE
% Date : November 2012
%
% Content : Euler Lagrange snake implementation
%
%%
clear all
close all

I=imread('test_connexite.tif');
I=double(I(:,:,1));
I2=I;
figure(1);
imagesc(I);
axis square;

%gaussian filtering
init = [1,1];
gauss = [1,1];
for i=1:11
    gauss = conv(gauss,init)/sum(2*gauss);
end
weight=gauss'*gauss;
for i = ceil(size(weight,1)/2) :1: size(I,1)-size(weight,1)+ceil(size(weight,1)/2)
    for j = ceil(size(weight,2)/2) :1: size(I,2)-size(weight,2)+ceil(size(weight,2)/2)
        convol=0;
        %compute convolution for the neighbourhood associated to the kernel
        for a = 1:size(weight,1)
            for b=1:size(weight,2) 
   
            convol = convol + (weight(a,b)*I2(i-a+ceil(size(weight,1)/2),j-b+ceil(size(weight,2)/2)));
            
            end
        end
        I(i,j)=convol;
    end
end     
I=I/max(I(:));                     



[Ix,Iy]=size(I);

%Get points of initialized contour and if necessary add more
[Vx,Vy]=select_points(I2);

%plot(Vx,Vy,'ok','LineWidth',2)

[X,Y]=size(Vx);
%closing contour
Vx(X)=Vx(1); 
Vy(X)=Vy(1);

%Second Order discrete derivation Matrix D_2
V=-2*ones(X,1);
W=1*ones(X-1,1);

D_2=diag(V)+diag(W,1)+diag(W,-1);
D_2(1,X)=1;
D_2(X,1)=1;

% Fourth Order discrete derivation Matrix D_4
V=6*ones(X,1);
W=-4*ones(X-1,1);
U=ones(X-2,1);

D_4=diag(V)+diag(W,1)+diag(W,-1)+diag(U,2)+diag(U,-2);
D_4(1,X)=-4;
D_4(X,1)=-4;
D_4(1,X-1)=1;
D_4(X-1,1)=1;
D_4(2,X)=1;
D_4(X,2)=1;

%Compute External Energy based on the gradient of the image gradient's norm
[TempX,TempY]=gradient(I);
nn=(TempX.^2+TempY.^2).^(1/2); %Compute norm
[Px,Py]=gradient(nn);   

% Build matrix A
dt=0.05; %Time step for temporal integration
alpha=1.0;
beta=1.0;
delta=900;
A=eye(X)+dt*(-alpha*D_2+beta*D_4);
A_inv=inv(A);

new_Vx=Vx;
new_Vy=Vy;

i=0;
% stopping criterion based on iterations
while(i~=650) 
    Inter_Px=interp2(Px,new_Vx,new_Vy); %contour interpolation
    Inter_Py=interp2(Py,new_Vx,new_Vy);
    old_Vx=new_Vx;
    old_Vy=new_Vy;
    new_Vx=A_inv*(old_Vx+dt*delta*Inter_Px); 
    new_Vy=A_inv*(old_Vy+dt*delta*Inter_Py); 
    i=i+1;
    NewI=[new_Vx,new_Vy];
    if(mod(i,10)==0) % speed up display
        clf
        imagesc(I2);
        hold on
        plot(new_Vx,new_Vy,'c','LineWidth',2);
        pause(0.3);
    end
end



