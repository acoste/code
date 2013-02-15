%%
% CS 6640 : Image Processing Project 2
%
% Author : Arthur COSTE
% Date : October 2012
%
% Content : Non linear warping using Radial basis function
%
%%
clear all
close all

sigma = 8;

I=imread('grid_2.tif');
I2=double(I(:,:,1));
newI = zeros(1*size(I2,1),1*size(I2,2));
interpol=1
figure(121)
[x,y] = select_points(I)
X=x(1:length(x)-1);
Y=y(1:length(y)-1);
dis(1)=sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2)
dis(2)=sqrt((X(3)-X(4))^2+(Y(3)-Y(4))^2)

%estimating momentum
for i=1:length(X)/2
    phi(i) = exp(-dis(i)/(sigma^2))
    K = [1,phi(i);phi(i),1];
    A = [K,zeros(length(X)/2);zeros(length(X)/2),K]
    %A=[1,phi(i),0,0;
       %phi(i),1,0,0;
       %0,0,1,phi(i);
       %0,0,phi(i),1]
    for j=1:1:2
        B(j)=[X(j+i*1)-X(j+(i-1)*1)]
        B(j+2)=[Y(j+i*1)-Y(j+(i-1)*1)]
    end
    param = pivot_gauss(A,B')
end

for i=1:1:size(newI,1)
        for j=1:1:size(newI,2)
            current_v=[i;j];
            dist = [j-X(1);i-Y(1)];
            d = sqrt(dist(1)^2+dist(2)^2);
            weight = exp(-d/(sigma^2))
            dist2 = [j-X(3);i-Y(3)];
            d2 = sqrt(dist2(1)^2+dist2(2)^2);
            weight2 = exp(-d2/(sigma^2))
            new_v=[i+weight*(param(1))+param(3)*weight2;j+weight*(param(2))+param(4)*weight2];
            if new_v(1) <=1
                new_v(1)=1;
            end
             if new_v(2) <=1
                new_v(2)=1;  
             end
             if new_v(1) > size(I2,1)
                new_v(1)=size(I2,1); 
             end
             if new_v(2) > size(I2,1)
                new_v(2)=size(I2,2);  
             end
             newI(i,j)=I2(round(new_v(1)),round(new_v(2)));

        end
end

figure(435)
imagesc(newI)
colormap(gray)
axis square
%B=[1,phi(i);1,phi(i)]
%param = pivot_gauss(X,Y)
