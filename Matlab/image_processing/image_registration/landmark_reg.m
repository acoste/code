%%
% CS 6640 : Image Processing Project 2
%
% Author : Arthur COSTE
% Date : October 2012
%
% Content : Calculation of affine transformation based on landmarks
%
%%
clear all
close all

I=imread('CTscan2.tif');
I2=double(I(:,:,1));
I3=imread('CTscan5.tif');
I4=double(I3(:,:,1));
newI = zeros(1*size(I2,1),1*size(I2,2));
interpol=1
gauss=0
figure(121)
[x,y] = select_points(I)
figure(122)
[x2,y2] = select_points(I3)

if gauss == 1
X=[x(1),y(1),1,0,0,0;
    0,0,0,x(1),y(1),1;
    x(2),y(2),1,0,0,0;
    0,0,0,x(2),y(2),1;
    x(3),y(3),1,0,0,0;
    0,0,0,x(3),y(3),1]

X2=[x2(1),y2(1),x2(2),y2(2),x2(3),y2(3)];
param = pivot_gauss(X,X2)

end

if gauss == 0
    X=[x(1),y(1),1,0,0,0;
    0,0,0,x(1),y(1),1];
    X2=[x2(1),y2(1)];
    
    for i=2:length(x)
        X=[X;x(i),y(i),1,zeros(1,3);zeros(1,3),x(i),y(i),1]
        X2 = [X2,x2(i),y2(i)]
    end
    param = X\X2'
end



sx = (x2(2)-x2(1))/(x(2)-x(1))
sy = (y2(2)-y2(1))/(y(2)-y(1))
transfo_matrix = [param(1),param(2), param(3);
                  param(4),param(5), param(6);
                  0,       0,       1       ]
% theta = asin(param(2))
% tx=param(3)
% ty=param(6)
%               
% transfo_matrix = [cosd(theta),sind(theta), tx; -sind(theta),cosd(theta),ty;0,0,1];


if interpol==0
    % Nearest Neighbor Interpolation
    
     for i=1:1:size(newI,1)
        for j=1:1:size(newI,2)
            current_v=[i;j;1];
            new_v=(transfo_matrix)\current_v;
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
            intensity =  I2(round(new_v(1)),round(new_v(2)));
            newI(i,j) = intensity;
        end
     end


    figure(121)
    imagesc(I2)
    colormap(gray)
    title('Input Image')
    axis square
    figure(122)
    imagesc(newI)
    colormap(gray)
    axis square
    title('Transformed Image')
end
 
if interpol==1
    % bilinear Interpolation
    
     for i=1:1:size(newI,1)
        for j=1:1:size(newI,2)
            current_v=[i;j;1];
            new_v=(transfo_matrix)\current_v;
            if new_v(1) <=1
                new_v(1)=1;
            end
             if new_v(2) <=1
                new_v(2)=1;  
             end
             if new_v(1) > size(I4,1)
                new_v(1)=size(I4,1); 
             end
             if new_v(2) > size(I4,1)
                new_v(2)=size(I4,2);  
             end
            
            % Compute the Q_ij
            neighbor=[floor(new_v(1)), ceil(new_v(1)), floor(new_v(2)),ceil(new_v(2)) ]
            % Compute coefficients
            b1 = I4(neighbor(1),neighbor(3));
            b2 = I4(neighbor(2),neighbor(3));%-I2(neighbor(1),neighbor(3));
            b3 = I4(neighbor(1),neighbor(4));%-I2(neighbor(1),neighbor(3));
            b4 = I4(neighbor(2),neighbor(4));%-I2(neighbor(2),neighbor(3))-I2(neighbor(1),neighbor(4))+I2(neighbor(2),neighbor(4));
%             b1 = I2(neighbor(1),neighbor(3))
%             b2 = I2(neighbor(2),neighbor(3))
%             b3 = I2(neighbor(1),neighbor(4))
%             b4 = I2(neighbor(2),neighbor(4))
            % compute new intensity
            %newint = b1+b2*(new_v(1))*(1-new_v(2))+b3*(new_v(2))*(1-new_v(1))+b4*(new_v(1)*new_v(2))
             newint = b1*(neighbor(2)-new_v(1))*(neighbor(4)-new_v(2))+b2*(new_v(1)-neighbor(1))*(neighbor(4)-new_v(2))+b3*(neighbor(2)-new_v(1))*(new_v(2)-neighbor(3))+b4*(new_v(1)-neighbor(1))*(new_v(2)-neighbor(3))
              
            intensity =  newint;
            newI(i,j) = intensity;
        end
     end

%     figure(121)
%     imagesc(I2)
%     colormap(gray)
%     title('Input Image')
%     axis square
    figure(124)
    imagesc(newI)
    colormap(gray)
    axis square
    title('Transformed Image')
%     hold on
%     plot(x(1:3), y(1:3), 'r-','linewidth',3);
%     plot(x(1:3), y(1:3), 'b+','linewidth',3);
 end