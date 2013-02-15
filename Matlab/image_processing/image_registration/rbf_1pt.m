%%
% CS 6640 : Image Processing Project 2
%
% Author : Arthur COSTE
% Date : October 2012
%
% Content : Non linear warping using Radial basis function
%
%%
%clear all
close all
color =['g','r'];
sigma = 0.010;
forward =0
I=imread('grid_2.tif');
I2=double(I(:,:,1));
%newI = zeros(1*size(I2,1),1*size(I2,2));
interpol=1
figure(122)
[x,y] = select_points(I)
X=x(1:length(x)-1);
Y=y(1:length(y)-1);

%estimating momentum 

    
% applying rbf transformation

if forward == 1
    alpha = [X(2)-X(1);Y(2)-Y(1)] % compute momentum Y0-X0

    for i=1:1:size(I2,2)
        for j=1:1:size(I2,1)
            current_v=[i;j;1];
            dist = [j-X(1);i-Y(1)];         %distance to X0
            d = sqrt(dist(1)^2+dist(2)^2);
            %weight = exp(-d/(sigma^2))
            weight = 1/(1+(sigma*d)^2)
            %weight = 1/sqrt(1+(sigma*d)^2)
            new_v=[i+weight*(alpha(1)),j+weight*(alpha(2))];
            newI(ceil(abs(new_v(1))),ceil(abs(new_v(2)))) = I2(i,j);
           
        end
    end
    figure(43623)
    imagesc(newI)
    colormap(gray)
    axis square
    title('RBF Forward Method')
end

if forward == 0
    for i=1:1:size(newI,2)
            for j=1:1:size(newI,1)
                alpha = [Y(2)-Y(1);X(2)-X(1)];
                %phi = [X(1)-i,Y(1)-j]
                dist = [j-X(1);i-Y(1)];
                d = sqrt(dist(1)^2+dist(2)^2);
                %weight = exp(-d/(sigma^2))
                %weight = (sigma^2)/(d^2+sigma^2)
                weight = 1/(1+(sigma*d)^2)
                %weight = 1/((d^2)*log(d))
                new_v=[i+weight*(alpha(1)),j+weight*(alpha(2))];

                 if new_v(1) <=1
                    new_v(1)=1;
                end
                 if new_v(2) <=1
                    new_v(2)=1;  
                 end
                 if new_v(1) > size(I2,1)
                    new_v(1)=size(I2,1); 
                 end
                 if new_v(2) > size(I2,2)
                    new_v(2)=size(I2,2);  
                 end
                newI(i,j)=I2(ceil(new_v(2)),ceil(new_v(1)));

            end
    end

    figure(43634)
    imagesc(newI)
    colormap(gray)
    hold on
    plot(x(1),y(1),'g+','linewidth',3);
    plot(x(2),y(2),'r+','linewidth',3);
    plot(x(1:2), y(1:2), color(1:2),'linewidth',3);
    axis square
    title('RBF Reverse Method')
end
