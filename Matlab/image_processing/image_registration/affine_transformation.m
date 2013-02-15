%%
% CS 6640 : Image Processing Project 2
%
% Author : Arthur COSTE
% Date : September 2012
%
% Content : Affine Transformation
%
%%
clear all
%close all
I=imread('CTscan.tif');
I2=double(I(:,:,1));
oldposition = zeros(size(I2,1),size(I2,2));
interpol=1
sx=2;
sy=1;
theta = 21;
tx = -20;
ty= 40;
transfo_mat = [sx*cosd(theta),sy*sind(theta), tx; -sx*sind(theta),sy*cosd(theta),ty;0,0,1];
newI = zeros(1*size(I2,1),1*size(I2,2));

%     for i=1:1:size(I2,1)
%         for j=1:1:size(I2,2)
%             current_v=[i;j;1];
%             current_intensity=I2(i,j);
%             new_v=transfo_mat*current_v;
%             newI(ceil(abs(new_v(1))+1),ceil(abs(new_v(2))+1)) = current_intensity;
%            
%         end
%     end
  
if interpol==0
    % Nearest Neighbor Interpolation
    
     for i=1:1:size(newI,1)
        for j=1:1:size(newI,2)
            current_v=[i;j;1];
            new_v=(transfo_mat)\current_v;
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
            new_v=(transfo_mat)\current_v;
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
            
            % Compute the Q_ij
            neighbor=[floor(new_v(1)), ceil(new_v(1)), floor(new_v(2)),ceil(new_v(2)) ]
            % Compute coefficients
            b1 = I2(neighbor(1),neighbor(3));
            b2 = I2(neighbor(2),neighbor(3));%-I2(neighbor(1),neighbor(3));
            b3 = I2(neighbor(1),neighbor(4));%-I2(neighbor(1),neighbor(3));
            b4 = I2(neighbor(2),neighbor(4));%-I2(neighbor(2),neighbor(3))-I2(neighbor(1),neighbor(4))+I2(neighbor(2),neighbor(4));
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


    figure(121)
    imagesc(I2)
    colormap(gray)
    title('Input Image')
    axis square
    figure(124)
    imagesc(newI)
    colormap(gray)
    axis square
    title('Transformed Image')
 end