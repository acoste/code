%%
% CS 6640 : Image Processing Project 4
%
% Author : Arthur COSTE
% Date : November 2012
%
% Content : Hough Transform 2 : using edge orientation
%
%%
function OutputIm = Hough_transform_2(filtering)
gradient=imread('piste_edges.tif');
%gradient=imread('edges-lines-orig_edges.tif');
%gradient=imread('slc2.tif');

I2=double(gradient(:,:,1));
%angle = imread('edges-lines-orig_angle.tif');
%angle = imread('slc_angle.tif');
angle = imread('angle_map2.tif');
I3=double(angle(:,:,1));
newI = zeros(600,362);
orig = [1,0];
figure(1)
imagesc(I2)
colormap(gray)
axis square
k=1;
tic
for i=1:1:size(I2,1)
    for j=1:1:size(I2,2)
        if I2(i,j)==255
            Y(k)=i;
            X(k)=j; 
            theta_im(k) = I3(i,j)*(360/255);
            theta_lower = round((theta_im(k))) - 8;
            if (theta_lower <=0) theta_lower = 1; end
            theta_upper = round((theta_im(k))) + 8;
            if (theta_upper >360) theta_upper = 360; end
            for theta=theta_lower:1:theta_upper
                rho = (X(k))*cosd(theta) + (Y(k))*sind(theta);
                if rho <=orig(1)
                    rho=orig(1)-rho+1;
                end
                 if rho > size(newI,1)
                    rho=size(newI,1); 
                 end
                 
            newI(round(rho),round(theta)) = newI(round(rho),round(theta)) + 1;
        
            end
            k=k+1;
        end
    
    end
end

if (filtering ==1)
weight=[1/16,2/16,1/16;2/16,4/16,2/16;1/16,2/16,1/16]
 for i = ceil(size(weight,1)/2) :1: size(newI,1)-size(weight,1)+ceil(size(weight,1)/2)
    for j = ceil(size(weight,2)/2) :1: size(newI,2)-size(weight,2)+ceil(size(weight,2)/2)
        convol=0;
        %compute convolution for the neighbourhood associated to the kernel
        for a = 1:size(weight,1)
            for b=1:size(weight,2) 
   
            convol = convol + (weight(a,b)*newI(i-a+ceil(size(weight,1)/2),j-b+ceil(size(weight,2)/2)));
            
            end
        end
        newI(i,j)=convol;
    end
end
end

temp=newI;
 figure(2326)
 imagesc(log(newI))
 title('Hough Transform')
 axis square
  maximum=max(max(newI));
seuil = maximum-2;
k=1;
 maxx=max(max(newI))
 for i=1:1:250
     for j=1:1:362
         if (newI(i,j)) >= maxx-26
             maxi_rho(k)=i;
             if j<180
             maxi_theta(k)=j;
             end
             if j>180
                 maxi_theta(k)=j+180;
             end
             k=k+1;
         end
     end
 end
toc
 figure(1)
 hold on
  for g=1:1:length(maxi_rho)
    x=1:1:450;
    y= (-(cosd(maxi_theta(g)))/sind((maxi_theta(g))))*x  + (maxi_rho(g))/sind((maxi_theta(g)));
    plot(x,y,'-r','LineWidth',2)
    hold on

 end