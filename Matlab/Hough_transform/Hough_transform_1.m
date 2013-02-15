%%
% CS 6640 : Image Processing Project 4
%
% Author : Arthur COSTE
% Date : November 2012
%
% Content : Hough Transform part 1
%
%%
function OutputIm = Hough_transform_1(decrement,filtering)
I=imread('edges-lines-orig_edges.tif');
I2=double(I(:,:,1));
newI = zeros(600,362);
orig = [1,0];
figure(1)
imagesc(I2);
colormap(gray)
axis square
k=1;
tic
%perform Hough Transform
for i=1:1:size(I2,1)
    for j=1:1:size(I2,2)
        if I2(i,j)==255
            Y(k)=i;
            X(k)=j; 
            for theta=1:1:360
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
temp=newI;
 figure(2326)
 imagesc((newI));
 title('Hough Transform')
 axis square
  maximum=max(max(newI));
  
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

if(decrement==1)
    for i=1:1:size(X,2)   
        for theta=1:1:360
            rho = X(i)*cosd(theta) + Y(i)*sind(theta);
                if rho <=orig(1)
                    rho=orig(1)-rho+1;
                end
                 if rho > size(newI,1)
                    rho=size(newI,1); 
                 end
 
            if  (newI(ceil(rho),round(theta))<=maximum-35)        
                newI(ceil(rho),round(theta)) = newI(round(rho),round(theta)) - 1;
                if newI(ceil(rho),round(theta)) <= 0
                    newI(ceil(rho),round(theta)) = 0;
                end
            end
        end
    end
end

% compute maximum on the accumulator and threshold
 k=1;
 maxx=max(max(newI))
 if (decrement==0)
    seuil = maxx-55;
end
if (decrement == 1)
    seuil = maxx-30;
end
if (filtering == 1)
    seuil = maxx-20;
end
if (filtering == 1 && decrement == 1)
    seuil = maxx-15;
end
 for i=1:1:250
     for j=1:1:360
         if (newI(i,j)) >= seuil
             maxi_rho(k)=i;
             if j<180
             maxi_theta(k)=j;
             end
             if j>=180
                 maxi_theta(k)=j-1;
             end
             k=k+1;
         end
     end
 end
 q=toc
 %display
 
 figure(2325)
 imagesc((newI));
 title('Hough Transform')
 axis square
 figure(1)
 hold on
 
  for g=1:1:length(maxi_rho)
    x=0:1:250;
    y= (-(cosd(maxi_theta(g)))/sind((maxi_theta(g))))*x  + (maxi_rho(g))/sind((maxi_theta(g)));
    plot(x,y,'-r','LineWidth',2)
    hold on

  end
 OutputIm=newI;