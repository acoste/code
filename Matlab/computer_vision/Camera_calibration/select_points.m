%%
% CS 6640 : Image Processing Project 3
%
% Author : Arthur COSTE
% Date : October 2012
%
% Content : select landmarks
%% 
function [x,y]=select_points(I)
  %function [x,y]=select_points(I)
% 
% I: image
  % [x,y]: points coordinates
%
% Select points on an image
% left click: select position
% right click: last point
% 
%

colormap(gray);
x=[];
y=[];

disp('left click: select position');
disp('right click: last point');


but = 1;
while but == 1
  clf
  imagesc(I);
  colormap(gray)
  hold on
%   plot(x, y, 'r-','linewidth',3);
  plot(x, y, 'b+','linewidth',3);
  axis square
  
  [s, t, but] = ginput(1);

x = [x;s];
y = [y;t];

end  

