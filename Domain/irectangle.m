function [fd,bbox,pfix]=irectangle()
% mesh of a specimen: initialization 
fd=@drect;
bbox=[0,0;10,15];
pfix=[0,0;2,0;4,0;6,0;8,0;10,0;
      0,15;2,15;4,15;6,15;8,15;10,15];