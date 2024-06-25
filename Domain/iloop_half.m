function [fd,bbox,pfix]=iloop_half()
% mesh of a specimen: initialization 
fd=@dloop_half;
bbox=[-50,-70;0,70];
pfix=[0,-70;0,-30;0,30;0,70;-50,0;-30,0;-20,-70;-20,70;-10,-70;-10,70;-5,-70;-5,70;-15,-70;-15,70];