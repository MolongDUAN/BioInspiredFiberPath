function fd=dloop_half(p)
% Signed distance function for a square with hole
r1=50; r2=20; r3=30; r4=5;
xc1=0; xc2=0;  xc3=0; xc4=0; xc5=0; xc6=0;
yc1=0; yc2=-50; yc3=50; yc4=0; yc5=-50; yc6=50;
dc1=dcircle(p,xc1,yc1,r1);
% dc2=dcircle(p,xc2,yc2,r2); dc3=dcircle(p,xc3,yc3,r2); 
dc4=dcircle(p,xc4,yc4,r3);
% dc5=dcircle(p,xc5,yc5,r4); dc6=dcircle(p,xc6,yc6,r4);
dr=drectangle(p,0,50,-70,70);
dr2=drectangle(p,-20,20,-70,-40);
dr3=drectangle(p,-20,20,40,70);

d1=min(dc1,min(dr2,dr3));
% d2=min(dr,min(dc4,min(dc5,dc6)));
d2=min(dr,dc4);
fd=max(d1,-d2);