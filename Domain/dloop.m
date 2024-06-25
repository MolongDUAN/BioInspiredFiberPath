function fd=dloop(p)
% Signed distance function for a square with hole
r1=50; r2=30; r3=16.73408;
xc1=0; xc2=0; xc3=-26.69363216; xc4=26.69363216; xc5=-26.69363216; xc6=26.69363216;
yc1=0; yc2=0; yc3=61.16279497; yc4=61.16279497; yc5=-61.16279497; yc6=-61.16279497;
dc1=dcircle(p,xc1,yc1,r1);
dc2=dcircle(p,xc2,yc2,r2);
dc3=dcircle(p,xc3,yc3,r3);
dc4=dcircle(p,xc4,yc4,r3);
dc5=dcircle(p,xc5,yc5,r3);
dc6=dcircle(p,xc6,yc6,r3);

dr1=drectangle(p,-20,20,40,60);
dr2=drectangle(p,-20,20,-60,-40);
dr3=drectangle(p,-10,10,50,75);
dr4=drectangle(p,-10,10,-75,-50);

d1=min(dc1,min(dr4,min(dr3,min(dr1,dr2))));
d2=min(dc2,min(dc3,min(dc4,min(dc5,dc6))));
fd=max(d1,-d2);