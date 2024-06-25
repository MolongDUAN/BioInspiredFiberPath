function fd=dtype1(p)
% Signed distance function for a square with hole
r1=76; r2=76; r3=76; r4=76;
xc1=-82.5; xc2=82.5;  xc3=-82.5; xc4=82.5;
yc1=-28.5; yc2=-28.5; yc3=28.5; yc4=28.5;
dc1=dcircle(p,xc1,yc1,r1);
dc2=dcircle(p,xc2,yc2,r2);
dc3=dcircle(p,xc3,yc3,r3); 
dc4=dcircle(p,xc4,yc4,r4);

dr=drectangle(p,-9.5,9.5,-82.5,82.5);
dr2=drectangle(p,-9.5,-6.5,-28.5,28.5);
dr3=drectangle(p,6.5,9.5,-28.5,28.5);

d2=min(dr2,min(dr3,min(dc1,min(dc2,min(dc3,dc4)))));
fd=max(dr,-d2);