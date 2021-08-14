//*******************************************
// Scilab script for visualisation of the
// Coriolis force.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (update)
//********************************************

clf(); scf(0); a=gcf(); a.figure_size= [600,600];

T = 24.*3600.; // period
fre = 2*%pi/T; // rotation rate
radius = 20; 
time = 0;

dt = T/200; // time step

xp = 6;
yp = 6;
up = 0.0; // speed in m/s
vp = 0.0;

// time loop starts here

factor = 0.9;
fac2 = factor*radius;

xpp(1:200)= %nan;
ypp(1:200)= %nan;

for n = 1:200 // run for 1 period

drawlater; clf;

// isoview scaling 
plot2d(0,0,-1,"030"," ",[-20,-20,20,20])

time = time+dt;

xset("color",2)
xfarc(-fac2,fac2,2*fac2,2*fac2,0,360*64)

xx = fac2*sin(-fre*time);
yy = fac2*cos(-fre*time);

x1 = [factor*xx xx];
y1 = [factor*yy yy];
xset("color",5)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

xx = fac2*cos(-fre*time);
yy = -fac2*sin(-fre*time);

x2 = [factor*xx xx];
y2 = [factor*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x2,y2,"lines",1)

xx = -fac2*sin(-fre*time);
yy = -fac2*cos(-fre*time);

x3 = [factor*xx xx];
y3 = [factor*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x3,y3,"lines",1)

xx = -fac2*cos(-fre*time);
yy = fac2*sin(-fre*time);

x4 = [factor*xx xx];
y4 = [factor*yy yy];
xset("color",0)
xset("thickness",6)
xpoly(x4,y4,"lines",1)

x1 = [0 0];
y1 = [factor*radius radius];
xset("color",5)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

x1 = [0 0];
y1 = [-radius -factor*radius];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

y1 = [0 0];
x1 = [factor*radius radius];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

y1 = [0 0];
x1 = [-radius -factor*radius];
xset("color",0)
xset("thickness",6)
xpoly(x1,y1,"lines",1)

//new ball position

alpha = cos(dt*2*fre);
beta = sin(dt*2*fre);

upn= alpha*(up-dt*2*fre*vp)-beta*dt*2*fre*up;
vpn= alpha*(vp+dt*2*fre*up)+beta*dt*2*fre*vp;

xp = xp+dt*upn/1000;
yp = yp+dt*vpn/1000;

xpn = xp*cos(-fre*time)+yp*sin(-fre*time);
ypn = yp*cos(-fre*time)-xp*sin(-fre*time);

xpp(n) = xpn;
ypp(n) = ypn;

up = upn;
vp = vpn;

xset("color",1)
xfarc(xpn-1,ypn+1,2.0,2,0,360*64)
xset("color",7)
xfarc(xpn-0.8,ypn+0.8,1.6,1.6,0,360*64)

xset('thickness',2)
plot2d(xpp,ypp,8,'000','',[-radius,-radius,radius,radius])

title("Fixed Frame of Reference","fontsize",3,"font_style",4);

xpause(5d2);

// creation of GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end
drawnow;

end;
