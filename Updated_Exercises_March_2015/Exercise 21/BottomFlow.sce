//*******************************************
// This is the Scilab script for Exercise 21.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 update
//********************************************

// This script produces an animation of currents
// and tracer distributions in the bottom layer.

clf; scf(0); a=gcf(); a.figure_size= [1000,450]; a.color_map = hsvcolormap(64);

// read input data
tin=read("t2.dat",-1,101); h0=read("h0.dat",-1,101); 
u1=read("u2.dat",-1,101); v1=read("v1.dat",-1,101);
e1=read("eta2.dat",-1,101);

x1 = (0:2:200)'; y1 = (0:2:100)'; // location vectors  
[ntot nx] =size(u1); ntot = floor(ntot/51);

for n = 1:ntot // animation loop

time = real(6*n)-6; 

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:101); v2 = v1(jtop:jbot,1:101);
t = tin(jtop:jbot,1:101); eta = e1(jtop:jbot,1:101);

// interpolate flow field onto coarser grid
u(1:26,1:51) = 0.; v(1:26,1:51) = 0.; // scaling 

for j = 1:25; for k = 1:50;
  j1 = 2*j-1; j2 = j1+1; k1 = 2*k-1; k2 = k1+1;
  uu = 0.; vv = 0.;
  for jstar = j1:j2; for kstar = k1:k2;
    uu = uu + u2(jstar,kstar); vv = vv + v2(jstar,kstar); 
  end; end;
  u(j,k) = uu/4.; v(j,k) = vv/4.;
end; end;

x = (2:4:202)'; y = (2:4:102)'; // location vectors  

drawlater; clf;

xset("fpf"," "); //suppress label output 
col = 1:40; contour2d(x1,y1,t',40,col); // tracer contours
champ(x,y,u',v',1.0);// vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;200,100];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 50 100 150 200], ["0" "50" "100" "150" "200"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 50 100], ["0" "50" "100"]);

xstring(80,100,"time = "+string(int(100*time/24)/100)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(115,-9.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-14,55,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow;

// save frames as GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
//end
 
end; // end of animation loop

