//*******************************************
// This is the Scilab script for Exercise 18.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [900,450];

// read input data
tin=read("T.dat",-1,101); h0=read("h.dat",51,101);
etain=read("eta.dat",-1,101); 
u1=read("u.dat",-1,101); v1=read("v.dat",-1,101);

x = (0:10:1000)'; y = (0:10:500)'; // location vectors  
[ntot nx] =size(u1); ntot = floor(ntot/51);

for n = 1:ntot// animation loop

time = real(60*n);

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:101); v2 = v1(jtop:jbot,1:101);
t = tin(jtop:jbot,1:101); eta = etain(jtop:jbot,1:101);

drawlater;clf();

xset("thickness",2);xset("fpf"," "); //suppress label output 
// sea-level contours
contour2d(x,y,eta',10,1:10);
b = gca(); b.font_size = 3; b.data_bounds = [0,0;1000,500];
b.margins = [0.125,0.125,0.125,0.125]
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 200 400 600 800 1000], ["0" "200" "400" "600" "800" "1000"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 100 200 300 400 500], ["0" "100" "200" "300" "400" "500"]);

uu = u2(1:5:51,1:5:101);
vv = v2(1:5:51,1:5:101);

xset("thickness",1); champ(x(1:5:101),y(1:5:51),uu',vv',1.5); // vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;1000,500];

xstring(400,500.5,"time = "+string(int(100*time/24)/100)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 4;
xstring(490,-40,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-80,230,"y (km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow;

// save frames as GIF files (optional)
// if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
// else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
// end
end;
