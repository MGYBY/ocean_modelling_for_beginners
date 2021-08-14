//*******************************************
// This is the Scilab script for Exercise 16.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************
clf; scf(0); a=gcf(); a.figure_size= [1000,350];

// read input data
tin=read("t.dat",-1,151); h0=read("h0.dat",-1,151); etain=read("eta.dat",-1,151); 
u1=read("u.dat",-1,151); v1=read("v.dat",-1,151);
x = (0:1:150)'; y = (0:1:50)'; // location vectors   
[ntot nx] =size(u1); ntot = floor(ntot/51);

for n = 1:ntot// animation loop

time = real(6*n); 

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
u2 = u1(jtop:jbot,1:151); v2 = v1(jtop:jbot,1:151); 
t = tin(jtop:jbot,1:151); eta = etain(jtop:jbot,1:151);

drawlater; clf();

xset("fpf"," "); //suppress label output 
col1(1:10) = 5; contour2d(x,y,h0',10,col1); // bathymetry contours

xset("fpf"," "); //suppress label output 
col2 = 1:20; contour2d(x,y,t',[0.0:0.05:1.0],col2);// tracer contours
b = gca(); b.font_size = 3; b.data_bounds = [0,0;150,50];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [2,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 30 60 90 120 150], ["0" "30" "60" "90" "120" "150"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 10 20 30 40 50], ["0" "10" "20" "30" "40" "50"]);

uu = u2(1:5:51,1:5:151); vv = v2(1:5:51,1:5:151);
xset("thickness",1); champ(x(1:5:151),y(1:5:51),uu',vv',2.0); // velocity vector plot
b = gca(); b.font_size = 3; b.data_bounds = [0,0;150,50];

xstring(110,50.5,"time = "+string(int(100*time/24)/100)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 4;
xstring(57,-9.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-15,21,"y (km)"); // add label for y-axis
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

end; // end of iteration loop

