//*******************************************
// This is the Scilab script for Exercise 11.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************
clf; scf(0); a=gcf(); a.figure_size= [550,550];

// read input data
h0=read("h0.dat",-1,51); B1=read("b.dat",-1,51);

x = (0:0.1:5)'; y = (0:0.1:5)'; // location vectors  
[ntot nx] =size(B1);
ntot = floor(ntot/nx); time = 0.0;

for n = 1:3*ntot/8+1 // animation loop

time = time+0.5; 
jtop = (n-1)*51+1; jbot = jtop+50; 
// grab respective data block
B = B1(jtop:jbot,1:51); 

drawlater; clf;

xset("thickness",1); xset("fpf"," "); //suppress label output 
col(1:11) = 1; contour2d(x,y,h0',[0:2:20],col); // bathymetry contours
b = gca(); b.thickness = 3;

bmax = max(B);
contour2d(x,y,B',[0:bmax/12:bmax],[1:13],"011"); 
xset("fpf"," "); //suppress label output 

b = gca(); b.font_size = 3; b.data_bounds = [0,0;5,5];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);

xstring(2.8,5.05,"time = "+string(int(time))+" hrs"); //add time
c = gce(); c.clip_state = "off"; c.font_size = 4;

xstring(0.2,5.05,string(int(bmax*100))+" per cent"); //add time
c = gce(); c.clip_state = "off"; c.font_size = 3;

xstring(2.2,-0.6,"x (km)"); // add label for x-axis
c = gce(); c.clip_state = "off"; c.font_size = 3;
xstring(-0.5,2.6,"y"); // add label for y-axis
c = gce(); c.clip_state = "off"; c.font_size = 3;
xstring(-0.6,2.2,"(km)"); // add label for y-axis
c = gce(); c.clip_state = "off"; c.font_size = 3;
xpause(4d5); 

drawnow();

// save frames as GIF files (optional)
//if n < 10 then
//  xs2gif(0,'ex100'+string(n)+'.gif')
//else
//  if n < 100 then
//    xs2gif(0,'ex10'+string(n)+'.gif')
//  else
//    xs2gif(0,'ex1'+string(n)+'.gif')
//  end
// end

end; //end of animation loop

