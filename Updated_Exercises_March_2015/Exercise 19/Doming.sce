//*******************************************
// This is the Scilab script for Exercise 19.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

clf; scf(0); a=gcf(); a.figure_size= [900,500];

// read input data
hin=read("h1.dat",-1,101); 

x1 = (0:10:1000)'; y1 = (0:10:500)'; // location vectors
[ntot nx] =size(hin); ntot = floor(ntot/51);

for n = 1:ntot // animation loop

time = real(60*n);

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
h = hin(jtop:jbot,1:101); 

drawlater; clf();

xset("thickness",2); xset("fpf"," "); //suppress label output 
rect = [0, 1000, 0, 500, -100, 100];
plot3d(x1,y1,(200-h'),-60,79,' ',[-4,1,4],rect); // surface plot
xset("thickness",1); 
contour(x1,y1,200-h'+6,15, -60,79,' ',flag=[0,2,4],rect); // contour plot

str = "time = "+string(int(100*time/24)/100)+" days";
xstring(400,950,str); //add time info
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(670,-400,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(1200,-100,"y (km)"); // add label for y-axis
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
// end

end; // end of animation loop
