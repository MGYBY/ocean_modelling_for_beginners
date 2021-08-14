//*******************************************
// This is the Scilab script for Exercise 12.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

// This code requires the bathymetry data from Exercise 10; that is, 
// the file "h0.dat" needs to be available in the current data folder.

clf; scf(0); a=gcf(); a.figure_size= [550,550];

// read input data
h0=read("h0.dat",-1,51); 
trx1=read("TRx.dat",-1,5000); try1=read("TRy.dat",-1,5000);

x = (0:0.1:5)'; y = (0:0.1:5)'; // location vectors  
[ntot ntra] =size(trx1);
time = 0.0;

for n = 1:ntot // animation loop

time = n*5/60; 
 
trax = trx1(n,1:ntra); // grab respective data line
tray = try1(n,1:ntra); // grab respective data line

drawlater; clf;

xset("thickness",2); xset("font size",3); xset("fpf"," "); //suppress label output 
contour2d(x,y,h0',[0:2:20],1:11); // bathymetry contours

//xset("thickness",2); xset("font size",2);

xset("mark size",0.2)
plot2d(trax,tray,-3,rect=[0,0,5,5],nax=[4,6,4,6]); 
b = gca(); b.font_size = 3; b.data_bounds = [0,0;5,5];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);

xstring(2.6,5.05,"time = "+string(int(10*time)/10)+" hrs"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 4;
xstring(2.2,-0.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.5,2.6,"y"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.6,2.2,"(km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

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

end; // end of animation loop

