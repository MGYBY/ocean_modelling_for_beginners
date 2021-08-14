//*******************************************
// This is the Scilab script for Exercise 12.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

// This code requires the bathymetry data from Exercise 10; that is, 
//  the file "h0.dat" needs to be available in the current data folder.

clf; scf(0); a=gcf(); a.figure_size= [550,550];

// read input data
h0=read("h0.dat",-1,51); 
trx1=read("TRx.dat",-1,5000); try1=read("TRy.dat",-1,5000); 

x = (0:0.1:5)'; y = (0:0.1:5)'; // location vectors  
[ntot ntra] =size(trx1); // extract dimension

drawlater; clf;

xset("thickness",2); xset("font size",3);
col(1:11)=1; contour2d(x,y,h0',[0:2:20],col); // bathymetry contours
xset("fpf"," "); //suppress label output 



ii = 0;

// I used only a limited number of floats (500 here)
for i = 1:500
  ii = ii+1;
  trax = trx1(1:ntot,i); // grab respective data line
  tray = try1(1:ntot,i); // grab respective data line
  xset("thickness",2); xset("font size",2);
  xset("mark size",0.2)
  if ii > 10; ii = ii-10; end;
  plot2d(trax,tray,ii,rect=[0,0,5,5],nax=[4,6,4,6]); // plot trajectory
end;

b = gca(); b.font_size = 3; b.data_bounds = [0,0;5,5];
b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3];
b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);
b.y_ticks = tlist(["ticks", "locations","labels"],..
  [0 1 2 3 4 5], ["0" "1" "2" "3" "4" "5"]);

xstring(2.2,-0.5,"x (km)"); // add label for x-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.5,2.6,"y"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;
xstring(-0.6,2.2,"(km)"); // add label for y-axis
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow();

