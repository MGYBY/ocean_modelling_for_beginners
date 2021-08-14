//*******************************************
// This is the Scilab script for Exercise 20.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (update)
//********************************************

// This script produces an animation of the evolution of the shape 
// of the density interface.

clf; scf(0); a=gcf(); a.figure_size= [500,400];

// read input data
hin=read("h1.dat",-1,51); 

x1 = (0:1:50)'; y1 = (0:1:50)'; // location vector;
[ntot nx] =size(hin); ntot = floor(ntot/51);

for n = 1:ntot // animation loop

time = real(2*n)-2; 
// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
h = hin(jtop:jbot,1:51); 

drawlater; clf();

rect = [0, 50, 0, 50, -100, 0];
plot3d(x1,y1,-h',-50,128,' ',[4,1,5],rect); // bathymetry contours

xstring(51,10,"time = "+string(int(10*time/24)/10)+" days"); //add time
b = gce(); b.clip_state = "off"; b.font_size = 3;

drawnow;

// create GIF files (optional)
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
