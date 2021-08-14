//*******************************************
// This is the Scilab script for Exercise 21.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************

// This script produces an animation of the evolution
// of interface displacements.

clf; scf(0); a=gcf(); a.figure_size= [800,400]; clf(gcf(), "reset"); // reset colormap

// read input data
hin=read("h1.dat",-1,101); 

x1 = (0:1:100)'; y1 = (0:1:50)'; // location vector;
[ntot nx] =size(hin); ntot = floor(ntot/51);

for n = 1:ntot // iteration loop

time = real(6*n)-6; 

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
h = hin(jtop:jbot,1:101); 

drawlater; clf();

rect = [0, 100, 0, 50, -300, 0];
plot3d(x1,y1,-h',-40,201,'',[2,1,5],rect); // interface displacements

xstring(120,-95,"time = "+string(int(100*time/24)/100)+" days"); //add time
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
