//*******************************************
// This is the Scilab script for Exercise 8.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************
clf; scf(0); a=gcf(); a.figure_size= [1000,700];

// read input data
eta=read("eta.dat",-1,51); eta0=read("eta0.dat",-1,51);

x = (1:1:51)'; y = (1:1:51)'; // location vectors  
ntot = 150; // total number of frames

for n = 1:ntot // animation loop

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
etac = eta(jtop:jbot,1:51); 

drawlater; clf;
plot3d(x,y,1000*(etac-eta0),alpha=50,theta=50,flag = [2,2,3]);
b = gca(); b.data_bounds = [1,1,-20;51,51,20];
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

