//*******************************************
// This is the Scilab script for Exercise 9.
//
// Use the help facility for more information 
// on individual functions used.
//
// Author: J. Kaempf, 2015 (updated)
//********************************************
clf; scf(0); a=gcf(); a.figure_size= [1500,500]; 

// read input data
eta=read("eta.dat",-1,201); eta0=read("eta0.dat",-1,201);
h0=read("h0.dat",-1,201);

x = (1:1:201)'; y = (1:1:51)'; // location vectors  
ntot = 100; // total number of frames

for n = 1:100 // animation loop

drawlater; clf;

// grab respective data block
jtop = (n-1)*51+1; jbot = jtop+50; 
etac = eta(jtop:jbot,1:201); 
etacc = etac-eta0;

// exclude unwanted data from plot
for j = 1:51; for k = 1:201;
  if h0(j,k) < 0; etacc(j,k) = %nan; end;
end; end;

plot3d(x,y,-0.6*h0',-70,50-5,'',[5,1,0]);
plot3d(x,y,15*etacc',-70,50-5,'',[4,5,3],ebox=[0,200,1,51,-20,20])

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
