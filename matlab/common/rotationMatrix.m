function [R]=rotationMatrix(u, theta)
%Rodrigues' rotation formula :http://en.wikipedia.org/wiki/Rotation_matrix

P = u*u'; 

I =[ 1 0 0; 0 1 0; 0 0 1];

Q =[ 0    -u(3)   u(2)
    u(3)   0     -u(1)
   -u(2)  u(1)    0   ]; 

R = P + (I-P) * cos(theta) + Q * sin(theta);
