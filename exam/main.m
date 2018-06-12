clear all;
close all;
clc;
A=[1 3 4 5;1 2 3 4;0 0 1 2;3 3 2 4];
[u,b,v]=gkl_bidiag(A);

b0=u'*A*v;
b0-b

