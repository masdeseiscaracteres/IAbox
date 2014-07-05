clc;
clear;
close all;

K=4;
nT=randi([5 15],1,K);
nR=randi([5 15],1,K);
d=randi([5 15],1,K);
nT=[2 2 3 2];
nR=[3 3 3 3];
d=ones(1,K);

aux1=[nT' nR' d']

system_vec2str(nT,nR,d,'KeepOrder')
system_vec2str(nT,nR,d,'aaa')