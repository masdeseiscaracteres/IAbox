clc;
clear;
close all;

[nT nR D options]=GetTestScenario([]);

H = GenerateChannel(nT,nR,options);

[U V]=RandomBeamforming(H,D,options)