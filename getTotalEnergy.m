function [ totalenergy totalEnergy_dB] = getTotalEnergy( input_dB )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tmpLinearinput = 10^(input_dB/20);
F0Energy = tmpLinearHarmonic.^2;                 
totalEnergy = sum(F0Energy);                          % total energy of the harmonics of this F0
totalEnergy_dB = 10*log10(totalF0Energy);             % total energy of the harmonics of this F0 (dB)
        
