clc
clear all
close all
warning('off')

load('Dataset\Data');

dict.dir='Dictionaries\';
dict.learn=0; %set to 1 to learn dictionary
dict.iter=25;
dict.S=2;
dict.errorFlag=0;
dict.errorGoal=1;
dict.totalerr=2;
dict.nn=1;
method.S=2;
kbest.sur=[3,3]; %surrounding points about each sensor
kbest.option=0; %set 1 to show pressure vs error 

subi=1;
legi=1;
[~,y,~,~,W]=testallfunc(data,subi,legi,dict,method,kbest);