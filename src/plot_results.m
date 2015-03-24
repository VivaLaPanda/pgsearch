clear all;
close all;
clc;


MeanTest = csvread('MeanCost.txt');
PGTest = csvread('PGCost.txt');
x1 = mean(MeanTest,2);
x2 = mean(PGTest,2);
plot(x1, 'r');
hold on
plot(x2);
legend('Mean Search', 'PG Search');