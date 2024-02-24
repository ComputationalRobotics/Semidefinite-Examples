clc; clear; close all;

syms x y;
I = [x*y^3 - x^2, x^3*y^2-y];
G1 = gbasis(I,'MonomialOrder','degreeLexicographic');

G2 = gbasis(I,'MonomialOrder','lexicographic');