# MECH105-Algorithms
These are the algorithms I created for MECH 105.
## Table of Contents
1. Homework 2. *Algorithm 1* for Solving a simple electrical circut.
2. Homework 2. *Algorithm 2* for Solving the degredation of Aqueous Bromide.
3. Homework 



### Algorithm 1
% Function parameters

q0 = 10; %
R = 60;
L = 9;
C = 0.00005;

% Use linspace to create an array of 100 points between 0 and 0.8

t = linspace(0,.8,100);

% Calculate the values of q

q = q0*exp((-R*t)/(2*L)).*cos((((1/(L*C))-(R/(2*L))^2)^.5)*t);

% Plot q vs t

subplot(1,2,1)
plot(t,q),title('Capacitor Charge'),xlabel('time'),ylabel('Charge')
% Make the capacitor 10x bigger
q2=q0*exp((-R*t)/(2*L)).*cos(sqrt((1/(L*(C*10)))-(R/(2*L))^2)*t);

% Plot q2 vs t

subplot (1,2,2)
plot(t,q2),title('Bigger Capacitor Charge'),xlabel('time'),ylabel('Charge')
### Algorithm 2
% Given experimental data
t_exp = 10:10:60;
c_exp = [3.4 2.6 1.6 1.3 1.0 0.5];

% Expected function
t_func = 0:.5:70
c_func = 4.84*(exp(-.034*(t_func)))

% Plot
hold on
plot(t_exp,c_exp,'rd')
plot(t_func,c_func,'g--')
xlabel('time in minutes')
ylabel('concentration in ppm')
legend('expected function','experimental data')
hold off
