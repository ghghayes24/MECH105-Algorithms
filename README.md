# MECH105-Algorithms
These are the algorithms I created for MECH 105. Remember to add appropriate multiplication signs as they had to be removed for coherance.
## Table of Contents
1. _Homework 2_. *Algorithm 1* for Solving a simple electrical circut.
2. _Homework 2_. *Algorithm 2* for Solving the degredation of Aqueous Bromide.
3. _Homework 3_. *Algorithm 3* for solving the volume of a tank.
4. _Homework 4_. *Algorithm 4* for creating special matrix m x n.
5. 



### Algorithm 1
% Function parameters

q0 = 10; %
R = 60;
L = 9;
C = 0.00005;

% Use linspace to create an array of 100 points between 0 and 0.8

t = linspace(0,.8,100);

% Calculate the values of q

q = q0exp((-Rt)/(2L)).cos((((1/(LC))-(R/(2L))^2)^.5)t);

% Plot q vs t

subplot(1,2,1)
plot(t,q),title('Capacitor Charge'),xlabel('time'),ylabel('Charge')
% Make the capacitor 10x bigger
q2=q0exp((-Rt)/(2L)).cos(sqrt((1/(L(C10)))-(R/(2L))^2)t);

% Plot q2 vs t

subplot (1,2,2)
plot(t,q2),title('Bigger Capacitor Charge'),xlabel('time'),ylabel('Charge')
### Algorithm 2
% Given experimental data

t_exp = 10:10:60;
c_exp = [3.4 2.6 1.6 1.3 1.0 0.5];

% Expected function

t_func = 0:.5:70
c_func = 4.84(exp(-.034(t_func)))

% Plot

hold on
plot(t_exp,c_exp,'rd')
plot(t_func,c_func,'g--')
xlabel('time in minutes')
ylabel('concentration in ppm')
legend('expected function','experimental data')
hold off
### Algorithm 3
% Specify the variables needed to solve this problem (ie. height of each section, diameter, radiaus, ...)
%   It is alwasy easier to work with variables (diameter_cyl = 25) than to use numbers everywhere, since a 
%   diameter indicates something specific but the number 25 could mean anything

dc=25;
hc=19;
dt=46;
rc=.5dc;
rt=.5dt;
% Specify the height of the water

h = 20
% You can comment / uncomment lines below for testing. This will overwrite the previous line for h = 20.
% For submission, make sure all of the following lines are commented out and h = 20! (OR IT IS MARKED AS WRONG)

%h = 5
%h = 19
%h = 47
%h = -1

rn=((h-19)/(14/(rt-rc))+rc);

% Now compute the volume. Using conditional statments you will want to first check the height makes sense,
% and then solve the volume depending on what portion of the tank has been filled.
% Make sure that your volume is stored in the variable v! (OR IT WILL BE MARKED AS WRONG)
% You may find it more convenient to move v around in you code, it is only given here to indicate what variable to use

if h<=19 && h>0
    v=(pirc^2h)
elseif 33>=h && h>19
    v=(pirc^2hc)+((pi(h-19)(rc^2+rn^2+(rcrn)))/3)
else h>33||h<0
    disp('error! Invalid Height Entered')
end
fprintf('v=%f',v)
### Algorithm 4
function [A] = specialMatrix(n,m)
% This function should return a matrix A as described in the problem statement
% Inputs n is the number of rows, and m the number of columns
% It is recomended to first create the matrxix A of the correct size, filling it with zeros to start with is not a bad choice

if nargin~=2
    error('Error. Invalid input')
end
A=zeros(n,m);
for i=1:n
    for j=1:m
       if i==1
           A(i,j)=j;
       elseif j==1
           A(i,j)=i;
       else
           A(i,j)=A(i-1,j)+A(i,j-1)
       end
    end
end
% Now the real challenge is to fill in the correct values of 

end
% Things beyond here are outside of your function
