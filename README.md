# MECH105-Algorithms
These are the algorithms I created for MECH 105. Remember to add appropriate multiplication signs as they had to be removed for coherance.
## Table of Contents
1. _Homework 2_. *Algorithm 1* for Solving a simple electrical circut.
2. _Homework 2_. *Algorithm 2* for Solving the degredation of Aqueous Bromide.
3. _Homework 3_. *Algorithm 3* for solving the volume of a tank.
4. _Homework 5_. *Algorithm 4* for creating special matrix m x n.
5. _Homework 6.5_. *Algorithm 5* for converting binary.
6. _Homework 11_. *Algorithm 6* for finding roots using false position.
7. _Homework 11_. *Algorithm 7* for finding roots using the Bisect matlab func.
8. _Homework 17_. *Algorithm 8* for LU factorization with partial pivoting (only partially correct)
9. _Homework 22_. *Algorithm 9* for the simpsons 1/3rd rule.



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
### Algorithm 5
function [base2] = binaryConverter(base10)
%binary A simple function to convert a base10 number to binary
%   Inputs:
%       base10 - A number in base10
%   Outputs:
%       base2 - The base10 number converted to binary

nb=2;
sum=base10;
remainder=0;
i=1;
base2=0;
while sum>0
    remainder=mod(sum,nb);
    sum=floor(sum/nb);
    base2(i)=remainder
    i=i+1;
end
base2=flip(base2)
end
### Algorithm 6
function [root, fx, ea, iter] = falsePosition(func, xl, xu, es, maxit, varargin)
%falsePosition finds the root of a function using false position method

if nargin<3
    error('Need at least 2 unputs')
end
test=func(xl,varargin{:})func(xu,varargin{:});
if test>0
    error('No change in sign')
end
if nargin<4 | isempty(es)
    es=.0001;
end
if nargin<5 | isempty(maxit)
    maxit=200;
end
iter=0;
xr=xl;
ea=100;
while (1)
    xrold=xr;
    xr=xu-(func(xu)(xl-xu))/(func(xl)-func(xu));
    iter=iter+1;
    if xr~=0
        ea=abs(((xr-xrold)/xr)100);
    end
    test=func(xl,varargin{:})func(xr,varargin{:});
    if test<0
        xu=xr;
    elseif test>0
        xl=xr;
    else
        ea=0;
    end
    if ea<=es | iter>maxit
        break
    end
end
root=xr;
fx=func(xr,varargin{:});
end
### Algorithm 7
% Define problem constants

format long
g = 9.81;
mu = 0.55;
F = 150;
m = 25;
es=.001;
maxit=200;
d=mugm
xl=0;
xu=90;
theta=[xl:5:xu];
x=theta;
% Define the function to be solved for. Is the angle specified in radians or degrees?

func=@(x)[((d./(cosd(x)+(mu.sind(x))))-F)]; % x is theta
% THINK, what makes range sense for angle?


% Plot your function. Does plotting give you an idea about where the root is?
% Finaly solve for the root using the bisection script given to you, which can be called as:

plot(x,func(x))
xlabel('theta')
ylabel('func')
upper_bound=xu
lower_bound=xl
[root, fx, ea, iter] = bisect(func, lower_bound, upper_bound);
angle=root
disp(func)
% These variables will be checked for your final answer:
%angle =  % the angle in degree that solves this problem
%fx =     % the function value at your solved angle
%ea =     % the bisection error estimate
%iter =   % the number of bisection iterations

fprintf('%d\n,angle')
fprintf('%d\n,fx')
fprintf('%d\n,ea')
fprintf('%d\n,iter')
### Algorithm 8
function [L, U, P] = luFactor(A)
if nargin<1
    error('Need more input arg.s');
end
[m,n]=size(A);
if m~=n
    error('Needs to be a matrix where rows=collumns');
else
    for row=1:m
        for collumn=1:n
            if A(row,collumn)==0
                error('Matrix cannot have "0" terms')
            end
        end
    end
end
% luFactor(A)

L=zeros(n,m);
U=A;
P=eye(m);
%	LU decomposition with pivoting
% inputs:
%	A = coefficient matrix
% outputs:
%	L = lower triangular matrix
%	U = upper triangular matrix
%       P = the permutation matrix

for collumn=1:n
    [pivot row]=max(abs(U(collumn:n,collumn)));
    row=row+collumn-1;
        if row~=collumn
            U1=U(collumn,:);
            U(collumn,:)=U(row,:);
            U(row,:)=U1;
            p=P(collumn,:);
            P(collumn,:)=P(row,:);
            P(row,:)=p;
              if k>=2
                L1=L(collumn,1:collumn-1);
                L(collumn,1:collumn-1)=L(row,1:collumn-1);
                L(row,1:collumn-1)=L1;
              end
        end
end
for collumn=1:n
    for s=collumn+1:row
          L(s,collumn)=U(s,collumn)/U(collumn,collumn);
          U(s,collumn:row)=U(s,collumn:row)-L(s,collumn)U(collumn,collumn:row);
          L(m/n,n./m)=1;
          L(m,n)=1;
          L(n-1,m-1)=1;
    end
end
end
### Algorithm 9
function [I] = Simpson(x, y)
% Numerical evaluation of integral by Simpson's 1/3 Rule
% Inputs
%   x = the vector of equally spaced independent variable
%   y = the vector of function values with respect to x
% Outputs:
%   I = the numerical integral calculated

[A,B]=size(x);
k=length(x);
if k~=length(y)
    error('Impossible. Check lengths of arrays x and y');
else test=linspace(x(1),x(k),k);
    for i=1:k
        if test(i)~=x(i)
            error('Length of x must be consistent with Simpsons rule')
        end
    end
end
h=(x(2)-x(1))
if k==3
    I=(h/(3))(y(1)+4(sum(y(2:2:k-2)))+2(sum(y(3:2:k-1))+y(k)));
elseif mod(k,2)==0
    warning('Trapezoidal rule must be used as there is an even length of x')
    if k==2
        I=(x(k)-x(1))((y(k)+y(1))/2);
    else
        a=(h/(3))(y(1)+4(sum(y(2:2:k-2)))+2(sum(y(3:2:k-1))+y(k)));
        b=(x(k)-x(k-1))((y(k)+y(k-1))/2);
        I=a+b;
    end
else
    a=(h/(3))(y(1)+4(sum(y(2:2:k-2)))+2(sum(y(3:2:k-1))+y(k)));
    b=(x(k)-x(k-1))((y(k)+y(k-1))/2);
    I=a+b;
end
end
