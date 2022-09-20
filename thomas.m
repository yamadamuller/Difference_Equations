%Thomas' Algorithm - Solving Tridiagonal Matrixes
%Based on https://folk.ntnu.no/leifh/teaching/tkt4140/._main040.html

A = [1 0 0 0; 
    0.44444 -1.1610 0.63492 0; 
    0 0.44444 -1.0622 0.57778;
    0 0 0 1];
s = [0.008; 
    0.;
    0.;
    0.003];

a = [0; A(2,1); A(3,2); A(4,3)]; %lower diagonal
b = diag(A); %main diagonal
c = [A(1,2); A(2,3); A(3,4); 0]; %lower diagonal
d = [s(1,1);s(2,1);s(3,1);s(4,1)]; %right hand side

N = length(b); %necessary iterations
x = zeros(N,1); %ui matrix

%elimination
for i=2:N
    q = a(i)/b(i-1);%qj=aj/bj-1
    b(i)=b(i)-c(i-1)*q;%bj=bj-qj*cj-1
    d(i)=d(i)-d(i-1)*q;%dj=dj-qj*dj-1
end

%substitution
x(N)=d(N)/b(N);%xN=dN/bN

for j=N-1:-1:N-3 %(N-1,N-2,N-3)
    x(j)=(d(j)-c(j)*x(j+1))/b(j); %xj = dj-cj*xj+1/bj
end

disp(x);








