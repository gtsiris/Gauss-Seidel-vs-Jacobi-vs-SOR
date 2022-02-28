%Georgios Tsiris, 1115201700173
clear; clc;
%disp('METHODOS SOR iterations');

prompt_a = 'Pick parameter a: ';
a=input(prompt_a)

prompt_b = 'Pick parameter b: ';
b=input(prompt_b)

prompt_n = 'Pick parameter n: ';
n=input(prompt_n)

ai(1:n-1)=-a;
bi(1:n)=4;
ci(1:n-1)=-b;
A=diag(ai,-1)+diag(bi)+diag(ci,1);

x=ones(n,1);
d=A*x;

x0=d;
x1=x0;

tol=(1/2)*10^-6;
maxits=100;

itcount=0;


CL=-tril(A, -1);
%disp('CL');disp(CL);

CU=-triu(A, 1);
%disp('CU');disp(CU);

I=eye(n);
%disp(I);
D=diag(diag(A));
%disp('D');disp(D);

D1=inv(D);
%disp('D1'); disp(D1);

L=D1*CL;
%disp('L');disp(L);

U=D1*CU;
%disp('U');disp(U);

B=L+U;
%disp('B');  disp(B);
%idiotimes tou pinaka B (Jacobi)
x=eig(B);     %disp(x);
rB=max(abs(x));
%m_low=min(abs(x));
%m_up=max(abs(x)); 
%disp('m_low'); disp(m_low); disp(' m_up'); disp( m_up);
disp(' rB'); disp(rB);
%omega=2.0/(1.0+sqrt(1-m_up*m_up));
omega=2.0/(1.0+sqrt(1-rB*rB));
disp('omega'); disp(omega);

t=tic;
while itcount<=maxits
    x0=x1;
    %disp('x0'); disp(x0);
    x1=inv(I-omega*U)*((1-omega)*I+omega*L)*x0+omega*inv(I-omega*U)*D1*d;
    %printf("step %d - x1:\n", itcount);disp(x1);
    nm=norm(x1-x0, Inf);
    %disp('nm'); disp(nm);
    if nm<tol
        iter=itcount;
        disp('siglisi se'); disp(iter); disp('epanalipseis');
        break;
    end
    itcount=itcount+1;
    %disp('itcount'); disp(itcount);
    %disp('x0'); disp(x0);
end% while
%if nm<tol
%disp('siglisi se'); disp(itcount); disp('epanalipseis');
%end
if nm>tol
    disp('oxi siglisi meta apo'); disp(maxits); disp('epanalipseis');
end
%disp('x1'); disp(x1);
e=toc(t);
%disp('e=cputime');disp(e);
