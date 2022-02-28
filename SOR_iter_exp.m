%Georgios Tsiris, 1115201700173
clear; clc;
%disp('METHODOS SOR iterations-Peiramatiko');

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

CU=-triu(A, 1);

I=eye(n);
D=diag(diag(A));

D1=inv(D);

L=D1*CL;

U=D1*CU;

B=L+U; 
%idiotimes tou pinaka B (Jacobi)
x=eig(B);
%disp('x'); disp(x);
rB=max(abs(x));  %m_up=rB
%disp(' rB'); disp( rB);

%ipologismos omega
disp('Methodos SOR');
omega=2.0/(1.0+sqrt(1-rB*rB));

disp('*******************    SOR     *********************');
fid611=fopen('C:\Users\maria\Desktop\demo\SOR_iter.txt','a+');
fprintf(fid611, 'w         iter    \n\n');

nm=0;
col=0;
for w=0.1:0.1:1.9
    nm=0;
    x0=d;
    x1=x0;
    tt=tic;
    itcount=0;
    while itcount<=maxits
        x0=x1;
        x1=inv(I-w*U)*((1-w)*I+w*L)*x0+w*inv(I-w*U)*D1*d;
        nm=norm(x1-x0, Inf);
        if nm<tol
            iter=itcount;
            %fprintf('siglisi se %d epanalipseis\n',iter);
            fprintf(fid611, '%f    %f    \n\n',w, iter);
            break;
        end
        itcount=itcount+1;
    end% while
    %apothikeusi apotelesmaton se pinakes gia dimiourgia
    %grafikis parastasis
    col=col+1;
    finw(col,1)=w;
    finit(col,1)= itcount;
    %disp(finw);  disp(finit);
    if nm>tol
        %disp('oxi siglisi meta apo'); disp(maxits); disp('epanalipseis');
    end
    e=toc(tt);
    %disp('e=cputime');disp(e);
end  %for w

% Dimiourgia grafikis parastasis
plot(finw, finit, 'r');
xlabel('\omega');              %  add axis labels and plot title
ylabel('iter');                %  add axis labels and plot title
title('Behavior of iterations with respect to \omega');

disp('Ypologismos elaxistou apo peiramatiko prosdiorismo');
[iterations, thesi]=min(finit);
wopt=finw(thesi,1); %disp('wopt'); disp(wopt);
% disp('iterations');disp(iterations);
%disp('thesi');disp(thesi);
%disp('ypologismos olon ton theseon me tin idia peiramatiki fasmatiki aktina');
thesi=find(finit==iterations);%disp('thesi');disp(thesi);
wopt=finw(thesi,1); %disp('wopt'); disp(wopt);


disp('Peiramatika apotelesmata');
disp('wopt'); disp(wopt);
disp('iterations');disp(iterations);
disp('thesi');disp(thesi);

disp('Beltistes parametroi');
disp('omega'); disp(omega);

fclose(fid611);   
