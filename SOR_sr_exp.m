%Georgios Tsiris, 1115201700173
clear; clc;
%disp('METHODOS SOR-Peiramatikos ipologismos fasmatikis aktinas');

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

CL=-tril(A, -1);

CU=-triu(A, 1);

I=eye(n);

D=diag(diag(A));

D1=inv(D);

L=D1*CL;

U=D1*CU;

B=L+U;
%disp('B');  disp(B);
%idiotimes tou pinaka B (Jacobi)
x=eig(B);     %disp(x);
rB=max(abs(x));
%disp('rB'); disp(rB);
omega=2.0/(1.0+sqrt(1-rB*rB));
disp('omega'); disp(omega);
U1=I-omega*U;
%disp(U1);
U1a=inv(U1);
%disp(U1a);
D_W=U1a*((1-omega)*I+omega*L);
%disp(D_W);
lambda=eig(D_W);
s_SOR=max(abs(lambda));
disp('theoretical_sr_SOR'); disp(s_SOR);

disp('*******************    SOR     *********************');
fid61=fopen('C:\Users\maria\Desktop\demo\SOR.txt','a+');
fprintf(fid61, 'w         r    \n\n');

col=0;
maria=0;
for omega=0.0:0.1:2.0
    U1=I-omega*U;
    %disp(U1);
    U1a=inv(U1);
    %disp(U1a);
    D_W=U1a*((1-omega)*I+omega*L);
    %disp(D_W);
    lambda=eig(D_W);
    s_SOR=max(abs(lambda));
    %disp('s_SOR'); disp(s_SOR);
    fprintf(fid61, '%f       %f    \n\n',omega, s_SOR);
    %apothikeusi apotelesmaton se pinakes gia dimiourgia
    %grafikis parastasis
    col=col+1;
    finw(col,1)=omega;
    fins(col,1)=s_SOR;
    %disp(finw);  disp(fins);
    plot(finw,fins,'r-');
    xlabel('\omega');              %  add axis labels and plot title
    ylabel('\rho(U_{\omega_0})');
    title('Behavior of \rho(U_{\omega_0}) with respect to \omega');
    disp('Ypologismos elaxistou apo peiramatiko prosdiorismo');
    [elaxS, thesi]=min(fins);
    wopt=finw(thesi,1); disp('wopt'); disp(wopt);
    disp('elaxS');disp(elaxS);
    disp('thesi');disp(thesi);
    %ypologismos olon ton theseon me tin idia peiramatiki fasmatiki
    %aktina
    %thesi=find(fins==elaxS);disp('thesi');disp(thesi);
    %wopt=finw(thesi,1); disp('wopt'); disp(wopt);
    maria=maria+1;
    %disp(maria);
end

fclose(fid61);
