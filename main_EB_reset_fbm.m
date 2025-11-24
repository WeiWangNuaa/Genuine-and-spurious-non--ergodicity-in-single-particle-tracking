
clc,clear

M=unique(ceil(logspace(0,3,20)));
dt=0.1;
t0=0.01;
n=1000;
r=1;
x0=0;
c=1;

H=0.2
for i=1:length(M)
    T=M(i);
    N=M(i)/dt;
    d=[1:1:N];
parfor k=1:n
[ X(k,:)]=samples_FBM_resetting_without_memory(T,H,dt,r,x0);
end
TAMSD=sum((X(:,(c+1):N)-X(:,1:N-c)).^2,2)./(N-c);
mean_TAMSD=mean(TAMSD,1);
EB1(i)=var(TAMSD,1);
clearvars X TAMSD
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=0.5
for i=1:length(M)
    T=M(i);
    N=M(i)/dt;
    d=[1:1:N];
parfor k=1:n
[ X(k,:)]=samples_FBM_resetting_without_memory(T,H,dt,r,x0);
end
TAMSD=sum((X(:,(c+1):N)-X(:,1:N-c)).^2,2)./(N-c);
mean_TAMSD=mean(TAMSD,1);
EB2(i)=var(TAMSD,1);
clearvars X TAMSD
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=0.8
for i=1:length(M)
    T=M(i);
    N=M(i)/dt;
    d=[1:1:N];
parfor k=1:n
[ X(k,:)]=samples_FBM_resetting_without_memory(T,H,dt,r,x0);
end
TAMSD=sum((X(:,(c+1):N)-X(:,1:N-c)).^2,2)./(N-c);
mean_TAMSD=mean(TAMSD,1);
EB3(i)=var(TAMSD,1);
clearvars X TAMSD
end



figure
loglog(M,EB1,'^r','markerfacecolor','r','MarkerSize',5)
hold on
loglog(M,EB2,'ob','markerfacecolor','b','MarkerSize',5)
loglog(M,EB3,'sg','markerfacecolor','g','MarkerSize',5)
d=[1e0 1e3]
loglog(d,0.1*d.^(-1),'k--')
% loglog(d,0.0001*d.^(4*alpha-6),'k--')
xlabel('$T$','Interpreter','latex','Fontsize',10)
ylabel('EB','Interpreter','latex','Fontsize',10) 
set(gca,'FontSize',16);







