
clc,clear

M=unique(ceil(logspace(0,3,20)));
dt=0.1;
t0=0.01;
n=200;
c=1;

alpha=0.7
for i=1:length(M)
    T=M(i)
    N=M(i)/dt;
    d=[1:1:N];
parfor k=1:n
[ X(k,:)]=generate_sample_RLMFBM_MEMORY_new(T,dt,alpha,d);
end
TAMSD=sum((X(:,(c+1):N)-X(:,1:N-c)).^2,2)./(N-c);
mean_TAMSD=mean(TAMSD,1);
EB(i)=var(TAMSD,1);
clearvars X TAMSD
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=1.4
for i=1:length(M)
    T=M(i)
    N=M(i)/dt;
    d=[1:1:N];
parfor k=1:n
[ X(k,:)]=generate_sample_RLMFBM_MEMORY_new(T,dt,alpha,d);
end
TAMSD=sum((X(:,(c+1):N)-X(:,1:N-c)).^2,2)./(N-c);
mean_TAMSD=mean(TAMSD,1);
EB(i)=var(TAMSD,1);
clearvars X TAMSD
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
loglog(M,EB1,'^r','markerfacecolor','r','MarkerSize',5)
hold on
loglog(M,EB2,'ob','markerfacecolor','b','MarkerSize',5)
d=[1e0 1e3]
loglog(d,0.1*d.^(-1),'k--')
loglog(d,0.0001*d.^(4*alpha-6),'k--')
xlabel('$T$','Interpreter','latex','Fontsize',10)
ylabel('EB','Interpreter','latex','Fontsize',10) 

set(gca,'FontSize',16);







