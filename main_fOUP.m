clc
clear

%% change Hurst exponent
H=0.5

T=500;
dt=0.01
n=1000;
N=T/dt;
t=10;
lambda=1;
sigma=1;
x0=0;

parfor j=1:n
X(j,:)=samples_fOUP(N,dt,lambda,sigma,T,H,x0);   
end

u=t/dt;
c=unique(ceil(logspace(-2,2,30)/dt));
e=c(find(c<=(N-u)));

parfor i=1:length(c)
    
TAMSD_X(:,i)=sum((X(:,1+c(i):N)-X(:,1:N-c(i))).^2,2)./(N-c(i));

end
mean_TAMSD_X=mean(TAMSD_X,1);
EAMSD_X=mean(X.^2,1);



Y=X(:,u+e)-X(:,u);
MSI_X=mean(Y.^2,1);

h=figure

%TAMSDs
loglog(c*dt,TAMSD_X(2:50,:),'m')
hold on
h4=loglog(c*dt,TAMSD_X(1,:),'m')
%EA-TAMSD
h1=loglog(c*dt,mean_TAMSD_X,'b-','LineWidth',2)
%MSD
h2=loglog(c*dt,EAMSD_X(c),'ko','MarkerFaceColor','k')
%% MSI
h3=loglog(e*dt,MSI_X,'gs','MarkerFaceColor','g','markersize',5)
hold on

%%%% plateau of MSD and MSI
yline(gamma(2*H+1)/2,'k:','LineWidth',1.5);
yline(gamma(2*H+1),'k:','LineWidth',1.5);

legend([h2 h3 h1 h4],{'MSD','MSI','mean TAMSD','TAMSDs'},'Interpreter','latex','Fontsize',13)
legend('boxoff')

xlabel('$t$','Interpreter','latex','Fontsize',13)
ylabel('MSD, MSI, TAMSD','Interpreter','latex','Fontsize',13) 
set(gca,'FontSize',16);
if H==0.8
ylim([1e-4 1e2])
end
if H==0.2
ylim([1e-2 1e1])
end
if H==0.5
ylim([1e-3 1e1])
end
xlim([1e-2 1e2])

%%Save data
filename=['fOUP_trajectory_H=' num2str(H) '.mat']
save(filename,'X','H','dt','T','t','sigma','lambda','x0','N','n')

