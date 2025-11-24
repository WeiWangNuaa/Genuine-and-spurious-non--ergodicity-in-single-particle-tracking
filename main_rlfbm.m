clc
clear

%parameter
alpha=1.45;


T=300;
dt=0.1;
N=T/dt;
n=500;
t=10;
V=[1:1:N];

parfor j=1:n
X(j,:)=generate_sample_RLMFBM_MEMORY_new(T,dt,alpha,V);   
end

c=unique(ceil(logspace(log10(dt),2,30)/dt));
parfor i=1:length(c)
    
TAMSD_X(:,i)=sum((X(:,1+c(i):N)-X(:,1:N-c(i))).^2,2)./(N-c(i));

end
mean_TAMSD_X=mean(TAMSD_X,1);
EAMSD_X=mean(X.^2,1);


u=t/dt;
e=c(find(c<=(N-u)));
Y=X(:,u+e)-X(:,u);
MSI_X=mean(Y.^2,1);

h=figure
loglog(c*dt,TAMSD_X(2:50,:),'m')
hold on
h4=loglog(c*dt,TAMSD_X(1,:),'m')
h1=loglog(c*dt,mean_TAMSD_X,'b-','LineWidth',2)
h2=loglog(c*dt,EAMSD_X(c),'ko','MarkerFaceColor','k')
hold on
box on

% MSI
h3=loglog(e*dt,MSI_X,'gs','MarkerFaceColor','g','markersize',7)
hold on

s=20;
loglog([dt s*dt],5*[dt s*dt].^(2*alpha-1),'--k','LineWidth',2)


xlabel('$\Delta$','Interpreter','latex','Fontsize',13)
ylabel('MSD, MSI, TAMSD','Interpreter','latex','Fontsize',13) 
set(gca,'FontSize',16);

legend([h2 h3 h1 h4],{'MSD','MSI','mean TAMSD','TAMSDs'},'Interpreter','latex','Fontsize',13)
legend('boxoff')
ylim([1e-5 1e5])
xlim([1e-1 1e2])

filename=['RLFBM_trajectory_a=' num2str(alpha) '.mat']
save(filename,'X','alpha','dt','T','t','N','n')




