clc
clear


% num=1:exponential waiting times;num=2:power-law waiting times
num=2;
a=1.5;

T=500;
dt=0.1;
n=2000;
N=T/dt;
t=10;
v=0;
tau=0.01;

parfor j=1:n
X(j,:)=generate_CTRW(T,dt,a,tau,num);
% X(j,:)=generate_subordinated_velocity(T,dt,a,ds,v);
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
% loglog(c*dt,TAMSD_X(2:50,:),'m')
hold on
h4=loglog(c*dt,TAMSD_X(1,:),'m')
h1=loglog(c*dt,mean_TAMSD_X,'b-','LineWidth',2)
h2=loglog(c*dt,EAMSD_X(c),'ko','MarkerFaceColor','k')
hold on
box on
xlim([dt 100])

% MSI
h3=loglog(e*dt,MSI_X,'gs','MarkerFaceColor','g','markersize',7)
hold on

% s=1000;
% k=sin(pi*a)/(pi*a*tau^a);
% k=1/gamma(a+1);
% loglog([dt s*dt],k*[dt s*dt].^(a),'--k','LineWidth',1.5)
% loglog([dt s*dt],k*[dt s*dt].^(1)/t^(1-a),'--k','LineWidth',1.5)
% loglog([dt s*dt],k*[dt s*dt].^(1)/T^(1-a),'--k','LineWidth',1.5)

% k1=(a-1)/(tau)*0.1;
% loglog([dt s*dt],k1*[dt s*dt].^(1),'--k','LineWidth',1.5)
% loglog([dt s*dt],k1*[dt s*dt].^(1),'--r','LineWidth',1.5)
% loglog([dt s*dt],k1*[dt s*dt].^(1),'--k','LineWidth',1.5)

% k=1.5;
% loglog([dt s*dt],k*[dt s*dt].^(1),'--k','LineWidth',1.5,'HandleVisibility','off')
% loglog([dt s*dt],k*[dt s*dt].^(1),'--k','LineWidth',1.5)
% loglog([dt s*dt],k*[dt s*dt].^(1),'--k','LineWidth',1.5)

% 
% mean_TAMSD_X./(k*[c].^(1)/T^(1-a))
% EAMSD_X(c)./(k*[c].^(a))
% EAMSD_Y./(k*[e].^(1)/t^(1-a))

% yline(EAMSD_theory_limted_X,'-.k','LineWidth',1.5);
% yline(2*EAMSD_theory_limted_X,'-.k','LineWidth',1.5);

xlabel('$\Delta$','Interpreter','latex','Fontsize',13)
ylabel('MSD, MSI, TAMSD','Interpreter','latex','Fontsize',13) 
set(gca,'FontSize',16);
legend([h2 h3 h1 h4],{'MSD','MSI','mean TAMSD','TAMSDs'},'Interpreter','latex','Fontsize',13)
legend('boxoff')
% ylim([1e-2 1e4])
% 
%%Save data
filename=['CTRW_trajectory_a=' num2str(a) '.mat']
save(filename,'X','a','dt','T','t','tau','N','n')

