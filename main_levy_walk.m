clc
clear


% num=1:exponential waiting times;num=2:power-law waiting times
num=2;
%parameter
a=1.2;

T=500;
dt=0.01
t0=0.01;
n=300;
N=T/dt;
t=100;
x0=0

parfor j=1:n
X(j,:)=generate_levy_walk(T,dt,a,t0,num,x0);   
end

c=unique(ceil(logspace(-2,2,30)/dt));
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

%%a<1
% s=10000;
% k=2*(a-1)/(3-a)/(2-a)*t0^(a-1);
% loglog([dt s*dt],k*[dt s*dt].^(3-a),'--k','LineWidth',1.5)
% loglog([dt s*dt],(1-a)*[dt s*dt].^2,'--k','LineWidth',1.5)
% loglog([dt s*dt],[dt s*dt].^(2),'--k','LineWidth',1.5)

%a>1
s=10000;
k=2*(a-1)/(3-a)/(2-a)*t0^(a-1);
loglog([dt s*dt],k*[dt s*dt].^(3-a),'--k','LineWidth',1.5)
msd=k/(a-1)*(-1-(2-a)*(1+c*dt/t).^(3-a)+(3-a)*(1+c*dt/t).^(2-a)+(c*dt/t).^(3-a))*t^(3-a);
loglog(c*dt,msd,'k--','LineWidth',2)
% yline(EAMSD_theory_limted_X,'-.k','LineWidth',1.5);
% yline(2*EAMSD_theory_limted_X,'-.k','LineWidth',1.5);

%%exponnetial
% s=10000;
% % loglog([dt s*dt],(1-a)*[dt s*dt].^2,'--k','LineWidth',1.5)
% loglog([dt s*dt],[dt s*dt].^(2),'--k','LineWidth',1.5)

xlabel('$\Delta$','Interpreter','latex','Fontsize',13)
ylabel('MSD, MSI, TAMSD','Interpreter','latex','Fontsize',13) 
set(gca,'FontSize',16);
legend([h2 h3 h1 h4],{'MSD','MSI','mean TAMSD','TAMSDs'},'Interpreter','latex','Fontsize',13)
legend('boxoff')
xlim([1e-2 1e2])
ylim([1e-6 1e4])

filename=['LW_trajectory_a=' num2str(a) '.mat']
save(filename,'X','a','dt','T','t','t0','x0','N','n')
% 

