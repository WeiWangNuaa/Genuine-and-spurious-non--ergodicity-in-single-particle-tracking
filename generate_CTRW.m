function y=generate_CTRW(T,dt,a,tau0,num)

beta=a;
N=T/dt;
x(1)=0;
n=1;


while n<=N+1

if num==2
u=unifrnd(0,1);
tau=tau0*((1-u).^(-1/beta)-1);
n_tau = round(tau / dt);

x(n+1:n+1+n_tau)=x(n);
n=n+1+n_tau;
x(n)=x(n-1)+randn;
end


if num==1
tau=exprnd(a);
n_tau=round(tau/dt);

x(n+1:n+1+n_tau)=x(n);
n=n+1+n_tau;
x(n)=x(n-1)+randn;
end

end

y=x(2:N+1);





end

