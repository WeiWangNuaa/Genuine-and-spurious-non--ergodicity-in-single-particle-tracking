function x=samples_FBM_resetting_without_memory(T,H,dt,r,x0)

N=ceil(T/dt);
num=unifrnd(0,1,N,1);
num=num<(dt*r); %1:resetting  0:non-resetting
num=find(num==1);%resetting points
t=num*dt; %resetting points time
n=length(t); % num of resetting points
x(num)=x0;

if n==0  %no resetting
    x=Wood_Chan_method(N,H,dt);
end


if n==1  % one resetting point
    
    if num(1)>1
       x(1:num(1)-1)=Wood_Chan_method((num(1)-1),H,dt);
    end

   if num(end)<N
    d=N-num(end);
    x(num(end)+1:N)=Wood_Chan_method((d),H,dt);
   end
   
end
    
if n>=2    % resetting points more than 1
    
  if num(1)>1
     x(1:num(1)-1)=Wood_Chan_method((num(1)-1),H,dt);
  end

  for i=2:n
    d=num(i)-num(i-1);
    if d>1
    x(num(i-1)+1:num(i)-1)=Wood_Chan_method((d-1),H,dt);
    end
  end

  if num(end)<N
    d=N-num(end);
    x(num(end)+1:N)=Wood_Chan_method((d),H,dt);
  end
  
end


end