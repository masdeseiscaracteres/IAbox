function P=function_power(H,V,U)

K=size(H,1); %number of users
P=zeros(K,1);
for u=1:K
P(u)=trace(U{u}'*H{u,u}*V{u}*V{u}'*H{u,u}'*U{u});
end