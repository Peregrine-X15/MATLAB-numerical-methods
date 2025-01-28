function phi_star =TriDiagS(a,c,n,d,S_phi)

for i=2:n
xmult=a(i-1)/d(i-1);
d(i)=d(i)-xmult*c(i-1);
S_phi(i)=S_phi(i)-xmult*S_phi(i-1);
end

phi_star(n)=S_phi(n)/d(n);

for i=n-1:-1:1
phi_star(i)=(S_phi(i)-c(i)*phi_star(i+1))/d(i);
end
end