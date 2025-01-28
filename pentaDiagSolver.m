function phi_star =pentaDiagSolver(e,a,d,c,f,n,S_phi)

for i=2:n-1
    const1=a(i-1)/d(i-1);
    d(i)=d(i)-const1*c(i-1);
    S_phi(i)=S_phi(i)-const1*Q(i-1);
    const2=e(i-1)/d(i-1);
    a(i)=a(i)-const2*c(i-1);
    d(i+1)=d(i+1)-const2*f(i-1);
    
end





end