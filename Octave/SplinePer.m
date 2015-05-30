function z = SplinePer (x, y) 
    n=length(x);
    twoes=2*ones(1,n);
##    twoes[1] = twoes[1] +2;
    h=diff(x);
    lambda=zeros(1,n-1); ##
    for i=1:(n-2) 
    	lambda(i+1)=h(i+1)/(h(i+1)+h(i));
    endfor
    mu=zeros(1,n-1);
    for i=1:(n-2)
    	mu(i)=mu(i)+h(i)/(h(i+1)+h(i));
    endfor
#    mu(n-1)=1;

    A=1*diag(mu,-1)+1*diag(twoes)+1*diag(lambda,1); 
    [fil col] = length(A);
    A(1, fil) = h(n-1)/(h(n-1)+h(n-2));
    
    dd1=zeros(n-1,1);
    for i=1:(n-1)
   		dd1(i)=(y(i+1)-y(i))/(x(i+1)-x(i)); #Primeras derivadas
    endfor
    dd2=zeros(n-2,1);
    for i=1:(n-2)
    	dd2(i)=(dd1(i+1)-dd1(i))/(x(i+2)-x(i)); #Segundas derivadas
    endfor
    
    A(1,:)=h(1)*A(1,:);   #Ajustes de la matriz
    A(1,1)=A(1,1)-h(1)/3;
    A(1,2)=A(1,2)+h(1)/6;
    A(n,:)=h(n-1)*A(1,:);
    A(n,1)=A(1,1)+h(1)/3;
    A(n,2)=A(1,2)-h(1)/6;
    gamma=zeros(n-1,1);
    gamma(1)=(h(1)-dd1(1))/h(1) + (dd1(1)-h(n))/h(n);
    for i=1:(n-2)
    	gamma(i+1)=6*dd2(i);
    endfor
    m=A\gamma;
    s = ppint(ppint(SplineLineal(x, m)));  

endfunction
