function p = normalpdf(x,mu,si) 
p = (1/sqrt(2*pi*si^2))*exp(-((x-mu).^2)./(2*si^2)); 
return
