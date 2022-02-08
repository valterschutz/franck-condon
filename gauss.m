function y = gauss(x,mu,sig)
y = exp(-1/2*((x-mu)./sig).^2);
end

