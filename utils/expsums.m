function [a, b, f] = expsums(n, method)
%EXPSUMS Approximate the inverse through exponential sums
%
% 1 / X ~ A(1) * EXP(-B(1) * X) + ... + A(N) * EXP(-B(N) * X), 
%
% where Re(X) > 0. 

if ~exist('method', 'var')
	method = 'sinc';
end

switch method
		
	case 'sinc'
		
		n = floor(n/2);
		l = (-n:n);
		
		alpha = 1;
		d = 2;
		h = sqrt(2*pi*d/alpha/n);
		a = h ./ (1 + exp(-l*h));
		b = log(1 + exp(l*h));
        
    case 'sinc2'
        n = floor(n/2);
        l = (-n:n);
        
        alpha = 1;
        d = 3*pi/4;
        h = sqrt(2*pi*d/alpha/n);
		a = h ./ (1 + exp(-l*h));
		b = log(1 + exp(l*h));
        
    case 'fourier'
        l = (1 : n);
        h = pi / (n+1);
        
        t = -pi/2 + l .* h;
        a = h ./ (1 + exp(tan(t))) ./ (cos(t).^2);
        b = log(1 + exp(-tan(t)));
        
    otherwise
        error('Unsupported expsum method');
		
end

f = @(x) expsum_evaluate(a, b, x);

	function r = expsum_evaluate(a, b, x)
		r = a(1) * expm(-b(1) * x);
		for jj = 2 : length(b)
			r = r + a(jj) * expm(-b(jj) * x);
		end
	end

end

