function divisor = gcd_vect(a, varargin)
    if ~isempty(varargin)
        a = [a, varargin{:}];
    elseif length(a) == 1
        divisor = a;
        return;
    end
    divisor = gcd(a(1), gcd_vect(a(2:end)));
end