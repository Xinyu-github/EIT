function out=complex_truncation(input)

% Gives the truncation of a complex number in (or a vector of complex numbers):
%   - If the real part of this number is positive, then it returns this
%   number
%   - otherwise it will return imag(input) (same imaginary part, real part =0)
out=input;
out(real(out)<0)=1i*imag(out(real(out)<0));