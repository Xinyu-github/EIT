function sig_mix = conductivity_mixture(sig1,sig2,alpha)
%     Input:
%       sig1:       conductivity of first component
%       sig2:       conductivity of second component
%       alpha:      ratio of two components. alpha = 1 for pure 2nd component

%     Output:
%       sig_mix:    conductivity of mixture 

    nominator = 2*sig2+sig1-2*alpha*(sig2-sig1);
    denominator = 2*sig2+sig1+alpha*(sig2-sig1);
    sig_mix = sig2*nominator/denominator;

end