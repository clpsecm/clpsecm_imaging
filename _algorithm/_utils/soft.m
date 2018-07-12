function y = soft(x,lda)
%SOFT(x,lda) Soft threshold of x with threshold lda
if lda<0 
    error('lda should be non-negative'); 
else
    y = sign(x).*(max(abs(x)-lda,0));
end




