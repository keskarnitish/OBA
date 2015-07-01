function output=soft_threshold_l1(y,lambda)
output=sign(y).*max(abs(y)-lambda,0);
end