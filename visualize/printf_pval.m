function p_str = printf_pval(pval)

if pval<0.001
    p_str = 'p<0.001';
else
   p_str = sprintf('p=%.3f',pval); 
end