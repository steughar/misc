for i = 1:100
  x(i) = 2*(i-50) - (i-50)^2/5 + 3;
  tx(i) = (i-50)
endfor

plot(tx,x)
