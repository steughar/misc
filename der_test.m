x = 2;
y = 3;

f = x^2 + y^3;
h = 0.001;

dF_x = (((x+h)^2 + y^3) - (x^2 + y^3))/h
dF_y = ((x^2 + (y+h)^3) - (x^2 + y^3))/h

gradX = 2*x
gradY = 3*y^2



