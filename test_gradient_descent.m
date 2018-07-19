function test_gradient_descent()


function [error] = compute_error(x,y,a,b,R)
  error = (sqrt((x-a)^2 + (y-b)^2) - R)^2;
endfunction

function [Result] = gradA(a,b,R,x,y)
  Result = 0;
  for i = 1:length(x)
      Numerator = -2 * (x(i) - a) * (sqrt((x(i) - a)^2 + (y(i) - b)^2) - R);
      Denominator = sqrt((x(i) - a)^2 + (y(i) - b)^2);
      Result = Result + Numerator / Denominator;
  end  
  Result;
endfunction


function [Result] = gradB(a,b,R,x,y)
  Result = 0;
  for i = 1:length(x)
      Numerator = -2 * (y(i) - b) * (sqrt((x(i) - a)^2 + (y(i) - b)^2) - R);
      Denominator = sqrt((x(i) - a)^2 + (y(i) - b)^2);
      Result = Result + Numerator / Denominator;
  end  
  Result;
endfunction


function [Result] = gradR(a,b,R,x,y)
  Result = 0;
  for i = 1:length(x)
    Result = Result + (-2*(sqrt((x(i) - a)^2 + (y(i) - b)^2) - R));
  end  
  Result;
endfunction


clear
clc
a = load("W:\\work\\LKARD\\GradienDescend\\eval_circle\\txt\\Q150F28000.txt")

average_a = 0;
average_b = 0;

x = a(:,1);
y = a(:,2);
f0 = a(:,3);

for i = 1:length(x)
  average_a = average_a + x(i)
  average_b = average_b + y(i)
end

average_a = average_a/length(x);
average_b = average_b/length(x);
average_R = 0;
sigma = 0.02;

for i = 1:length(x)
  average_R = average_R + sqrt((x(i) - average_a)^2 + (y(i) - average_b)^2);
end

average_R = average_R/length(x)

a = average_a;
b = average_b;
R = average_R;

for i = 1:1000
  error = 0;
  for i = 1:length(x)
    error = error + compute_error(x(i),y(i),a,b,R);
  end
  
  a = a - sigma * (gradA(a,b,R,x,y)/length(x));
  b = b - sigma * (gradB(a,b,R,x,y)/length(x));
  R = R - sigma * (gradR(a,b,R,x,y)/length(x));
end
a
b
R
error

for i = 1:length(x)
  x_ideal_circle(i) = cos(2*pi/(length(x)-1) * i)*R;
  y_ideal_circle(i) = sin(2*pi/(length(x)-1) * i)*R;
end

figure(1)
plot(x,y)
hold on;
plot(x_ideal_circle,y_ideal_circle, "--")
hold off;


average_a = average_a - a
average_b = average_b - b

x = x - a
y = y - b

##phi = 0;
##rotation_matrix = [cos(phi), sin(phi); -sin(phi), cos(phi)]
##
##for i = 1:length(x)   
##  dot_matrix = [x(i),y(i)]
##  result_matrix = dot_matrix*rotation_matrix;
##  x(i) = result_matrix(1);
##  y(i) = result_matrix(2);
##end
##
##for i = 1:length(x) 
##  x(i) = x(i)*cos(phi) - y(i)*sin(phi);
##  y(i) = y(i)*cos(phi) + x(i)*sin(phi); 
##end

figure(2)
plot(average_a,average_b, "+")
hold on;
plot(0,0,"*")
plot(x,y, "--");
hold off;
average_a
##------------------------------------------------------------------------------

function [Result] = getAlpha(Q,f,f0)
  Result = 2*atan(Q*(f/f0 - f0/f));
endfunction

function [Result] = WeightingFunction(Q,f,f0,x,y,R)
  Result = 0;
  for i = 1:length(x) # form new dots
    alpha = getAlpha(Q,f,f0(i));
    new_x(i) = cos(alpha)*R;
    new_y(i) = sin(alpha)*R;
  end
  
  for i = 1:length(x) # compute error
    Result = Result + sqrt((x(i) - new_x(i))^2 + (y(i) - new_y(i))^2);
  end
endfunction

function [Result] = gradQ(Q,f,f0,x,y,R)
  delta = 0.00001;
  Result = (WeightingFunction(Q+delta,f,f0,x,y,R) - WeightingFunction(Q,f,f0,x,y,R)) / delta;
endfunction

function [Result] = gradF(Q,f,f0,x,y,R)
  delta = 0.00001;
  Result = (WeightingFunction(Q,f+delta,f0,x,y,R) - WeightingFunction(Q,f,f0,x,y,R)) / delta;
endfunction

function [new_x, new_y] = formNewDots(Q,f,f0)
    for i = 1:length(x) # form new dots
    alpha = getAlpha(Q,f,f0(i));
    new_x(i) = cos(alpha)*R;
    new_y(i) = sin(alpha)*R;
    end
endfunction
  


dots = load("40_dots_after_c.txt")

c_x = dots(:,1);
c_y = dots(:,2);

##figure(10)
##plot(x,y)
##hold on;
##plot(c_x,c_y,"-*-")
##hold on;

x = c_x;
y = c_y;

Q = 6;
f = 29000;
sigma = 0.01;

for counter = 1:1000

  error = WeightingFunction(Q,f,f0,x,y,R)
  counter
  gradQ_v = gradQ(Q,f,f0,x,y,R)
  gradF_v = gradF(Q,f,f0,x,y,R)
  Q = Q - sigma*gradQ_v
  f = f - sigma*gradF_v
  
##  figure(3);
##  plot(x,y, "--*--");
##  hold on;
##  plot(new_x, new_y, "+")
##  hold off;
  
end

  for i = 1:length(x) # form new dots
    alpha = getAlpha(Q,f,f0(i));
    new_x(i) = cos(alpha)*R;
    new_y(i) = sin(alpha)*R;
  end
  
  figure(1)
  plot(new_x, new_y, "-*-")
  hold on;
  plot(x,y)
  hold off;
  
  figure(2)
  plot(x,y)
  hold on;
  Q = 49.979410848621114;
  f = 28964.626153654517
  formNewDots(Q, f, f0)
  plot(new_x, new_y, "--*--")
  
  error = WeightingFunction(Q,f,f0,x,y,R)
  
endfunction