clear
a = load("40_dots.txt")

average_a = 0;
average_b = 0;

x = a(:,2)
y = a(:,3)

for i = 1:length(a(:,1))
  average_a = average_a + x(i)
  average_b = average_b + y(i)
end

average_a = average_a/length(x)
average_b = average_b/length(x)

plot(x,y)
hold on;
plot(average_a, average_b, "+")  
hold on;


%------------------------------------------------


a_ideal_circle = 22.0696
b_ideal_circle = 8.5345 
R_ideal_circle = 22.5548  
plot(a_ideal_circle, b_ideal_circle, "*");
hold on;

a = a_ideal_circle
b = b_ideal_circle
R = R_ideal_circle

for i = 1:length(x)
  x_ideal_circle(i) = cos(2*pi/(length(x)-1) * i)*R + a;
  y_ideal_circle(i) = sin(2*pi/(length(x)-1) * i)*R + b;
end

plot(x_ideal_circle, y_ideal_circle)
hold off;

x = x - a
y = y - b
plot(x,y)
hold on;

phi = pi/3
rotation_matrix = [cos(phi), -sin(phi); sin(phi), cos(phi)]

for i = 1:length(x)   
  dot_matrix = [x(i),y(i)]
  result_matrix = dot_matrix*rotation_matrix;
  x(i) = result_matrix(1);
  y(i) = result_matrix(2);
end

plot(x,y,'----')



%--------------------- util functions
function [gradQ] = gradQ(Q, f, f0)
  numeratorQ = 2*(f/f0 - f0/f);
  denominatorQ = Q^2*(f/f0 - f0/f)^2 + 1;
  gradQ = -(numeratorQ/denominatorQ);
endfunction

function [gradF] = gradF(Q, f, f0)
  numeratorF = 2*Q*(f0/f^2 - 1/f0);
  denominatorF = Q^2*(f/f0 - f0/f)^2 + 1;
  gradF = -(numeratorF/denominatorF);
endfunction

function[x_circle, y_circle] = update_circle(x_circle, y_circle, N, R, f0) 
  for i = 1: N
    alpha = phi - 2*atan(Q*(f/f0(i) - f0(i)/f));
    x_circle(i) = R*cos(alpha) + A;
    y_circle(i) = R*sin(alpha) + B;
  end
endfunction

function [x1] = x1(R, Q, f, f0)
  x1 = R*cos(2*atan(Q*(f/f0 - f0/f)))
endfunction

function [y1] = y1(R, Q, f, f0)
  y1 = R*sin(2*atan(Q*(f/f0 - f0/f)))
endfunction


%---------------------------------------










