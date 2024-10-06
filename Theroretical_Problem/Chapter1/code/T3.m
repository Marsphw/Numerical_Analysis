x = zeros(10000, 1);
x(1) = -1;
times = 2;
while abs(x(times - 1) - x(times)) >= 1e-4
    x(times) = x(times - 1) - (4*power(x(times - 1), 3) - 2*power(x(times - 1), 2) + 3)/(12*power(x(times - 1), 2) - 4*x(times - 1));
    times = times + 1;
end
fprintf("Result: %.6f\n", x(times - 1));