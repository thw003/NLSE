y1 = @(x) (-x.^2);
y2 = @(x) (exp(-10*x.^2));
figure(1)
ezplot(y1);
grid on
figure(2)
ezplot(y2);
grid on
