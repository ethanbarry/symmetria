set terminal pngcairo size 1920,1200
set xrange [0:1]
set yrange [0:1]

# Define our integrand.
f(x) = x * x * x * x / sqrt(2 * (1 + x * x))


# Evaluation points
x1 = 0.0
x2 = 1.0
x3 = 1/2.0 - (1/(2 * sqrt(5)))
x4 = 1/2.0 + (1/(2 * sqrt(5)))

# Evaluate them.
f1 = f(x1)
f2 = f(x2)
f3 = f(x3)
f4 = f(x4)

plot f(x1) w l, f(x2) w l, f(x3) w l, f(x4) w l, f(x) with filledcurves y=0 fc rgbcolor "#33658a" fs transparent solid 0.5
