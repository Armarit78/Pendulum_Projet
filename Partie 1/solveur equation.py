from sympy import symbols, Eq, sin, cos, solve

m1, m2, l1, l2, g, th1, th2, dth1, dth2, d2th1, d2th2 = symbols(
    'm1 m2 l1 l2 g th1 th2 dth1 dth2 d2th1 d2th2')

eq1 = Eq(((m1 + m2) * l1 * d2th1 ) + (m2 * l2 * d2th2 * cos(th1 - th2)) + (m2 * l2 * (dth2)**2 * sin(th1 - th2)) + ((m1 + m2) * g * sin(th1)), 0)

eq2 = Eq((l1 * d2th1 * cos(th1 - th2)) + (l2*d2th2) - (l1 * (dth1)**2 * sin(th1 - th2)) + (g * sin(th2)), 0)


solutions = solve([eq1, eq2], (d2th1, d2th2), simplify=True, rational=True)
print(solutions[d2th1])
print(solutions[d2th2])