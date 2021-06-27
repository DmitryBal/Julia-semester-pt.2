# Задача 1. Реализовать метод Ньютона, написав функцию со следующим заголовком => newton(r::Function, x; ε_x=1e-8, ε_y=1e-8, nmaxiter=20)
function newton(f::Function,a=1e-8, derrirative=20)
    x0 = a
    xn = f(x0)
    y = xn - f(xn) / derrirative(xn) # метод Ньютона(касательной): y - f(x0) = f'(x0)(x-x0). В качестве исходной точки x0 выбирается тот конец интервала, которому отвечает ордината того же знака, что и знак f''(x)
    while abs(y-xn) > math.pow(10,-3)
        xn = y
        y = xn - f(xn) / derrirative(xn)
    end

    return y

end

#Задача 2. Решить с помощью этой функции (newton) уравнение $\cos(x)=x$.
newton(x->(x-cos(x))/(1+sin(x)), 0.5)

#Задача 3. Реализовать ещё один метод функции newton со следующим заголовком
newton(ff::Tuple{Function,Function}, x; ε_x=1e-8, ε_y=1e-8, nmaxiter=20) = newton(x->ff[1](x)/ff[2](x), x; ε_x, ε_y, nmaxiter)

#Задача 4. Решить с помощью этого варианта функции newton уравнение $\cos(x)=x$.
newton((x->x-cos(x), x->1+sin(x)), 0.5) 

#Задача 5. Реализовать еще один метод функции newton со следующим заголовком
newton(ff, x; ε_x=1e-8, ε_y=1e-8, nmaxiter=20) = newton(x->(y=ff(x); y[1]/y[2]), x; ε_x, ε_y, nmaxiter)

#Задача 6. Решить с помощью последнего варианта функции newton уравнение $\cos(x)=x$.
newton(x->(x-cos(x),1+sin(x)), 0.5) 