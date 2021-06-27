#Задача 1. Написать функцию, вычисляющую n-ую частичную сумму ряда Тейлора функции 
#cos(x)=1-\frac{x^2}{2!}+\frac{x^4}{4!}-\frac{x^6}{6!}+... для заданного значения аргумента x. 
#Вычислительная сложность алгоритма должна иметь оценку O(n).

function cos(n,x)
    sum = 1 # сумма ряда Тейлора для cos(x)
    frac = 1 # коэффициент перед параметром
    for i in 1:n
        frac*=i 
        if (i%2==0)
            sum+=(((-1)^((i/2)%2))*(x^i))/frac # применяем формулу Тейлора для  cos(x) = sum([(-1)^n/(2n)!] * x^(2n+1))
        end
    end
    return sum
end

# Задача 2. Написать функцию, вычисляющую значение суммы ряда Тейлора функции cos(x) в заданной точке с машинной точностью.

function cos(x,e) # где e - Эпсилон
    sum=1
    frac=1
    i=1
    a=1
    while(abs(a)>e)
        frac*=i
        if (i%2==0)
            a=(((-1)^((i/2)%2))*(x^i))/frac
            sum+=a
        end
    i+=1
    end
    return sum
end

#Задача 3. Построить семейство графиков n-ых частичных сумм ряда Тейлора функции cos(x) на одном её периоде, для n=2,4,8,16.
function cos_graphs()
    f1(x) = 1-(x^2)/2  # частичная сумма ряда Тейлора функции cos(x) при n = 1
    f2(x) = f1(x) + (x^4)/factorial(4)  # частичная сумма ряда Тейлора функции cos(x) при n = 2
    f3(x) = f2(x) - (x^6)/factorial(6) + (x^8)/factorial(8) # частичная сумма ряда Тейлора функции cos(x) при n = 3
    f4(x) = f3(x) - (x^10)/factorial(10) + (x^12)/factorial(12) - (x^14)/factorial(14) + (x^16)/factorial(16) # частичная сумма ряда Тейлора функции cos(x) при n = 4
    p=plot(f1)  # строим семейство графиков
    plot!(f2) 
    plot!(f3)
    plot!(f4)
end
#Задача 5. Следующий степенной ряд определяет семейство так называемых функций Бесселя 1-го рода порядка $m$ ($m=0,1,2,...$) 
#J_m(x)=\Big(\frac{x}{2}\Big)^m\sum_{k=0}^{\infty}\frac{(-1)^k}{k!(k+m)!}\Big(\frac{x}{2}\Big)^{2k} 
#Написать функцию besselj(m,x), вычисляющую функцию Бесселя 1-го рода порядка $m$ в точке $x \in \mathbb{R}$ 
#с машинной точностью, и построить семейство графиков для $m=0,1,2,3,4,5$ (вид требуемых графиков см., например, здесь).

function bessel(m,x)
    sum=1/factorial(m)
    i=1
    a=1
    while(abs(a)>ε)
        a*=((-1)/(i*(i+m)))*(x/2)*(x/2)
        sum+=a
        i+=1
    end
    sum*=(x/2)^m
    return sum
end


# Задача 6. Написать функцию linsolve(A,b), получающую на вход невырожденную квадратную верхнетреугольную матрицу A (матрицу СЛАУ, приведенную к ступенчатому виду), 
# вектор-столбец b (правую часть СЛАУ), и возвращающую решение соответствующей СЛАУ.
# Пояснение. Матрица называется верхнетреугольной, если все её элементы, стоящие ниже главной диагонали равны нулю.

function linsolve(a,b)
    x = collect(length(a))

    for i in length(a):1
        for j in length(a)-1:i
            for k in length(a)-2:i
            if a[j][k]>=0 & a[i][i]<=0 || a[j][k]<=0 & a[i][i]>=0
                a[j][i] = a[j][i] + a[i][i]*a[j][i]/a[i][i]
            else
                a[j][i] = a[j][i] - a[i][i]*a[j][i]/a[i][i]
            end
        end
    end

    for i in eachindex
        x[i] = b[i]/a[i][i]
    end

    return x
end

#Задача 7. Написать функцию convert!(A), получающую на вход прямоугольную матрицу (например, - расширенную матрицу СЛАУ) 
#и пробразующую эту матрицу к ступенчатому виду с помощью элементарных преобразований строк.

function convert!(A)
    for i in 1:length(A)
        for j in length(A):-1:i+1
            A[j,:]-=(A[j][i]/Ab[i][i])*Ab[i,:] # Делаем срез по столбцу
        end
    end
    return A
end

function det(A)
    m = lengtn(A)
    n = length(A[0]) # n - минор элемента A[0][0]
    if m != n
        return None
    if n == 1
        return A[0][0]
    detA = 0
 
    for j in range(n)
        detA += A[0][j] * pow(-1,j) * det(minor(A, 0, j)) # формула детерминанта по первой строке: det = sum((-1)^k * a[1][k] * M[1][k]), где k = 1:n, M - минор матрицы
        signum *= -1
    end
    return detA
end

#Задача 9. Написать функцию inv(A), получающую на вход квадратную матрицу, и возвращающую обратную матрицу, 
#если матрица обратима, или - значение nothing, в противном случае.

function inv(A)
    if (det_(A)==0) # если детерминант равен 0 , то не сущ. обратной матрицы
        return nothing 
    end
    B=Matrix{Float64}(undef,length(A),2*length(A)) # создаем новую вещественную матрицу
    for i in 1:length(A)
        for j in 1:(2*length(A))
            if (i <= length(A))
                B[i][j]=A[i][j]
            elseif (i==j)
                B[i][j]=1
            else
                B[i][j]=0
            end
        end
    end
    convert!(B)     # convert!() - cм. Задача 7/практика 9
    linsolve(B,zeros(length(A))) #linsolve() - cм. Задача 6/практика 9
    for i in 1:length(A)
        B[i][:] = (1/B[i][i])*B[i][:]
    end
    return B[:,(length(A)+1):(2*length(A))] # возвращаем обратную матрицу
end

#Задача 10. Написать функцию rang(A), получающую на вход матрицу (вообще говоря, прямоугольную), и возвращающую её ранг.

function rang(A)
    B=copy(A)
    rang = length(B)
    while (det_(B)==0 && rang!=1)
        rang-=1
        B=B[1:rang]
    end
    return rang
end
