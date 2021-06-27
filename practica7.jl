#Задача 1. Написать функцию pow(a, n::Integer), возвращающую значение a^n, и реализующую алгоритм быстрого возведения в степень.

function binpow(a, n::Integer)
    k = n
    t = 1
    p = a
    while k > 0
        if (mod(k,2)==0) # если степень четная
            k/=n
            p*=p
        else            # если степень нечетная
            k-=1
            t*=p
        end   
    end
    return t # возвращаем t = a^n
end

#Задача 2. Написать функцию fibonacci(n::Intrger), возвращающую n-ое число последовательности Фибоначчи, 
#имеющую оценку алгоритмической сложности O(log(n)), и не используя известную формулу Бине.

function fibonacci(n::Integer)
    a = b = c = 1
    f1,f2,f3 = 0,1,1
    fib_n = 0
    counter = 3
    while (n>counter)
        f1 = f2 + f3
        f2 = f3 + f1
        f3 = f1 + f2
        counter += 3
    end
    if n == counter
        fib_n = f3
    elseif n == (couner - 1)
        fib_n = f2
    else
        fib_n = f1      
    end
    return fib_n
end

#Задача 3. Написать функцию log(a::Real,x::Real,ε::Real), реализующую приближенное вычисление логарифма по основанию a>1 числа x>0 
#с максимально допустимой погрешностью e>0 (без использования разложения логарифмической функции в степенной ряд).

function log(a::Real, x::Real, e::Real) 
    z, t, y = x, 1, 0
    #ИНВАРИАНТ: a^y * z^t == x (=const)
    while z > a || z < 1/a || t > e
        if z > a
            z /= a
            y += t # т.к. z^t = (z/a)^t * a^t
        elseif z < 1/a
            z *= a
            y -= t # т.к. z^t = (z*a)^t * a^-t
        else # t > ε
            t /= 2
            z *= z # т.к. z^t = (z*z)^(t/2)
        end
    end
    # y: |log_a(x)-y| <= e
end

#Задача 4. Написать функцию isprime(n)::Bool, возвращающую значение true, если аргумент есть простое число, и - значение false, 
#в противном случае. При этом следует иметь ввиду, что число 1 простым не считается.

function isprime(n::Int)::Bool
    d=2
    while (d*d<=n)  #проверяем все целые числа, чьи квадраты меньше n
        if (n%d==0) 
            return false
        end
        d+=1
    end
    return true
end

#Задача 5. Написать функцию eratosphen(n), возвращающую вектор всех простых чисел, не превосходящих заданного натурального числа n.

function eratosphen(n::Int)
    a = zeros(Int, n)  # создаем нулевой вектор с длинной n
    for i in (1:n)  # заполням пустой массив нутуральными числами < n
        a.push(i)
    end
    b[]
    b[0] = 2
    b[1] = 3
    for i in 2:n
        if isprime(i) # Если число простое, то добавляем его динамический массив
            b[i] = i
        end
    end
    return b
end

# Задача 6. Написать функцию factor(n), получающую некоторое натуральное число n, и возвращающую кортеж, состоящий из вектора его простых делителей (в порядке возрастания) и вектора кратностей этих делителей, 
# т.е. выполняющую факторизацию заданного числа. Оценка вычислительной сложности алгоритма должна быть $O(\sqrt{n})$.

function factor(n)
    if (isprime(n)) #если число простое
        return n,1
    end

    dividers = []  # динамический вектор его простых делителей 
    multiplicity = [] # динамический вектор кратности его делителей
    
    divid = 2 # наименьшие возможный делитель
    n_1 = n 
    k = 0
    
    while (divid*divid<=n && n_1>1)
        if (isprime(divid))      # если делитель простой
            if (n_1%divid == 0)  # если проверяемое число является делителем числа n
                push!(dividers,divid)
                while (n_1%divid == 0) # если данный делитель встречатся больше одного раза
                    n_1 /= divid
                    k += 1
                end
                push!(multiplicity,k)
                k = 0
            end 
        end
        divid += 1
    end
    if (n_1 != 1)
        push!(dividers,Int(n_1))
        push!(multiplicity,1)
    end
    return dividers,multiplicity
end

# Задача 7. Написать фунуцию, получающую натуральный аргумент n, и возвращающую для него значение функции Эйлера.

function euler_function(n)
    if (isprime(n)) # Если число простое то fi(n) = n-1
        return n-1
    else
        a,b=factor(n) 
        result = 1
        for i in 1:length(a)
            if (b[i]==1)  # если делитель встречается 1 раз, то fi(n) = n - 1
                a[i]-=1
            else           # если делитель встречается  более 1-го раза, то fi(n) = n^k - n^(k-1), где k - кратность рассматриваемого делителя
                a[i] = a[i]^b[i] - a[i]^(b[i]-1)
            end
            result*=a[i] 
        end
        return result
    end
end

#Задача 8. Самостоятельно написать подобную функцию, реализующую расширенный алгоритм Евклида.

function ext_gcd(a, b)
    if a == 0
        return (b, 0, 1)
    else
        new_gcd  = gcd(b % a, a)
        div = new_gcd[1]
        x = new_gcd[2]
        y = new_gcd[3] 
    end
    return (div, y - (b // a) * x, x)
end

function gcd(a,b) # при этом a>b
    if b == 0
        return a
    end
    return gcd(b, a%b)
end

#Задача 9. Написать функцию inv(m::Integer,n::Integer) возвращающий обратный элемент к значению m в кольце вычетов по модулю n 
#(см. лекцию 3). При этом, если значение n не обратимо, то должно возвращаться значение nothing.

function inv(m::Integer,n::Integer)
    if (ext_gcd(m,n)>1)
        return nothing
    else
        return ext_gcd(m,n)
    end
end

#Задача 10. Написать функцию zerodivisors(m), возвращающую все делители нуля кольца вычетов по заданному модулю n.
function zerodivisors(n::Integer)
    for i in 1:n-1
        if gcd(i,n) != 1
            print(i)
        end
    end
end

#Задача 11. Написать функцию nilpotents(n), для заданного n возвращающую диапазон 
#(т.е. значение типа StepRange{Int64,Int64}), содержащий все не тривиальные нильпотенты кольца вычетов по модулю n. 
#Возвращаемый диапазон должен быть пустым, если нетривиальных нильпотентов в кольце нет.

function allnilpotents(n::Integer)::StepRange{Int64,Int64}
    result = []
    dividers,_ = factor(n) #в данном случае нам понадобиться только вектор  простых делителей 
    m_group = 1
    for i in 1:length(dividers)
        m_group *=dividers[i] #мультипликативная группа состоит из чисел взаимнопростых с основанием кольца
    end
    count = Int(n/m_group)
    for i in 1:count-1
        push!(result,m_group*i)
    end
    return result
end

#Залача 12. Написать функцию ord(a,p), возвращающую порядок заданного элемента a мультипликативной группы 
#кольца вычетов по заданному простому модулю p.

function ord(a,p)
    fi = p-1    # функция Эйлера
    result = 0
    for i in 1:(p-1)
        if (fi%i==0 && (a^i)%p==1)
            result = i
            break
        end
    end
    return result
end
