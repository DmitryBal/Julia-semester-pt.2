struct Poly{T}
    coeff::Vector{T}
end

import Base.+
import Base.*
import Base.%
import Base.÷
import Base./

function +(p::Poly,q::Poly)::Poly
    if (length(p.coeff)<length(q.coeff))
        p,q=q,p
    end
    v = copy(p.coeff) 
    for i in 1:length(q.coeff)
        v[i]=v[i]+q.coeff[i]
    end
    return Poly{Int}(v)
end

function *(q::Int,p::Poly)::Poly
    v = copy(p.coeff)
    for i in 1:length(p.coeff)
        v[i]*=q
    end
    return Poly{Int}(v)
end

function *(p::Poly,q::Poly)::Poly
    v = Array{Int}(undef,length(p.coeff)+length(q.coeff)-1)
    for i in 1:(length(p.coeff)+length(q.coeff)-1)
        v[i]=0
    end
    for i in 1:(length(p.coeff)+length(q.coeff))
        for j in 0:(i-1)
            if (j+1<=length(p.coeff) && i-j<=length(q.coeff))
                v[i]+=p.coeff[j+1]*q.coeff[i-j]
            end                
        end
    end
    return Poly{Int}(v)
end

#Задача 1. Аналогично построить индуктивную функцию и реализовать соответствующий программный код, 
#вычисляющий значение второй производной многочлена в точке.

function diff_2(x,A)
    Q′′=0
    Q′=0
    Q=0
    for a in A
        Q′′=Q′′*x+2*Q′
        Q′=Q′*x+Q
        Q=Q*x+a
    end
    return Q′′
end

#Задача 2. Сделать то же самое для 3-ей производной многочлена, заданного своими коэффициентами

function diff_3(x,A)
    Q′′′=0
    Q′′=0
    Q′=0
    Q=0
    for a in A
        Q′′′=Q′′′*x+3*Q′′
        Q′′=Q′′*x+2*Q′
        Q′=Q′*x+Q
        Q=Q*x+a
    end
    return Q′′′
end

#Задача 3. Сделать то же самое для k-ой производной многочлена, заданного своими коэффициентами 
#(здесь предполагается, что k - это параметр функции evaldiffpoly)

function diff_k(x,A,k)
    arr = zeros(k+1)
    diff = Poly{Int}(arr)
    add = 0
    for a in A
        for i in length(diff.coeff):-1:1  
            if (i!=1)
                diff.coeff[i]=diff.coeff[i]*x+(i-1)*diff.coeff[i-1]
            else
                diff.coeff[i]=diff.coeff[i]*x+a
            end
        end
    end
    return diff.coeff[k+1]
end

#Задача 4. Определить функцицию с именем diff, для пользовательского типа Polynom, так чтобы ее можно было бы использовать так:
#при этом по умолчанию именованный аргумент ord сделать равным 1.

function diff(p,x; ord = 1)  # ord - порядок производной, p=Polynom([....])
    return diff_k(x,p.coeff,ord)
end

#Задача 5. Реализовать функцию divrem(a::AbstractVector,b::AbstractVector), осуществляющую деление многочлена на многочлен, 
#получающую на вход два массива с коэффициентами этих многочленов, и возвращающую другие два массива с коэффициентами 
#частного и остатка (имя divrem выбрано не случайно - имеется аналогичная встроенная функция с тем же именем, возвращающая 
#частное и остаток от деления целых чисел).

function divrem(a::AbstractVector,b::AbstractVector)
    c_a = copy(a)
    d = length(a)-length(b)  #разница между степенями многочленов
    res = zeros(length(a)-length(b)+1)
    for i in 1:length(b)
        res[i]=c_a[i]/b[1]
        for j in 1:length(b)
            c_a[i+j-1]=c_a[i+j-1]-res[i]*b[j]
        end
        c_a[i] = 0
    end
    return res,c_a #c_a - остаток
end

#Задача 6. Соответствующим образом переопределить операции % и ÷ 
#(по аналогии с соответствующими целочисленными операциями) для пользовательского типа Polynom.

function %(a::Poly,b::Poly)
    none,res = divrem(a.coeff,b.coeff)
    return res
end

function ÷(a::Poly,b::Poly)
    res,none = divrem(a.coeff,b.coeff)
    return res
end

#Задача 7. Для нашего пользовательского типа Polynom пределить ещё две функции, diff и integral, 
#реализующие опрерации вифференцирования и интегрирования, соответственно.

function diff(p)
    res = zeros(length(p.coeff)-1)
    for i in 1:length(res)
        res[i]=(length(p.coeff)-i)*p.coeff[i]
    end
    return Poly{Real}(res)
end

function integral(p)
    res = zeros(length(p.coeff))
    for i in 1:length(res)
        res[i]=p.coeff[i]/(length(p.coeff)-i+1)
    end
    return Poly{Real}(res)
end

#Задача 8. Пусть задан массив данных a=[a[1],...,a[N]]. Написать функцию currenstd, получающую на вход некоторую последовательность, 
#и возвращающую массив оценок "текущего" значения стандартного отклонения std (std=sqrt(D)), получаемых для по первым n 
#членам последовательности a (n=1,2,...N). Данный алгоритм должен быть онопроходным.

function currentstd(series)
    S¹ = eltype(series)(0)
    S² = eltype(series)(0)
    D=0
    M=0
    std = zeros(length(series))
    for (n,a) in enumerate(series)
        S¹ += a
        S² += a^2
        M = S¹/n
        D = S²/n-M^2
        std[n] = sqrt(D)
    end
    return std
end

#Задача 9. Построить соответствующее индуктивное расширение реализовать однопроходный алгоритм, 
#вычисляющий наибольшее значение обобщенной частичной суммы числовой последовательности.

function sub_sum(series)
    elem = series[1]
	sum = 0
	max_sum = 0
for i in 1:length(series)
	sum += series[i]
	elem = min(elem, sum - max_sum)
	max_sum = max(max_sum, sum)
end
return max_sum
end

#Задача 10. Реализовать однопроходный алгоритм, возвращающий диапазон индексов элементов заданной числовой последовательности, 
#соответствующего наибольшему значению обобщенной частичной суммы числовой последовательности.

function index_max_sub_sum(series)
    elem = series[1]
	elem_l = 1
	elem_r = 1
	sum = 0
	max_sum = 0
	max_pos = 0
for i in 1:length(series)
	sum += series[i]
	cur = sum - max_sum
	if (cur > elem)
		elem = cur
		elem_l = max_pos + 1
		elem_r = i
    end
	if (sum < max_sum)
		max_sum = sum
		max_pos = i
    end
end
return (elem_l,elem_r)
end