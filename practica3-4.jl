#Задача 1. Написать функцию с заголовком findallmax(A::AbstractVector)::AbstractVector{Int}, 
#возвращающую вектор индексов всех элементов массива A, имеющих максимальное значение. Алгоритм должен быть однопроходным, 
#т.е. иметь асимптотическую оценку вычислительной сложности O(n).

function find_all_max(a)
    i_max=[firstindex(a)]
    for i in firstindex(a)+1:lastindex(a)
        if A[i]>A[i_max[end]] # # Если вдруг нашли число большее нашего текущего максимума
            i_max=[i]
        elseif A[i]==A[i_max[end]] # Если нашли еще один максимум
            push!(i_max, i)
        end
    end
    return i_max #возвращающую вектор индексов всех элементов массива
end

#Задачи 2-3. Реализовать эти две разновидности пузырьковой сортировки.

function shenker!(a)
    len = length(a) # длина массива  
    flag = true     #создаем флаг
    first_index = 1 # перый индекс массива
    end_index = len - 1 # последний индекс массива

    while (flag == true)
        flag = false  

        # проход сначала 
        for i in first_index : end_index
            if (a[i] > a[i + 1])                  # Если выбранные элементы стоят не по возрастанию
                a[i] , a[i + 1] = a[i + 1] , a[i] # то меняем их местами
                flag = true
            end
        end

        if (not(flag)) # Если flag = false (то есть массив изначально был отсортирован)
            break      # то выходим из цикла
        end

        flag = false

        end_index = end_index - 1 

        # проход c конца
        for i in (end_index - 1):-1:(first_index - 1)
            if (a[i] > a[i + 1])                         # Если выбранные элементы стоят не по возрастанию
                a[i] , a[i + 1] = a[i + 1] , a[i]       # то меняем их местами
                flag = true
            end
        end

        first_index+= 1
    end

    return a # возвращаем видоизмененный массив
end

# Задача 4. Реализовать сортировку Шелла и проверить работоспособность алгоритма при следующих наборах значений промежутков
function shellsort!(a)
    n=length(a)
    distseries=(n÷2^i for i in 1:Int(floor(log2(n)))) 
    # distseries - это ГЕНЕРАТОР последовательности промежутков (т.е. члены этой последовательности будут вычисляться в процессе выполнения следующего цикла, заранее в памяти они не размещаются, как было бы при использовании массива или кортежа вместо генератора)
    for d in distseries
        for i in firstindex(a):lastindex(a)-d
            j=i
            while j>=firstindex(a) && a[j]>a[j+d]
                a[j],a[j+d] = a[j+d],a[j]
                j-=d
            end
        end
    end
    return a
end

# Задача 5. Написать функцию с заголовком => slice(A::AbstractVector, p::AbstractVector{<:Integer})
function slice(A::Vector{T},p::Vector{Int}) :: Vector{T} where T
    b :: Vector{T}
    s = 0
    for i in p
        b[s] = a[i]
        s += 1
    end
    return b
end

# Задча 6. Пусть perm - это некоторый вектор перестановок индексов одномерного массива A. 
# Написать свою реализацию встроенной функции permute!(A, perm), реализующую соответствующее премещение элементов массива A
# на месте (in-plice), т.е. без копирования их в новый массив . (Cвой вариант этой функции можно назвать permute_!).

function permute_!(A::Vector{T},perm::Vector{Int})::Vector{T} where T
    return slice(A,perm)
end

#Задача 7. Реализовать встроенные функции вставки/добавления (deleteat!, insert!) элемента массива

function deleteat!(b,arr) # удаление элементов
    flag = false
    for  i in 0:count-1
        if (arr[i] == value)
            index = i
            flag = true
        else
            index = -1
    end
    count = length(b)
    for  i in index:count-1
        b[i] = b[i + 1]
        count -= 1
    end
end

function insert!(b,x,index) # вставка элемента по индексу
    for  i in (count - 1):-1:index
        b[i + 1] = b[i]
    end
    b[index] = x
    count += 1
end

# Задача 8. Реализовать встроенные функции
#Для функции unique! обеспечить асимптотическую оценку вычислительной сложности O(n*log(n)).

function unique!(a) #unique! (удаляет из исходного массива повторяющиеся элементы, оставляя каждый элемент в единственном экземпляре), 
    a = sort!(a)                #сортируем массив
    for i in 0:length(a) - 2
        if a[i]==a[i+1]             
            deleteat!(a,a[i+1]) #удаляем повторяющиеся элементы
        end
    end
    return a
end


function allunique(a) #allunique (проверяет, состоит ли данный массив только из уникальных элементов). 
    a = sort!(a)                 #сортируем массив
    for i in 0:length(a) - 2
        if a[i]==a[i+1]
            return false # массив состоит не только из уникальных элементов
        end
    end
    return true # массив состоит только из уникальных элементов
end


# Задача 9. Реализовать встроенную функцию reverse!, переставляющую элементы в обратном порядке в самом массиве, т.е. "на месте"

function reverse!(a)
    len = length(a)
    n = len
    for i in 0:len/2
        a[i] , a[i+n] = a[i+n], a[i] # меняем элементы местами, начиная с середины массива
        n -= 2
    end
end

# Задача 10. Написать функцию, осуществляющую циклический сдвиг массива на m позиций "на месте", т.е. без использования дополнительного массива.

function shift(a,m)
    count = length(a) + m
    for i in (count - n) : -1 : 0
        a[mod((i + n),n)] = a[mod(i,n)]
    end
end

# Задача 11. Реализовать функцию, аналогичную встроенной функции transpose, с использованием вспомогательного массива

function transpose!(a::Matrix) # с использование вспомогательного массива 
    x = length(a) # строки 
    y = length(a[1]) # столбцы
    b :: Matrix
    for i in 1:a
        for j in 1:y 
            b[j][i] = a[i][j] # транспонируем матрицу
        end
    end
    return b # возвращаем новую матрицу равную транспонированной матрицы a
end
# Задача 12. Реализовать функцию, аналогичную встроенной функции transpose, осуществляющую транспонирование матрицы "на месте" (без использования вспомогательного массива)

function transpose!(a::Matrix) # без использование вспомогательного массива
    x = length(a) # строки 
    y = length(a[1]) # столбцы
    b :: Matrix
    for i in 1 : a/2+1      # транспонируем исходную матрицу без использование вспомогательного массива
        for j in 1:y/2+1    # а используя сортировку: "пузырек"
            t = a[j][i]  
            a[j][i] = a[i][j]
            a[i][j] = t
        end
    end
    return a  # возвращаемтранспонированной матрицу
end