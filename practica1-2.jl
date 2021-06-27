#Задача 1: Реализовать функцию, аналогичную встроенной функции reverse!, 
#назвав её, например, reverse_user!, для следующих случаев: 

# a) аргумент функции - вектор
function reverse_user!(array) 
    n = length(array)
    array = array[n:-1:1] 
end

function reverse_user(array)
    return reverse_user!(deepcopy(array))
end

# б) аргумент функции - матрица (2-мерный массив)
function reverse_user!(matrix, dim)
    for i = 1 : dim
        matrix[i, :] = reverse_user!(matrix[i, :])  #i-ая строка
    end
end

function reverse_user(matrix, dim)
    return reverse_user!(deepcopy(matrix), dim)
end

#Задача 2: Аналогично, реализовать функцию, аналогичную встроенной функции copy для следующих случаев: 
# a) аргумент функции - вектор

function copy_user!(array)
    New_vector_zero = zero(length(array)) #Обнуляем вектор
    for i in 1:length(array)
        New_vector_zero[i] = array[i] #копируем значения в новый вектор
    end
    return New_vector_zero #возвращаем новый вектор
end
# б) аргумент функции - матрица (2-мерный массив)

function copy_user!(matrix, dim)
    New_Matrix_zero = zero(matrix[i, :]) #Обнуляем матрицу
    for i in 1: dim
        New_Matrix_zero[i] = matrix[i, :] #копируем значения в новую матрицу
    end
    return New_Matrix_zero # возвращаем новую матрицу
end

#Задача 3 Реализовать алгоритм сортировки методом пузырька, написав следующие 4 обобщенные функции: bubblesort, bubblesort!, bubblesortperm, bubblesortperm!, по аналогии со встоенными функциями
# sort!, sort, sortperm!, sortperm, ограничившись только случаем, когда входной параметр есть одномерный массив (вектор).

#sort(x) - х не меняет, возвращает копию х отсортированную
function bubblesort!(a)
    n = length(a)

    for i = 1:n-1
        for j = 1:n-i
            if a[j]>a[j+1]
                a[j],a[j+1]=a[j+1],a[j] #меняем соседние элементы местами
            end
        end
    end
    return a #возвращаем видоизмененный вектор
end

bubblesort(a)=bubblesort!(deepcopy(a))

function bubblesort(a) #sort!(x) - отсортированный х, он сам меняется
    n = length(a)
    a=deepcopy(a)
    for i = 1:n-1
        for j = 1:n-i
            if a[j]>a[j+1]
                a[j],a[j+1]=a[j+1],a[j] #меняем соседние элементы местами
            end
        end
    end
    return a #возвращаем видоизмененный вектор
end

function bubblesortperm!(a) #sortperm!(x) - вектор индексов перестановки х
    n = length(a)
    b = collect(1:length(a))  #массив индексов
    for i = 1:n-1
        for j = 1:n-i
            if a[j]>a[j+1]
                a[j],a[j+1]=a[j+1],a[j] #меняем соседние элементы местами  
                b[j],b[j+1]=b[j+1],b[j] #меняем соседние индексы местами
            end
        end
    end
    return b #возвращаем видоизмененный массив индексов
end

bubblesortperm(a)=bubblesortperm!(deepcopy(a))

function bubblesortperm(a) #sortperm(x) - массив сам не меняется, только копия
    n = length(a)

    b = collect(1:length(a)) #массив индексов
    b = deepcopy(b)

    for i = 1:n-1
        for j = 1:n-i
            if a[j]>a[j+1]
                a[j],a[j+1]=a[j+1],a[j] #меняем соседние элементы местами  
                b[j],b[j+1]=b[j+1],b[j]  #меняем соседние индексы местами
            end
        end
    end
    return b #возвращаем видоизмененную копию массива индексов
end

#Задача 4: На основе разработанных в пункте 1 функций, сотрирующих одномерный массив, написать соответствующие функции, 
#которые бы могли получать на вход матрицу, и сортировать каждый из ее столбцов по отдельности. 
#Имена функций оставить прежними, что были и в пункте 1, воспользовавшись механизмом множественной диспетчеризации языка Julia.

    function bubblesort!(M::Matrix)::Matrix
        array = view(M,:,1)
        n = length(array)
        for i in size(i,n-1) 
            bubblesort!(view(M,:,i)) #сортируем выбранный столбец
        end
        return M  #возвращаем матрицу
    end

    bubblesort(M)=bubblesort!(deepcopy(M))

    function bubblesortperm!(M::Matrix)::Matrix
        index = Matrix{double}(undef,size(M))  #создание матрицы типа double
        array = view(M,:,1)
        n = length(array)
        for i in size(A,n-1)
            index[:,i] = bubblesortperm!(view(M,:,i)) 
        end
        return index
    end
    
    bubblesortperm(M)=bubblesortperm!(deepcopy(M))

# Задача 5 : отсортировать столбцы матрицы в порядке возрастания их суммы

function matrix_sort_a!(a::Matrix)
    x = view(a,:,1)
    n = length(x)

    for i = 1:n-1
        for  j = 1:n-i
            B = view(a,:,j) # делаем срез по выбранному столбцу
            Bsum = sum(B) # суммируем значения в выбранном столбце
            C = view(a,:,j+1) # делаем срез по выбранному столбцу
            Csum = sum(C) # суммируем значения в выбранном столбце
            if Bsum > Csum # Если общая cумма левого столбца больше общей суммы правого столбца
                B,C=C,B     # то меняем данным столбцы местами
            end
        end
    end
    return a # возвращаем видоизмененную матрицу
end

# Задача 6 : отсортировать столбцы матрицы в порядке возрастания кол-ва нулей в них
function matrix_sort_b!(a::Matrix)
    x = view(a,:,1)
    n = length(x)

    for i = 1:n-1
        for j = 1:n-i
            B = view(a,:,j)                 # делаем срез по выбранному столбцу
            len = length(findall((x==0),B))  # считаем кол-во нулей в выбранном столбце
            C = view(a,:,j+1)               # делаем срез по выбранному столбцу
            m = length(findall((x==0),C)) # считаем кол-во нулей в выбранном столбце
            if  len > m  # Если общее кол-во нулей в левом столбце больше общего кол-ва нулей в правом столбце
                B,C=C,B  # то меняем данным столбцы местами
            end
        end
    end
    return a  # возвращаем видоизмененную матрицу
end

#Задача 7 Написать функцию sortkey!(a, key_values), получающую на вход некоторый вектор a, и соответствующий вектор keyvalues ключевых значений элементов вектора a, осуществляющую сортировку вектора a по ключевым значениям его элементов, и возвращающую ссылку на вектор a. (Для сортировки вектора ключевых значений можно востпользоваться одной из разработанных в пункте 1 функций, или соответствующей встроенной функцией).
function sortkey!(key_values, a)
    index = sort_bubble_perm!(key_values)
    return view(a[index])
end

#Задача 8: Написать функцию calcsort, реализующую сортировку методом подсчета числа значений. 
#Рассмотреть 2 варианта функции (2 метода - в терминологии Julia): в первом варианте возможные значения элементов 
#сортируемого массива задаются некоторым диапазоном, во втором - некоторым отсортированным массивом (вектором).

function calcsort!(a, value)
    New_value_zero = zeros(Int, size(value)) #обнуляем матрицу
    for i in a
        New_value_zero[indexvalue(i,value)] += 1
    end
    s=1
    for i in eachindex(value)
        for j in 1:New_value_zero[i]
            a[k] = value[i]
            s+=1
        end
    end
    return a
end

calcsort(a, value) = calcsort!(deepcopy(a), deepcopy(value))


#Задача 10 Написать функции insertsort!, insertsort, insertsortperm, insertsortperm! (по аналогии с пунктом 1) реализующие алгоритм сортировки вставками

function insertsort!(a)
    len = length(a)
    for i in 2:len
        j = i - 1
        while j > 0 && a[j] > a[j+1]
            a[j+1],a[j] = a[j],a[j+1]
            j -= 1
        end 
    end
    return a    
end

insertsort(a) = insertsort!(deepcopy(a))

function insertsortperm!(a)
    len = length(a)
    ptr = []
    for i in 1:len
        push!(ptr, i)
    end
    for i in 2:len
        j = i - 1
        while j > 0 && a[j] > a[j+1]
            a[j+1], a[j] = a[j], a[j+1]
            ptr[j+1],ptr[j] = ptr[j],ptr[j+1]
            j -= 1
        end 
    end
    return ptr
end

insertsortperm(a) = insertsortperm!(deepcopy(a))

#Задача 11: Реализовать ранее написанную функцию insertsort! с помощью встроенной функции reduce 

insertsort!(a)=reduce(1:length(a))do _, k # при выполнении операции не используется первый аргумент
    while s>1 && a[k-1] > a[k]
        a[k-1], a[k] = a[k], a[k-1]
        s-=1
    end
    return a
end

#Задача 12: Дополнить функцию insertsort! процедурой "быстрого поиска"

function binary_search(num::Int,ptr::Array)
    left = 1
    right = length(ptr)
    while (left<=right)
        middle = div(left + right,2)
        if (num < ptr[middle])
            right = middle - 1
        elseif (num > ptr[middle])
            left = middle + 1
        else
            return middle
        end
    end
    return -1
end

#Задача 13: Написать обобщенную функцию (nummax), получающую на вход итерируемый объект, 
#содержащий некоторую последовательность (элементы которой можно сравнивать по величине), 
#и возвращающую число максимумов этой последовательности.

function nummax(a)
    len = lastindex(a)
    max = 0
    max_el= a[1]
    for i in 2:len
        if max_el == a[i] # Если нашли еще один максимум
            max += 1
        elseif max_el < a[i] # Если вдруг нашли число большее нашего текущего максимума
            max = 1
            max_el = a[i]
        end
    end
    return max # возвращаем число максимумов этой последовательности
end

#Задача 14: Написать обобщенную функцию (findallmax), получающую на вход итерируемый объект, содержащий некоторую последовательность 
#(элементы которой можно сравнивать по величине), и возвращающую вектор, составленный из индексов элементов входной последовательности, 
#имеющих максимальное значение.

function findallmax(a)
    max = maximum(a) #находим максимальное значение
    v = [] 
    for i in 1:length(a)
        if a[i] == max
            push!(v, i) 
        end
    end
    return v # возвращаем вектор индексов с макс значениями
end

#Задача 15: Написать обобщенную функцию findallmax высшего порядка, получающую на вход итерируемый объект, 
#содержащий некоторую последовательность и некоторую функцию (значение типа ::Function), 
#и возвращающую вектор, составленный из индексов элементов входной последовательности, 
#на которых заданная функция достигает максимального значения (речь идет о сужении заданной функции на заданной последовательности).

function findallmax(a, f::Function)
    b = copy(a)
    for i in 1:length(b)
        b[i] = f(b[i])
    end
    return findallmax(b)
end