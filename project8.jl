#Задача 1. Реализовать рекурсивный вариант сортировки слияниями с оценками сложности и объема потребляемой памяти O(nlog(n)).

function mergesort!(a)
    if length(a) == 1
        return a
    end
    a1 = mergesort!(a[begin:end÷2])  #  begin:end делает то же самое, что и макрос @view (однако работает только в версиях Julia 1.1.5 или выше)
    a2 = mergesort!(a[end÷2+1:end])
    merge(a1, a2)   # merge() - cм. Задача 1/практика 6
end

#Задача 2. Реализовать оба варианта нерекурсивной сортировки слияниями.

function mergesort_1!(a)
    i = 1
    j = 1
    while (2*i<length(a))
        if (j>=length(a))
            j=1
            i*=2
        end
        el = j+2i-1
        if (j+2i-1>length(a))
            el = length(a)
        end
        b=merge(a[j:j+i-1], a[j+i:el])   # merge() - cм. Задача 1/практика 6
        a[j:el]=b
        j+=2*i
    end
    return a
end

# Задача 4. Реализовать алгоритм быстрой сортировки Хоара двумя способами: основываясь на двух,
# указанных выше возможных вариантах частичной сортировки (partsort!, или - partsort2!).

#Задача 5. Написать функцию, реализующую вычисление k-ой порядковой статистики, и - (на её основе) функции, вычисляющую медиану.
function print_k_statistics(ptr::Array)::Array
    mediana = None
    num = None
    if !length(ptr)
        print("Empty")
    else
        num = length(ptr)
        num_1 = (num - 1) // 2
        if num % 2 == 0
            mediana = (ptr[num_1] + ptr[num_1 + 1]) / 2 
        else
            mediana = ptr[num_1]
        end
    end
        print(f"Длина - " << num << "Среднее - " << sum(ptr) / num )
        print(f"Минимальное значение -" << min(ptr) << ", Максимальное значение - " << max(ptr))
        print(f"Медиана - " << mediana)
end