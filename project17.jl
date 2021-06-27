#Задача 1. Написать и протестировать функцию, получающую на вход список смежностей некоторого графа, и возвращающую вектор индексов 
#его вершин, полученных в порядке поиска в глубину.

function dfsearch(startver::T, graph::ConnectList{T}) where T
    mark = zeros(Bool, length(graph))
    stack  = [startver]
    mark[startver] = 1
    visited = Int64[]
    while !isempty(stack)
        v = pop!(stack)
        push!(visited,v)
        for u in graph[v]
            if mark[u] == 0
                push!(stack,u)
                mark[u] = 1
            end
        end
    end
    return visited
end

#Задача 2. Написать и протестировать функцию, получающую на вход список смежностей некоторого графа, и возвращающую вектор индексов 
#его вершин, полученных в порядке поиска в ширину.

function bfsearch(startver::T, graph::ConnectList{T}) where T
    mark = zeros(Bool, length(graph))
    queue  = [startver]
    mark[startver] = 1
    visited = Int64[]
    while !isempty(queue)
        v = popfirst!(queue)
        push!(visited,v)
        for u in graph[v]
            if mark[u] == 0
                push!(queue, u)
                mark[u] = 1
            end
        end
    end
    return visited
end
# Задача 3. Написать и протестировать функцию, получающую на вход список смежностей некоторого графа, и возвращающую
#вектор валентностей его вершин по выходу.

function valence(graph::ConnectList{T}) where T
    val = zeros(size(graph,1))
    for i in 1:size(graph,1)
        val[i]=length(graph[i])
    end
    return val
end

#Задача 4. Написать и протестировать функцию, получающую на вход список смежностей некоторого графа, и возвращающую вектор 
#валентностей его вершин по входу.

function bfs_valence(graph::ConnectList{T}) where T
    mark = zeros(Int64, length(graph))
    queue  = [1]
    mark[1] = 1
    while !isempty(queue)
        v = popfirst!(queue)
        println(v)
        for u in graph[v]
            if mark[u] == 0
                push!(queue, u)
            end
            mark[u] = mark[u] + 1
        end
    end
    mark[1] -= 1
    return mark
end