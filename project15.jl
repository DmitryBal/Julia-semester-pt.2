# Задача 1. Написать функцию convert_to_nested(tree::ConnectList{T}, root::T) where T, получающую на вход дерево, 
# представленное списком смежностей tree и индексом его корня root, 
# и возвращающая представление того же дерева в виде вложенных векторов.
function convert_to_nested(tree::ConnectList{T},root::T) where T
    nested_tree = []
    for subroot in tree[root]
        push!(nested_tree, convert(tree, subroot))
    end
    push!(nested_tree, root)
    return nested_tree
end
#Задача 6 Написать функцию, получающую на вход имя некоторого типа (встоенного или пользовательского) языка Julia (тип этого аргумента - Type) и распечатывающая список всех дочерних типов в следующем формате:
function alltypes(type)
    for i in subtypes(type)
        println(i)
        alltypes(i)
    end
end