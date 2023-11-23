function firstpassagetime(domain, PATH, args...)
    a, b = domain
    counter = 1
    while true
        t, path = PATH(2^counter, args...)
        index = findfirst(path) do x
            x<=a || x>= b
        end
        !isnothing(index) && return t[index]
        counter += 1
        counter == 60 && error("NOT PASS!")
    end
end

function occupationtime(domain, PATH, T, args...)
    return 0
end