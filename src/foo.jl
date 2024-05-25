closest_hiatus_unique = [10, 17, 35]

for h in rand(closest_hiatus_unique, length(closest_hiatus_unique))
    if model_agesₚ[h-1] == model_agesₚ[h]
        n = findclosestunequal(i)
        if n < h
            model_agesₚ[n:h-1] .= model_agesₚ[n]
        elseif n > h
            model_agesₚ[h:n] .= model_agesₚ[n]
        end
    end
end
x = 1:100
