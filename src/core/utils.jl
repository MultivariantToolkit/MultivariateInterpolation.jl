# FIXME: starts with the last element, not the largest
function leja!(S::AbstractVector)
    n = length(S)
    maxi = argmax(eachindex(S))
    S[1], S[maxi] = S[maxi], S[1]

    for i in 2:n
        maxi = argmax(i:n) do j
            S_j = S[j]
            prod(S_k -> abs(S_j - S_k), view(S, 1:(i - 1)))
        end
        S[i], S[maxi] = S[maxi], S[i]
    end

    return S
end
