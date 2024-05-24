function get_dp3(z, H, μ)
    q₁, q₂, p₁, p₂ = z
    r₁³= ((q₁-μ)^2   + q₂^2)^1.5 # distance to m1, LARGER MASS
    r₂³= ((q₁+1-μ)^2 + q₂^2)^1.5 # distance to m2, smaller mass
   
    ∂p₃_∂q₁ = 2( p₂ - (1-μ)*(q₁-μ)/r₁³ - μ*(q₁+1-μ)/r₂³)
    ∂p₃_∂q₂ = 2(-p₁ - (1-μ)*    q₂/r₁³ - μ*      q₂/r₂³)
    ∂p₃_∂p₁ = -2q₂ - 2p₁
    ∂p₃_∂p₂ =  2q₁ - 2p₂
    Dp₃ = [∂p₃_∂q₁, ∂p₃_∂q₂, ∂p₃_∂p₁, ∂p₃_∂p₂]

    if q₂ < 0 # sign convention for p₃, since the √ operator is always positive
        if q₁ < (μ-1) # L₂
            Dp₃ = -Dp₃
        end
    elseif q₂ > 0
        if q₁ > (μ-1) # L₁
            Dp₃ = -Dp₃
        end
    end

    return Dp₃
end

function get_p3(x, H, μ)
    H = params.H
    μ = params.μ
    q₁, q₂, q₃, p₁, p₂ = x
    r₁ = sqrt((q₁-μ)^2 + q₂^2 + q₃^2)
    r₂ = sqrt((q₁+1-μ)^2 + q₂^2 + q₃^2)
    p₃ = √(2(H - q₂*p₁ + q₁*p₂ + (1-μ)/r₁ + μ/r₂) - p₁^2 - p₂^2)
    if q₂ < 0 # sign convention for p₃, since the √ operator is always positive
        if q₁ < (μ-1) # L₂
            p₃ = -p₃
        end
    elseif q₂ > 0
        if q₁ > (μ-1) # L₁
            p₃ = -p₃
        end
    end
    
    x[6] = p₃
end

function γ(z, H, μ) # changes z∈ℝ⁴  into x∈ℝ⁶
    x = zeros(6)
    x[[1,2,4,5]] = z # note q₃ = 0
    get_p3(x, H, μ) # update x[6] = p₃ from the Hamiltonian
end

Dγ(z, H, μ) = [1    0    0    0;
                0    1    0    0;
                0    0    0    0;
                0    0    1    0;
                0    0    0    1;
                get_dp3(z, H, μ)']

function γ⁻¹(x) # changes x∈ℝ⁶ into z∈ℝ⁴
    z = x[[1,2,4,5]]
end

Dγ⁻¹(x) = [1 0 0 0 0 0; 
            0 1 0 0 0 0; 
            0 0 0 1 0 0; 
            0 0 0 0 1 0]

