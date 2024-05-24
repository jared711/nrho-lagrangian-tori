import ThreeBodyProblem.rotz

function rv2pq(rv)
    x,y,z,ẋ,ẏ,ż = rv
    # R = [-1  0  0;
    #       0 -1  0;
    #       0  0  1] # rotate by 180 degrees (π radians) because we use the barcelona convention
    # E = [ 0 -1  0;
    #       1  0  0;
    #       0  0  0] 
    # r = rv[1:3]
    # v = rv[4:6]
    # q = R*r # rotate by 180 degrees (π radians) because we use the barcelona convention
    # p = E*R*r + R*v
    q = rotz(π)*[  x,   y, z] # rotate by 180 degrees (π radians) because we use the barcelona convention
    p = rotz(π)*[ẋ-y, ẏ+x, ż]
    return vcat(q,p)
end

function pq2rv(x)
    q₁, q₂, q₃, p₁, p₂, p₃ = x
    r = rotz(π)*[q₁, q₂, q₃]
    v = rotz(π)*[p₁ + q₂, p₂ - q₁, p₃]
    # R = [-1  0  0;
    #       0 -1  0;
    #       0  0  1] # rotate by 180 degrees (π radians) because we use the barcelona convention
    # E = [ 0 -1  0;
    #       1  0  0;    
    #       0  0  0]  
    # q = x[1:3]
    # p = x[4:6]
    # r = R*q
    # v = R*(p-E*R*r)
    return vcat(r,v)
end

function computeH(x, μ)
    q₁, q₂, q₃, p₁, p₂, p₃ = x
    r₁ = sqrt((q₁-μ)^2 + q₂^2 + q₃^2)
    r₂ = sqrt((q₁+1-μ)^2 + q₂^2 + q₃^2)
    H = 0.5*(p₁^2 + p₂^2 + p₃^2) + q₂*p₁ - q₁*p₂ - (1-μ)/r₁ - μ/r₂
end