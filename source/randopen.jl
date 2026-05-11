function randopen() # random number generator preventing 0 as a output
    rn = 0.0
    while rn < 1e-40
        rn = rand()
    end
    return rn
end
