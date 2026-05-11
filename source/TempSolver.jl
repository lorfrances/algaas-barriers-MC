function TempProps(Mat, LocNode)
    # this function prints the simulation results required for evaluating the temperature

    # read source
    pDataG  = readdlm("GPhononCountData.dat",  skipstart=20)
    pDataGL = readdlm("GLPhononCountData.dat", skipstart=20)
    pDataGX = readdlm("GXPhononCountData.dat", skipstart=20)
    pDataLG = readdlm("LGPhononCountData.dat", skipstart=20)
    pDataLX = readdlm("LXPhononCountData.dat", skipstart=20)
    pDataLL = readdlm("LLPhononCountData.dat", skipstart=20)
    pDataXG = readdlm("XGPhononCountData.dat", skipstart=20)
    pDataXL = readdlm("XLPhononCountData.dat", skipstart=20)
    pDataXX = readdlm("XXPhononCountData.dat", skipstart=20)
    sdot_ave = vec(mean(pDataGL, dims=1) .+ mean(pDataLG, dims=1) .+ mean(pDataG,dims=1) .+
              mean(pDataGX, dims=1) .+ mean(pDataXG, dims=1) .+
              mean(pDataLX, dims=1) .+ mean(pDataXL, dims=1) .+
              mean(pDataLL, dims=1) .+ mean(pDataXX, dims=1))
    data = zeros(NumBin,3)
    LatTemp = "TempParam" *suff *".dat"
    data[:,1] = LocNode
    data[:,2] = Mat.ThermCond
    data[:,3] = sdot_ave
    writedlm(LatTemp, data)
end