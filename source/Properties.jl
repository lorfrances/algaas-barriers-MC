function BinID(Posx)
    # assigns particle to its bin 1D x-direction (assumes uniform grid)
    RawBin = floor(Int, Posx/NodeLen) + 1
    Bin = clamp(RawBin, 1, NumBin)
    return Bin
end

function InstantaneousProperties(ParPos, ParStat, ParWvk, ParValley, ParPhysValley, WeightedCharge, Dist, Mat, LocNode, BinLen, BarHeight, Elip)
    WeightedCharge .= 0.0
    Dist.Bin .= 0
    Dist.Vel .= 0.0
    Dist.Ene .= 0.0
    Dist.ParEne .= 0.0
    Dist.ParVel .= 0.0
    Dist.Valley .= 0
    Dist.Current .= 0.0
    Dist.CurrentSpectrum .= 0.0
    for i in 1:NumPar
        if ParStat[i] == 0 # active
            Valley = ParValley[i]
            PhysValley = ParPhysValley[i]
            Posx = ParPos[i,1]
            Bin = BinID(Posx)
            kx = ParWvk[i,1]
            ky = ParWvk[i,2]
            kz = ParWvk[i,3]
            @views MIn = Elip.InvM[Bin, PhysValley, :]
            @views Tr = Elip.TrF[Bin, PhysValley, :]
            alpha = Mat.NonParabolic[Valley, Bin]/e_c # non-parabolic parameter (constant in space)
            w_x = kx*Tr[1] + ky*Tr[4] + kz*Tr[7]
            w_y = kx*Tr[2] + ky*Tr[5] + kz*Tr[8]
            w_z = kx*Tr[3] + ky*Tr[6] + kz*Tr[9]       
            if alpha > 1e-6
                Ene = 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0)))
            else
                Ene = hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0)
            end
            Ene_eV = Ene/e_c
            Velx = hbar/(m0*(1 + 2*alpha*Ene)) * (kx*MIn[1] + ky*MIn[4] + kz*MIn[7])  
            Dist.Bin[Bin] = Dist.Bin[Bin] + 1
            Dist.Vel[Bin] = Dist.Vel[Bin] + Velx
            Dist.Ene[Bin] = Dist.Ene[Bin] + Ene_eV
            Dist.Valley[Valley, Bin] = Dist.Valley[Valley, Bin] + 1

            # create various energy indices
            CBMValley = argmin(Mat.Eg[:,Bin])   # local global valley minimum (probably gamma valley)
	        CBMValley_o = argmin(Mat.Eg[:,1]) 
            AddValley = (Mat.Eg[Valley,Bin]-Mat.Eg[CBMValley,Bin]) # valley seperation from minimuum
	        AddValley_o = (Mat.Eg[Valley, 1]-Mat.Eg[CBMValley_o, 1]) # valley seperation from minimuum, in bulk

            EnInx = round(Int,Ene_eV/EnergyStepProp) # kinetic energy compoennt only
            EnInx_valley = round(Int,(Ene_eV+AddValley)/EnergyStepProp) # valley energy and kinetic energy
            
            EnInx = clamp(EnInx, 1, NumEnergyLevelProp) # clamp to data 
            EnInx_valley = clamp(EnInx_valley, 1, NumEnergyLevelProp)

            VelInx = round(Int, Velx/VelStep) + NumVelLevel + 1
            VelInx = clamp(VelInx, 1, NumVelLevel*2 + 1)

            if (Posx >= HBPosL) && (Posx <= HBPosR)
                AddBar = BarHeight[Valley] + AddValley_o + Dist.NodePotential[Bin]
            else
                AddBar = 0.0 + AddValley_o + Dist.NodePotential[Bin]
            end
            AddShift = round(Int, (Ene_eV + AddBar)/EnergyStepProp)
            AddShift = clamp(AddShift , 1, NumEnergyLevelProp)

            Dist.ParEne[Valley, EnInx, Bin] = Dist.ParEne[Valley, EnInx, Bin] + 1
            Dist.ParEne[4, EnInx_valley, Bin] = Dist.ParEne[4, EnInx_valley, Bin] + 1
            Dist.ParVel[VelInx, Bin] = Dist.ParVel[VelInx, Bin] + 1
            Dist.CurrentSpectrum[EnInx_valley, Bin] = Dist.CurrentSpectrum[EnInx_valley, Bin] + Velx
            if PoiScheme == 2
                if Posx > LocNode[Bin] # right
                    NextNode = Bin + 1 # the next nearest node is the next one on the right
                    if NextNode > NumBin
                        NextNode = Bin # if we are at the last section then there is no next nearest bin
                    end
                elseif Posx < LocNode[Bin]
                    NextNode = Bin - 1
                    if NextNode < 1
                        NextNode = 1
                    end
                else
                    NextNode = Bin
                end
                WeightN = 1 - (abs(Posx - LocNode[Bin])/NodeLen)
                WeightNN = 1 - WeightN
                WeightedCharge[Bin] = WeightedCharge[Bin] + WeightN*ParCharge
                WeightedCharge[NextNode] = WeightedCharge[NextNode] + WeightNN*ParCharge
            elseif PoiScheme == 1
                WeightedCharge[Bin] = WeightedCharge[Bin] + 1*ParCharge
            end
        end
    end
    Dist.Current .= (ParCharge./BinLen) .* Dist.Vel .* e_c
    Dist.CurrentSpectrum .= Dist.CurrentSpectrum .* e_c .* (ParCharge./BinLen)'

    return Dist, WeightedCharge
end

function UpdateAverages(Dist, Totals, step)
    Totals.CurrentBin .= Totals.CurrentBin .+ Dist.Current
    Totals.BinCount .= Totals.BinCount .+ Dist.Bin
    Totals.VelDist .= Totals.VelDist .+ Dist.Vel
    Totals.EneDist .= Totals.EneDist .+ Dist.Ene
    Totals.NodeField .= Totals.NodeField .+ Dist.NodeField
    Totals.NodePotential .= Totals.NodePotential .+ Dist.NodePotential
    Totals.ValleyDist .= Totals.ValleyDist .+ Dist.Valley

    if step >= stepStatisticsBegin
        Totals.ParEneDist .= Totals.ParEneDist .+ Dist.ParEne
        Totals.ParVelDist .= Totals.ParVelDist .+ Dist.ParVel
        Totals.CurrentSpectra .= Totals.CurrentSpectra .+ Dist.CurrentSpectrum
    end
    return Totals
end

function WriteToFileMain(Totals, OutFile, step, Mat, Dist, BinLen, BarHeight, LocNode)
    # writes averages to file and resets them
    AveBinDist = Totals.BinCount ./ FreqWriteOut
    AveVelDist = Totals.VelDist ./ FreqWriteOut ./ AveBinDist
    AveEneDist = Totals.EneDist ./ FreqWriteOut ./ AveBinDist
    AveCurrDist = Totals.CurrentBin ./ FreqWriteOut

    AveGDist = Totals.ValleyDist[1,:] ./ FreqWriteOut
    AveLDist = Totals.ValleyDist[2,:] ./ FreqWriteOut
    AveXDist = Totals.ValleyDist[3,:] ./ FreqWriteOut
    AveField = Totals.NodeField ./ FreqWriteOut
    AvePotential = Totals.NodePotential ./ FreqWriteOut

    # perfrom the phonon decomposition
    GNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[2,:] - Dist.ScattCountBin[1,:]).*e_c.*Mat.EpO)./(WriteOutTime)
    GLNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[7,:] - Dist.ScattCountBin[6,:]).*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
    GXNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[9,:] - Dist.ScattCountBin[8,:]).*e_c.*Mat.EpIV[2,:])./(WriteOutTime)
    LNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[11,:] - Dist.ScattCountBin[10,:]).*e_c.*Mat.EpO)./(WriteOutTime)
    LGNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[16,:] - Dist.ScattCountBin[15,:]).*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
    LXNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[18,:] - Dist.ScattCountBin[17,:]).*e_c.*Mat.EpIV[3,:])./(WriteOutTime)
    LLNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[20,:] - Dist.ScattCountBin[19,:]).*e_c.*Mat.EpIV[4,:])./(WriteOutTime)
    XNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[22,:] - Dist.ScattCountBin[21,:]).*e_c.*Mat.EpO)./(WriteOutTime)
    XGNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[27,:] - Dist.ScattCountBin[26,:]).*e_c.*Mat.EpIV[2,:])./(WriteOutTime)
    XLNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[29,:] - Dist.ScattCountBin[28,:]).*e_c.*Mat.EpIV[3,:])./(WriteOutTime)
    XXNetPhononEmissionBin = ((ParCharge./BinLen).*(Dist.ScattCountBin[31,:] - Dist.ScattCountBin[30,:]).*e_c.*Mat.EpIV[5,:])./(WriteOutTime)

    GLAbsorbTot = ((ParCharge./BinLen).*Dist.ScattCountBin[6,:].*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
    GLEmitTot = ((ParCharge./BinLen).*Dist.ScattCountBin[7,:].*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
    LGAbsorbTot  = ((ParCharge./BinLen).*Dist.ScattCountBin[13,:].*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
    LGEmitTot = ((ParCharge./BinLen).*Dist.ScattCountBin[14,:].*e_c.*Mat.EpIV[1,:])./(WriteOutTime)
   
    open(OutFile.BinCountAve,"a") do io
        writedlm(io, (AveBinDist)', ' ')
    end

    open(OutFile.VelDistAve,"a") do io
        writedlm(io, (AveVelDist)', ' ')
    end   

    open(OutFile.CurrentDataBinAve,"a") do io
        writedlm(io, (AveCurrDist)', ' ')
    end 

    open(OutFile.EnDistAve,"a") do io
        writedlm(io, (AveEneDist)', ' ')
    end  

    open(OutFile.GValleyData,"a") do io
        writedlm(io, (AveGDist)', ' ')
    end  

     open(OutFile.LValleyData,"a") do io
        writedlm(io, (AveLDist)', ' ')
    end  

    open(OutFile.XValleyData,"a") do io
        writedlm(io, (AveXDist)', ' ')
    end 

    open(OutFile.FieldTimeData,"a") do io
        writedlm(io, (AveField)', ' ')
    end 

    open(OutFile.PotentialTimeData,"a") do io
        writedlm(io, (AvePotential)', ' ')
    end  

    # phonon decomposition
    open(OutFile.GPhononCountData, "a") do io
        writedlm(io, (GNetPhononEmissionBin)', ' ')
    end

    open(OutFile.GLPhononCountData, "a") do io
        writedlm(io, (GLNetPhononEmissionBin)', ' ')
    end
    open(OutFile.GXPhononCountData, "a") do io
        writedlm(io, (GXNetPhononEmissionBin)', ' ')
    end
    open(OutFile.LGPhononCountData, "a") do io
         writedlm(io, (LGNetPhononEmissionBin)', ' ')
    end
    open(OutFile.LXPhononCountData, "a") do io
         writedlm(io, (LXNetPhononEmissionBin)', ' ')
    end
    open(OutFile.LLPhononCountData, "a") do io
        writedlm(io, (LLNetPhononEmissionBin)', ' ')
    end
    open(OutFile.XGPhononCountData, "a") do io
        writedlm(io, (XGNetPhononEmissionBin)', ' ')
    end
    open(OutFile.XLPhononCountData,"a") do io
        writedlm(io, (XLNetPhononEmissionBin)', ' ')
    end
    open(OutFile.XXPhononCountData,"a") do io
        writedlm(io, (XXNetPhononEmissionBin)', ' ')
    end
    open(OutFile.LPhononCountData,"a") do io
        writedlm(io, (LNetPhononEmissionBin)', ' ')
    end
    open(OutFile.XPhononCountData,"a") do io
        writedlm(io, (XNetPhononEmissionBin)', ' ')
    end
    open(OutFile.GLTotAbsorbData,"a") do io
        writedlm(io, (GLAbsorbTot)', ' ')
    end
    open(OutFile.GLTotEmitData,"a") do io
        writedlm(io, (GLEmitTot)', ' ')
    end
    open(OutFile.LGTotAbsorbData,"a") do io
        writedlm(io, (LGAbsorbTot)', ' ')
    end
    open(OutFile.LGTotEmitData,"a") do io
        writedlm(io, (LGEmitTot)', ' ')
    end

    BandProfile(AvePotential, LocNode, BarHeight, Mat)
    #kPtsValley(ParWvk, ParStat, ParPhysValley, Elip)
    Totals.BinCount .= 0
    Totals.VelDist .= 0.0
    Totals.EneDist .= 0.0
    Totals.ValleyDist .= 0
    Totals.NodeField .= 0.0
    Totals.NodePotential .= 0.0
    Totals.CurrentBin .= 0.0
    Dist.ScattCountBin .= 0
    if step > stepStatisticsBegin
        AveParEneDist = Totals.ParEneDist[4,:,:] ./ (step - stepStatisticsBegin)
        AveParEneDistG = Totals.ParEneDist[1,:,:] ./ (step - stepStatisticsBegin)
        AveParEneDistL = Totals.ParEneDist[2,:,:] ./ (step - stepStatisticsBegin)
        AveParEneDistX = Totals.ParEneDist[3,:,:] ./ (step - stepStatisticsBegin)  
        AveParVelDist = Totals.ParVelDist ./ (step - stepStatisticsBegin)
        AveCurrentSpectra = Totals.CurrentSpectra ./ (step-stepStatisticsBegin) ./ 1e6

        writedlm(OutFile.ParEneDistData, AveParEneDist)
        writedlm(OutFile.ParEneDistDataG, AveParEneDistG)
        writedlm(OutFile.ParEneDistDataL, AveParEneDistL)
        writedlm(OutFile.ParEneDistDataX, AveParEneDistX)
        writedlm(OutFile.ParVelDistData, AveParVelDist)
        writedlm(OutFile.CurrentSpectra, AveCurrentSpectra)
    end

    return Totals, Dist
end



function BandProfile(AveNodePotential, LocNode, BarHeight, Mat)
    BandDiag = "BandDiag" * suff * ".dat"
    # RETURNS VISUALIZATION OF BANDPROFILE, HAS ITS OWN FREQUENCY if desired
    BandDiagDat = zeros(NumBin,4) 
    Bands = zeros(NumBin)
    CBMValley = argmin(Mat.Eg[:,1])

    for i in 1:NumBin
        Pot_i = 0.0
        SlfCons = AveNodePotential[i] 
        for j in 1:NumValley
            if LocNode[i] >= HBPosL && LocNode[i] < HBPosR && HB == 1
                Pot_i = SlfCons + BarHeight[j] + (Mat.Eg[j,1]-Mat.Eg[CBMValley,1])

            else
                Pot_i = SlfCons + (Mat.Eg[j,1]-Mat.Eg[CBMValley,1])
            end
            BandDiagDat[i,1+j] = Pot_i 
        end
    end
    BandDiagDat[:,1] = LocNode
    writedlm(BandDiag, BandDiagDat)
end   

function kPtsValley(ParWvk, ParStat, ParPhysValley, Elip)
    filename = "kPtsValley" * suff * ".dat"
    open(filename, "w") do io
        for i in 1:NumPar
            if ParStat[i] == 0 # active
                if ParPhysValley[i] == 1
                    println(io, ParWvk[i,1], " ", ParWvk[i,2], " ", ParWvk[i,3])
                end
                if ParPhysValley[i] == 2
                    L_point = pi/(5.65e-10)
                    l_max = sqrt(3)*pi/(5.65e-10)
                    x = ParWvk[i,1] + L_point
                    y = ParWvk[i,2] + L_point
                    z = ParWvk[i,3] + L_point
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 3
                    L_point = pi/(5.65e-10)
                    l_max = sqrt(3)*pi/(5.65e-10)
                    x = ParWvk[i,1] - L_point
                    y = ParWvk[i,2] + L_point
                    z = ParWvk[i,3] + L_point
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 4
                    L_point = pi/(5.65e-10)
                    l_max = sqrt(3)*pi/(5.65e-10)
                    x = ParWvk[i,1] + L_point
                    y = ParWvk[i,2] - L_point
                    z = ParWvk[i,3] + L_point
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 5
                    L_point = pi/(5.65e-10)
                    l_max = sqrt(3)*pi/(5.65e-10)
                    x = ParWvk[i,1] + L_point
                    y = ParWvk[i,2] + L_point
                    z = ParWvk[i,3] - L_point
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 6
                    X_point = 2*pi/(5.65e-10)
                    l_max = 2*pi/(5.65e-10)
                    x = ParWvk[i,1] + X_point
                    y = ParWvk[i,2]
                    z = ParWvk[i,3]
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 7
                    X_point = 2*pi/(5.65e-10)
                    l_max = 2*pi/(5.65e-10)
                    x = ParWvk[i,1]
                    y = ParWvk[i,2] + X_point
                    z = ParWvk[i,3]
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
                if ParPhysValley[i] == 8
                    X_point = 2*pi/(5.65e-10)
                    l_max = 2*pi/(5.65e-10)
                    x = ParWvk[i,1] 
                    y = ParWvk[i,2]
                    z = ParWvk[i,3] + X_point
                    @views RoT = Elip.Rot[1, ParPhysValley[i], :]
                    l = x*RoT[1] + y*RoT[4] + z*RoT[7]
                    t1 = x*RoT[2] + y*RoT[5] + z*RoT[8]
                    t2 = x*RoT[3] + y*RoT[6] + z*RoT[9]
                    if l >= l_max
                        l_new = -2*l_max + l
                        x = l_new*RoT[1] + t1*RoT[2] + t2*RoT[3]
                        y = l_new*RoT[4] + t1*RoT[5] + t2*RoT[6]
                        z = l_new*RoT[7] + t1*RoT[8] + t2*RoT[9]
                    end
                    println(io, x, " ", y, " ", z)
                end
            end
        end
    end
end