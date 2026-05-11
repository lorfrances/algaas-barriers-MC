function main(ParPos, ParWvk, ParTau, ParValley, ParPhysValley, ParStat, LocNode, BinLen, Mat, Dist, Totals, OutFile, WeightedCharge, BarHeight, Elip, Contact, ScattTable, ScattCount, ScattIndex)

    time = 0.0
    step = 0
    CountWrite = 0

    while time < TotTime
        ParWvk, ParPos, ParValley, ParPhysValley, ParTau, ParStat, Dist = LoopParticles(ParPos, ParWvk, ParTau, ParStat, ParValley, ParPhysValley, Mat, BarHeight, Dist, LocNode, Elip, step, ScattTable, ScattCount, ScattIndex) 
        ParPos, ParWvk, ParStat, ParValley, ParPhysValley, ParTau = InjectOhmic(ParPos, ParWvk, ParStat, ParValley, ParPhysValley, ParTau, Dist, Mat, Elip, Contact) # inject fresh electrons, update ParPos, Wvk, Stat, Valley
        Dist, WeightedCharge = InstantaneousProperties(ParPos, ParStat, ParWvk, ParValley, ParPhysValley,  WeightedCharge, Dist, Mat, LocNode, BinLen, BarHeight, Elip)
        Totals = UpdateAverages(Dist, Totals, step)
        ElCharge, ElDop = ChargeDopDist(WeightedCharge, Mat, BinLen)
        Dist = PoissonSolver(ElCharge, ElDop, Dist, Mat)

        # all MC steps complete: advance counters
        time = time + Dt
        step = step + 1
        CountWrite = CountWrite + 1

        # print data occasionally
        if CountWrite == FreqWriteOut
            println("Time: ", time/1e-12, " ps")
            Totals, Dist = WriteToFileMain(Totals, OutFile, step, Mat, Dist, BinLen, BarHeight, LocNode)
            if step > stepStatisticsBegin*2
                TempProps(Mat, LocNode)
            end
            CountWrite = 0
        end


    end
end



function LoopParticles(ParPos, ParWvk, ParTau, ParStat, ParValley, ParPhysValley, Mat, BarHeight, Dist, LocNode, Elip, step, ScattTable, ScattCount, ScattIndex) 
    # loops through every particle, and drifts/scatters them in one fixed timestep in both 
    # real and kappa space - the primary MC loop

    Wvk = zeros(3)
    Pos = zeros(3)
    Vel_i = zeros(3)
    Pos_i = zeros(3)
    Wvk_i = zeros(3)
    Vel_f = zeros(3)
    EField_i = zeros(3)
    EField_f = zeros(3)
    initWvk = zeros(3)
    initPos = zeros(3)

    # above, we establish the particle property vectors for particle 1, such that we can continually recycle them for each particle (always overwritten) for effuciency
    for i in 1:NumPar
        CarryTime = 0.0
        if ParStat[i] == 0 # particle is active (not out of bounds)
            RemainTime = Dt # Begin the drift for one timestep
            loop = 1 # loop counter
            Valley = ParValley[i] # identify the current valley the electron occuppies (can only change by scatteirng, not drift, but affects drift)
            PhysValley = ParPhysValley[i] # identify the current 
            CrossFlag = 0
            while RemainTime > 0.0 # iterate until remaining drift time reaches 0.0. Because, when a particle scatters, is drift is sliced into two.
                # extract particle attributes
                initWvk[1] = ParWvk[i,1] # initial particle attributes coming into the loop.
                initWvk[2] = ParWvk[i,2]
                initWvk[3] = ParWvk[i,3]
                #######################
                initPos[1] = ParPos[i,1]
                initPos[2] = ParPos[i,2]
                initPos[3] = ParPos[i,3]
                initBin = BinID(initPos[1])

                if loop == 1 # first loop
                    NewTau = ParTau[i] # during the first loop, we have to remember to conclude to remaining drift until scatter from the previous loop, which is stored in ParTau
                else
                    NewTau = -TauMean*log(randopen()) # if scattering has occured, loop>1, and a new drift time is assigned
                end

                if NewTau <= RemainTime # the drift concludes within the remaining drift time by scattering
                    DriftTime = NewTau
                    RemainTime = RemainTime - NewTau # time remaining since the scattering
                    Wvk, Pos, CrossFlag = Drift(DriftTime, initWvk, initPos, Dist.NodeField, Wvk, Pos, LocNode, initBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, Valley, PhysValley, BarHeight, Mat, Elip, step)
                    if (CrossFlag == 1) || (CrossFlag == 2) # left the domain!
                        ParStat[i] = 1 # break the particle loop
                        break
                    end
                    ScatterBin = BinID(Pos[1])
                    Wvk, GlobalMech, Valley, PhysValley = Scatter(Wvk, ScattTable, ScattCount, ScattIndex, ScatterBin, Mat, Valley, PhysValley, Elip)
                    if GlobalMech != 0
                        Dist.ScattCountBin[GlobalMech, ScatterBin] = Dist.ScattCountBin[GlobalMech, ScatterBin] + 1
                    end
                else
                    DriftTime = RemainTime
                    CarryTime = NewTau - RemainTime
                    RemainTime = 0.0
                    Wvk, Pos, CrossFlag = Drift(DriftTime, initWvk, initPos, Dist.NodeField, Wvk, Pos, LocNode, initBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, Valley, PhysValley, BarHeight, Mat, Elip, step)
                    if (CrossFlag == 1) || (CrossFlag == 2)
                        ParStat[i] = 1
                        break
                    end

                end
                ParWvk[i,1] = Wvk[1]
                ParWvk[i,2] = Wvk[2]
                ParWvk[i,3] = Wvk[3]
                #####################
                ParPos[i,1] = Pos[1]
                ParPos[i,2] = Pos[2]
                ParPos[i,3] = Pos[3]
                ParValley[i] = Valley
                ParPhysValley[i] = PhysValley
                
                loop = loop + 1
            end
            ParTau[i] = CarryTime
        end 
    end

    #println("crossed: ",StepCross)
    return ParWvk, ParPos, ParValley, ParPhysValley, ParTau, ParStat, Dist
end

