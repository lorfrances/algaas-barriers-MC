function InjectOhmic(ParPos, ParWvk, ParStat, ParValley, ParPhysValley, ParTau, Dist, Mat, Elip, Contact)
    # function injects particle to the simulation based on maintaining a constant carrier concentration in first bin
    NumContactL, NumContactR = CountContacts(ParPos, ParStat)
    for j in 1:NumOhmic
        if j == 1
            TargetBin = 1
            Dir = 1
            InjPos = DomainL + 1.0e-18 #randopen()*BufferLen + DomainL
            TargetPop = div(RefBinNum,2)
            CurrentPop = NumContactL
        elseif j == 2
            TargetBin = NumBin
            Dir = -1
            InjPos = LenX - 1.0e-18 #LenX - randopen()*BufferLen
            TargetPop = div(RefBinNum,2)
            CurrentPop = NumContactR
        end
        MissingPar = max(0, TargetPop - CurrentPop)
        if MissingPar > 0
            CountInj = 0
            AvailPar = count(==(1), ParStat)
            if MissingPar > AvailPar
                println("Warning: insufficient particles for Ohmic injection: ", AvailPar-MissingPar)
            end
            for i in 1:NumPar
                # update the particle status 
                # check available particle count
                if ParStat[i] != 0 
                    # randomly choose a valley based on equilibrium distributions
                    RandomNum = randopen()
                    if RandomNum < Contact.probValley[1]
                        Valley = 1
                        PhysValley = 1
                    elseif RandomNum < Contact.probValley[2]
                        Valley = 2
                        PhysValley = floor(Int,randopen()*4) + 2
                    else
                        Valley = 3
                        PhysValley = floor(Int,randopen()*3) + 6
                    end

                    ParPos[i,1] = InjPos
                    # select a device-velocity weighted electron wavevector
                    Wvk = ValleyInjection(Valley, PhysValley, Mat, Elip, Contact, 1, Dir)
                    ParWvk[i,1] = Wvk[1]
                    ParWvk[i,2] = Wvk[2]
                    ParWvk[i,3] = Wvk[3]
                    ParStat[i] = 0
                    ParValley[i] = Valley
                    ParPhysValley[i] = PhysValley
                    ParTau[i] = -TauMean*log(randopen())
                    CountInj = CountInj + 1
                end
                if CountInj >= MissingPar
                    break
                end
            end
        end
    end
    return ParPos, ParWvk, ParStat, ParValley, ParPhysValley, ParTau
end

function CountContacts(ParPos, ParStat)
    NumContactL = 0
    NumContactR = 0
    for i in 1:NumPar
        if ParStat[i] == 0
            Posx = ParPos[i,1]
            Bin = BinID(Posx)
            if Bin == 1
                NumContactL = NumContactL + 1
            elseif Bin == NumBin
                NumContactR = NumContactR + 1
            end
        end
    end
    return NumContactL, NumContactR
end

function ValleyInjection(Valley, PhysValley, Mat, Elip, Contact, Bin, Dir)
    Wvk = zeros(3)
    # using the randomly selected valley, we need to select a wavevector corresponding a velocity biased toward inward device direction. We use a Monte Carlo sampling as this is complex.
    # select a random electron energy based on flux-weighted energy distribution (Erlang - sum of independent exponentials):
    E_r = -kB*T_pConst*log(randopen()*randopen())
    # create w space sphere
    alpha_np = Mat.NonParabolic[Valley, Bin]/e_c
    w_p = sqrt(((2 * m0 * E_r*(1 + alpha_np*E_r)))/hbar^2)
    
    # pick a random azimuth
    beta = 2*pi*randopen() # (0, 2\pi)
    
    # now, electrons are biased into the device direction (this is the a vector in w space). So, the polar angle is relative to it, and the "peak velocity" occurs when theta=0.
    alpha = acos(sqrt(1-randopen())) # biased theta angle
    a = Dir*Contact.aVec[PhysValley,:]
    sin_alpha = sin(alpha)
    cos_alpha = cos(alpha)
    cos_beta = cos(beta)
    sin_beta = sin(beta)

    wxp = w_p*sin_alpha*cos_beta
    wyp = w_p*sin_alpha*sin_beta
    wzp = w_p*cos_alpha
    
    # identical procedure to scattering, we rotate our coordinate system so that w_z is aligned with a^ vector, The new components
    # in that coordinate system were calcuated wxp, wyp, wzp. Now we find the angles in that rotation.

    a_p = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    a_xy = sqrt(a[1]^2 + a[2]^2)
    cos_theta = a[3]/a_p
    sin_theta = a_xy/a_p
    if a_xy > 1e-14
        cos_phi = a[1]/a_xy
        sin_phi = a[2]/a_xy
    else
        cos_phi = 1.0
        sin_phi = 0.0
    end

    w_x = wxp*cos_phi*cos_theta - wyp*sin_phi + wzp*cos_phi*sin_theta
    w_y = wxp*sin_phi*cos_theta + wyp*cos_phi + wzp*sin_phi*sin_theta
    w_z= -wxp*sin_theta + wzp*cos_theta
    
    # return to Herring-Vogt 
    @views TrInv = Elip.InvTrF[Bin, PhysValley, :]
    @views MIn = Elip.InvM[Bin, PhysValley, :]  
    Wvk[1] = w_x * TrInv[1] + w_y * TrInv[4] + w_z * TrInv[7]
    Wvk[2] = w_x * TrInv[2] + w_y * TrInv[5] + w_z * TrInv[8]
    Wvk[3] = w_x * TrInv[3] + w_y * TrInv[6] + w_z * TrInv[9]
    Vel_i = hbar/(m0*(1 + 2*alpha_np*E_r)) * (Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7]) 
    if sign(Vel_i) != sign(Dir)
        println("Injection error!")
    end
    return Wvk
end