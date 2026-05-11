function Drift(DriftTime, initWvk, initPos, NodeField, Wvk, Pos, LocNode, initBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, Valley, PhysValley, BarHeight, Mat, Elip, step)
    # the ith electron drifts for drift time
    RemainDriftTime = DriftTime # now, instead of scattering, the drift can be interuppted by barrier interface transitions
    checkInt = 0
    CrossFlag= 0
    while RemainDriftTime > 0.0 # the drift is split into parts if interface is crossed
        
        ##################FIRST APPROXIMATION#######################################
        AddFieldPoi_i = IdentifyFieldPoi(initBin, initPos, NodeField, LocNode) # Poisson field at initial position      
        EField_i[1] = AddFieldPoi_i
        EField_i[2] = 0.0
        EField_i[3] = 0.0
        ###################  
        @views MIn = Elip.InvM[initBin, PhysValley, :]           # extract inverse effective mass tensor
        @views Tr = Elip.TrF[initBin, PhysValley, :]             # extract transformation matrix

        # compute components in transformed isotropic space, for energy calculation
        w_x = initWvk[1]*Tr[1] + initWvk[2]*Tr[4] + initWvk[3]*Tr[7]
        w_y = initWvk[1]*Tr[2] + initWvk[2]*Tr[5] + initWvk[3]*Tr[8]
        w_z = initWvk[1]*Tr[3] + initWvk[2]*Tr[6] + initWvk[3]*Tr[9]       
        alpha = Mat.NonParabolic[Valley, initBin]/e_c # non-parabolic parameter (constant in space)
        
        if alpha > 1e-6
            Ene_i = 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0)))
        else
            Ene_i = hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0) # energy calculation
        end

        Vel_i[1] = hbar/(m0*(1 + 2*alpha*Ene_i)) * (initWvk[1]*MIn[1] + initWvk[2]*MIn[4] + initWvk[3]*MIn[7])    # velocity based on initial wavevector                          # intiial velocity approximation 
        Vel_i[2] = hbar/(m0*(1 + 2*alpha*Ene_i)) * (initWvk[1]*MIn[2] + initWvk[2]*MIn[5] + initWvk[3]*MIn[8])        
        Vel_i[3] = hbar/(m0*(1 + 2*alpha*Ene_i)) * (initWvk[1]*MIn[3] + initWvk[2]*MIn[6] + initWvk[3]*MIn[9])        
        ##################
        Pos_i[1] = initPos[1] + Vel_i[1]*RemainDriftTime                           # position update based on initial velocity
        Pos_i[2] = initPos[2] + Vel_i[2]*RemainDriftTime        
        Pos_i[3] = initPos[3] + Vel_i[3]*RemainDriftTime                           
       
        ###################
        Wvk_i[1] = initWvk[1] + RemainDriftTime*((-e_c*EField_i[1]/hbar))          # wavevector after drift - first approximation
        Wvk_i[2] = initWvk[2] + RemainDriftTime*((-e_c*EField_i[2]/hbar))
        Wvk_i[3] = initWvk[3] + RemainDriftTime*((-e_c*EField_i[3]/hbar))
        ##################
        Bin_i = BinID(Pos_i[1])   


        ####################SECOND APPROXIMATION######################################
        AddFieldPoi_f = IdentifyFieldPoi(Bin_i, Pos_i, NodeField, LocNode)     # Poisson field at the initial updated position approximation
        ################
        EField_f[1] = AddFieldPoi_f
        EField_f[2] = 0.0
        EField_f[3] = 0.0    # total electric field at the first position approximation after drift 
        ###############

        # components in transformed isotropic space for energy calculation
        w_x_i = Wvk_i[1]*Tr[1] + Wvk_i[2]*Tr[4] + Wvk_i[3]*Tr[7]
        w_y_i = Wvk_i[1]*Tr[2] + Wvk_i[2]*Tr[5] + Wvk_i[3]*Tr[8]
        w_z_i = Wvk_i[1]*Tr[3] + Wvk_i[2]*Tr[6] + Wvk_i[3]*Tr[9] 

        if alpha > 1e-6
            Ene_f = 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x_i^2+w_y_i^2+w_z_i^2)/(2*m0)))
        else
            Ene_f = hbar^2*(w_x_i^2+w_y_i^2+w_z_i^2)/(2*m0)
        end  

        Vel_f[1] = hbar/(m0*(1 + 2*alpha*Ene_f)) * (Wvk_i[1]*MIn[1] + Wvk_i[2]*MIn[4] + Wvk_i[3]*MIn[7])                              # intiial velocity approximation 
        Vel_f[2] = hbar/(m0*(1 + 2*alpha*Ene_f)) * (Wvk_i[1]*MIn[2] + Wvk_i[2]*MIn[5] + Wvk_i[3]*MIn[8])        
        Vel_f[3] = hbar/(m0*(1 + 2*alpha*Ene_f)) * (Wvk_i[1]*MIn[3] + Wvk_i[2]*MIn[6] + Wvk_i[3]*MIn[9])   

        #####################RK2 APPROXIMATION########################################
        Pos[1] = initPos[1] + 0.5*RemainDriftTime*(Vel_i[1] + Vel_f[1])
        Pos[2] = initPos[2] + 0.5*RemainDriftTime*(Vel_i[2] + Vel_f[2])     
        Pos[3] = initPos[3] + 0.5*RemainDriftTime*(Vel_i[3] + Vel_f[3])
        ###############                              # RK2 final position
        Wvk[1] = initWvk[1] + RemainDriftTime*((-e_c*0.5*(EField_i[1] + EField_f[1])/hbar))   # RK2 final wavevector
        Wvk[2] = initWvk[2] + RemainDriftTime*((-e_c*0.5*(EField_i[2] + EField_f[2])/hbar))  
        Wvk[3] = initWvk[3] + RemainDriftTime*((-e_c*0.5*(EField_i[3] + EField_f[3])/hbar))  

        ################
        HBCrossCondition, CrossFlag = CrossID(initPos, Pos) # did the particle cross some interface in the process of drifting. 
        # question: why can we use the RK2 to detect crossing? Crossing depends only on position, which depends only on Vel_i and Vel_f. Vel_f only depends on Wvk_i which is based on initial position only! 
        # so, no post-barrier properties which may corrupt the drift prediction is used, assuming effective mass is not changed (which we did not do)

        if CrossFlag == 0 # no crossing event - do not split drift
            UsedTime = RemainDriftTime
        elseif (CrossFlag == 1) || (CrossFlag == 2)# boundary cross                       
            UsedTime = RemainDriftTime        # we will immediately exit the loop
        elseif (CrossFlag == 3) # then HB was crossed, split the drift into parts to handle the interface
            initWvk, initPos, initBin, TimeToInt = DriftToInterface(HBCrossCondition, initPos, initWvk, Pos, Wvk, Wvk_i, Pos_i, Vel_i, Vel_f, RemainDriftTime, NodeField, LocNode, Mat, initBin, EField_i, EField_f, Valley, PhysValley, BarHeight, Elip, step)
            UsedTime = TimeToInt
            if TimeToInt > RemainDriftTime
                println("ERROR: The time to reach an interface has surpassed the drift time!")
            end
            checkInt = checkInt + 1 
        end
        RemainDriftTime = RemainDriftTime - UsedTime   # should only be nonzero due to HB crossing
    end
    if checkInt > 20
        println("WARNING: Particle crossed more than 20 interface crossed in a single timestep!  ", checkInt)
    end
    return Wvk, Pos, CrossFlag
end


function IdentifyFieldPoi(initBin, initPos, NodeField, LocNode)
    # for a given particle position (and bin), identifies the Poisson field there, based either on nearest grid-point or weighted grid (linear) approximation
    if PoiScheme == 2 
        Posx = initPos[1]
        WeightedField = PosWeightedGrid(Posx, LocNode, NodeField, initBin)
        AddFieldPoi = WeightedField 
    elseif PoiScheme == 1
        AddFieldPoi = NodeField[initBin]
    end         
    return AddFieldPoi
end

function CrossID(initPos, Pos)
    # returns whether a boundary was crossed, and if so, the nature of it
    CrossFlag = 0

    if Pos[1] > LenX        # particle leaves through the right boundary
        CrossFlag = 1
    elseif Pos[1] < DomainL # particle leaves through the left boundary
        CrossFlag = 2
    end
    
    ############################################ passing HB  part ########################################################################
    HBCrossL_in = 0
    HBCrossL_out = 0
    HBCrossR_out = 0
    HBCrossR_in = 0
    if HB == 1 && (CrossFlag != 1) && (CrossFlag != 2)
        if (initPos[1] < HBPosL) && (Pos[1] > HBPosL) # crossed barrier movin riht
            HBCrossL_in = 1
            CrossFlag = 3
        elseif (initPos[1] > HBPosL) && (Pos[1] < HBPosL) # crossed barrier moving left
            HBCrossL_out = 1
            CrossFlag = 3
        end

        if (initPos[1] < HBPosR) && (Pos[1] > HBPosR) # crossed barrier movin riht
            HBCrossR_out = 1
            CrossFlag = 3
        elseif (initPos[1] > HBPosR) && (Pos[1] < HBPosR) # crossed barrier moving left
            HBCrossR_in = 1
            CrossFlag = 3
        end
    end
    HBCrossCondition = (HBCrossL_in, HBCrossL_out, HBCrossR_out, HBCrossR_in)

    return HBCrossCondition, CrossFlag
end

function DriftToInterface(HBCrossCondition, initPos, initWvk, Pos, Wvk, Wvk_i, Pos_i, Vel_i, Vel_f, RemainDriftTime, NodeField, LocNode, Mat, initBin, EField_i, EField_f, Valley, PhysValley, BarHeight, Elip, step)
    # Find the time required to reach the interface, and then calculate the then recalculate the drift up to that point
    # identify the interface crossed
    if (HBCrossCondition[1] == 1) || (HBCrossCondition[2] == 1)
        IntPos = HBPosL
    elseif (HBCrossCondition[3] == 1) || (HBCrossCondition[4] == 1)
        IntPos = HBPosR
    end 

    # calculate distance to travel until the interface
    DistToInt = IntPos - initPos[1]
    if DistToInt >= 0.0 
        MovingDir = 1 # moving right
    else
        MovingDir = -1 # moving left
    end

    # calculate the time to reach the interface
    initVelo = Vel_i[1]
    finalVelo = Vel_f[1]
    Accel = (finalVelo - initVelo)/RemainDriftTime
    # solve quadratic equation
    if abs(Accel) > 1e2
        TimeSol1 = (-sqrt(2*Accel*DistToInt + initVelo^2) - initVelo)/Accel
        TimeSol2 = (sqrt(2*Accel*DistToInt + initVelo^2) - initVelo)/Accel
        # the solution is the smallest positive root
        if TimeSol1 >= 0 && TimeSol2 >= 0
            TimeToInt = min(TimeSol1, TimeSol2)
        elseif TimeSol1 >= 0 && TimeSol2 < 0
            TimeToInt = TimeSol1
        elseif TimeSol2 >= 0 && TimeSol1 < 0
            TimeToInt = TimeSol2
        else
            println("TIME ERROR!")
        end
    else
        TimeToInt = DistToInt/initVelo
    end

    if TimeToInt < 0 
        println("ERROR: The time to reach an interface is negative!")
    end


    AddFieldPoi_i = IdentifyFieldPoi(initBin, initPos, NodeField, LocNode) # interpolation can be activated: Poisson field at initial position
    ExpectedPos = IntPos - MovingDir*1.0e-18
    Bin_i = BinID(ExpectedPos) # we must be in the last bin
   
    EField_i[1] = AddFieldPoi_i
    EField_i[2] = 0.0
    EField_i[3] = 0.0  
    
    #################
    @views MIn = Elip.InvM[initBin, PhysValley, :]                                       # extract local effective mass (note, we always use the initial effective mass)
    @views Tr = Elip.TrF[initBin, PhysValley, :]                                         # transformation matrix for energy calc.
    
    # components in transformed isotropic space
    w_x = initWvk[1]*Tr[1] + initWvk[2]*Tr[4] + initWvk[3]*Tr[7]
    w_y = initWvk[1]*Tr[2] + initWvk[2]*Tr[5] + initWvk[3]*Tr[8]
    w_z = initWvk[1]*Tr[3] + initWvk[2]*Tr[6] + initWvk[3]*Tr[9]       
    alpha = Mat.NonParabolic[Valley, initBin]/e_c # non-parabolic parameter (constant in space)
    if alpha > 1e-6
        Ene= 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0)))
    else
        Ene= hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0)
    end
    Vel_i[1] = hbar/(m0*(1 + 2*alpha*Ene)) * (initWvk[1]*MIn[1] + initWvk[2]*MIn[4] + initWvk[3]*MIn[7])                              # intiial velocity approximation 
    Vel_i[2] = hbar/(m0*(1 + 2*alpha*Ene)) * (initWvk[1]*MIn[2] + initWvk[2]*MIn[5] + initWvk[3]*MIn[8])        
    Vel_i[3] = hbar/(m0*(1 + 2*alpha*Ene)) * (initWvk[1]*MIn[3] + initWvk[2]*MIn[6] + initWvk[3]*MIn[9])

    #################   
    Pos_i[1] = initPos[1] + Vel_i[1]*TimeToInt
    Pos_i[2] = initPos[2] + Vel_i[2]*TimeToInt   
    Pos_i[3] = initPos[3] + Vel_i[3]*TimeToInt 
    
    #################  
    Wvk_i[1] = initWvk[1] + TimeToInt*((-e_c*EField_i[1]/hbar))          # wavevector after drift - first approximation
    Wvk_i[2] = initWvk[2] + TimeToInt*((-e_c*EField_i[2]/hbar))
    Wvk_i[3] = initWvk[3] + TimeToInt*((-e_c*EField_i[3]/hbar))
    
    ####################SECOND APPROXIMATION######################################
    AddFieldPoi_f = IdentifyFieldPoi(Bin_i, ExpectedPos, NodeField, LocNode)                                                # we never smear the field at the barrier interface 
    EField_f[1] = AddFieldPoi_f
    EField_f[2] = 0.0
    EField_f[3] = 0.0    
    #####################                                                  # total electric field at the first position approximation after drift  # velocity at final position based on first approximation of wavevector
   
    #####################RK2 APPROXIMATION#########################################
    Pos[1] = ExpectedPos
    Pos[2] = Pos_i[2]  
    Pos[3] = Pos_i[3]                      # RK2 final position

    Wvk[1] = initWvk[1] + TimeToInt*((-e_c*0.5*(EField_i[1] + EField_f[1])/hbar))   # RK2 final wavevector
    Wvk[2] = initWvk[2] + TimeToInt*((-e_c*0.5*(EField_i[2] + EField_f[2])/hbar))  
    Wvk[3] = initWvk[3] + TimeToInt*((-e_c*0.5*(EField_i[3] + EField_f[3])/hbar))    # RK2 final wavevector
    
    PropVel = Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7]  # velocity SIGN check
    ############# check that wvk did not reverse traveling toward interface#########
    if sign(PropVel) == sign(MovingDir)
        # now, we consider exactly what the particle does at the interface and update the properties for continued drift
        Wvk, Pos= UpdateInt(HBCrossCondition, Wvk, Pos, Mat, Valley, PhysValley, BarHeight, Bin_i, Elip, step)
    else
        Pos[1] = Pos_i[1] # use the reversal point
        println("Particle reversed Wvk on the way to an interface, drift adjusted") # never have I observed this error
    end

    initWvk[1] = Wvk[1]
    initWvk[2] = Wvk[2] # initial values for the continued drift 
    initWvk[3] = Wvk[3]
    initPos[1] = Pos[1]
    initPos[2] = Pos[2]
    initPos[3] = Pos[3]
    initBin = Bin_i
    return initWvk, initPos, initBin, TimeToInt
end

function UpdateInt(HBCrossCondition,  Wvk, Pos, Mat, Valley, PhysValley, BarHeight, Bin, Elip, step)
    #using the "just before interface" properties, update the post-interface properties
    # compute particle energy in eV  
    
    alpha = Mat.NonParabolic[Valley, Bin] # initial position non-parabolic parameter (constant in space) eV-1
    Tr = @views Elip.TrF[Bin, PhysValley, :]  # transformation matrix for energy calc.
    MIn = @views Elip.InvM[Bin, PhysValley, :]

    # components in transformed isotropic space # for energy calculation
    w_x = Wvk[1]*Tr[1] + Wvk[2]*Tr[4] + Wvk[3]*Tr[7]
    w_y = Wvk[1]*Tr[2] + Wvk[2]*Tr[5] + Wvk[3]*Tr[8]
    w_z = Wvk[1]*Tr[3] + Wvk[2]*Tr[6] + Wvk[3]*Tr[9]

    if alpha > 1e-6
        Ene = 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0*e_c)))
    else
        Ene = hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0*e_c)
    end

    if HBCrossCondition[1] == 1 # particle is going to hit the barrier but may transmit - conservation of lateral momentum
        alpha_p = Mat.NonParabolic[Valley, BinID(HBPosL+1.0e-18)] #  position non-parabolic parameter (constant in space) eV-1
        NewE = Ene - BarHeight[Valley]
        # with NewE, we can find new k_x using full expansion of the kinetic energy and conseved k_y, k_z
        # it is a qudratic equation
        @views MIn_p = Elip.InvM[BinID(HBPosL+1.0e-18), PhysValley, :] # inverse mass tensor after HB
        A = MIn_p[1]
        B = 2*Wvk[2]*MIn_p[4] + 2*Wvk[3]*MIn_p[7]
        C = Wvk[2]^2*MIn_p[5] + Wvk[3]^2*MIn_p[9] + 2*Wvk[3]*Wvk[2]*MIn_p[8] - 2*m0*e_c*NewE*(1 + alpha_p*NewE)/hbar^2
        Disc = B^2 - 4*A*C
        FinWvk_x = NaN # arbitrary rejection
        v_x_approach = hbar/(m0*(1 + 2*alpha*Ene)) * (Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7])  
        if Disc >= 0.0
            kxSol1 = (-B + sqrt(Disc))/(2*A)
            kxSol2 = (-B - sqrt(Disc))/(2*A)
            vxSol1 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol1*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            vxSol2 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol2*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            if (vxSol1 > 0.0) && (vxSol2 > 0.0) # require smallest for solution for kxSol
                println("Two roots found with valid velocity, smallest one selected!") # rare event (never have I observed this)
                if vxSol1 >= vxSol2
                    FinWvk_x = kxSol2
                else
                    FinWvk_x = kxSol1
                end
            elseif vxSol1 > 0.0
                FinWvk_x = kxSol1
            elseif vxSol2 > 0.0
                FinWvk_x = kxSol2
            end
        end

        if (NewE >= 0) && (Disc >= 0.0) && (!isnan(FinWvk_x))
            #particle is  transmitted over barrier
            Wvk[1] = FinWvk_x
            Pos[1] = HBPosL + 1.0e-18
        else
            # particle is reflected by barrier
            v_x_reflected = -v_x_approach # velocity is reversed in x direction
            ## solve for new kx given v_x reversal
            Wvk[1] = (v_x_reflected/(hbar/(m0*(1 + 2*alpha*Ene)))  - Wvk[2]*MIn[4] - Wvk[3]*MIn[7])/MIn[1]
            Pos[1] = HBPosL-1.0e-18 # should not feel the field again
        end
        
    elseif HBCrossCondition[2] == 1 # particle falls down the barrier in the other direction
        NewE = Ene + BarHeight[Valley]
        @views MIn_p = Elip.InvM[BinID(HBPosL-1.0e-18), PhysValley, :] # inverse mass tensor after barrier transition
        alpha_p = Mat.NonParabolic[Valley, BinID(HBPosL-1.0e-18)] #  position non-parabolic parameter (constant in space) eV-1
        A = MIn_p[1]
        B = 2*Wvk[2]*MIn_p[4] + 2*Wvk[3]*MIn_p[7]
        C = Wvk[2]^2*MIn_p[5] + Wvk[3]^2*MIn_p[9] + 2*Wvk[3]*Wvk[2]*MIn_p[8] - 2*m0*e_c*NewE*(1 + alpha_p*NewE)/hbar^2
        Disc = B^2 - 4*A*C
        FinWvk_x = NaN # arbitrary rejection
        v_x_approach = hbar/(m0*(1 + 2*alpha*Ene)) * (Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7])  
        if Disc >= 0.0
            kxSol1 = (-B + sqrt(Disc))/(2*A)
            kxSol2 = (-B - sqrt(Disc))/(2*A)
            vxSol1 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol1*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            vxSol2 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol2*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])  
            if (vxSol1 < 0.0) && (vxSol2 < 0.0) # require smallest for solution for kxSol
                println("Two roots found with valid velocity, smallest one selected!")
                if vxSol1 >= vxSol2
                    FinWvk_x = kxSol1
                else
                    FinWvk_x = kxSol2
                end
            elseif vxSol1 < 0.0
                FinWvk_x = kxSol1
            elseif vxSol2 < 0.0
                FinWvk_x = kxSol2
            end
        end
        if (NewE >= 0) && (Disc >= 0.0) && (!isnan(FinWvk_x))
            #particle is  transmitted over barrier
            Wvk[1] = FinWvk_x
            Pos[1] = HBPosL - 1.0e-18
        else
            v_x_reflected = -v_x_approach # velocity is reversed in x direction
            ## solve for new kx given v_x reversal
            Wvk[1] = (v_x_reflected/(hbar/(m0*(1 + 2*alpha*Ene)))  - Wvk[2]*MIn[4] - Wvk[3]*MIn[7])/MIn[1]
            Pos[1] = HBPosL + 1.0e-18 # should not feel the field again
        end

    elseif HBCrossCondition[3] == 1 
        NewE = Ene + BarHeight[Valley]
        @views MIn_p = Elip.InvM[BinID(HBPosR+1.0e-18), PhysValley, :]
        alpha_p = Mat.NonParabolic[Valley, BinID(HBPosR+1.0e-18)] # initial position non-parabolic parameter (constant in space) eV-1
        A = MIn_p[1]
        B = 2*Wvk[2]*MIn_p[4] + 2*Wvk[3]*MIn_p[7]
        C = Wvk[2]^2*MIn_p[5] + Wvk[3]^2*MIn_p[9] + 2*Wvk[3]*Wvk[2]*MIn_p[8] - 2*m0*e_c*NewE*(1 + alpha_p*NewE)/hbar^2
        Disc = B^2 - 4*A*C
        FinWvk_x = NaN # arbitrary rejection
        v_x_approach = hbar/(m0*(1 + 2*alpha*Ene)) * (Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7])  
        if Disc >= 0.0
            kxSol1 = (-B + sqrt(Disc))/(2*A)
            kxSol2 = (-B - sqrt(Disc))/(2*A)
            vxSol1 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol1*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            vxSol2 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol2*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            
            if (vxSol1 > 0.0) && (vxSol2 > 0.0) # require smallest for solution for kxSol
                println("Two roots found with valid velocity, smallest one selected!")
                if vxSol1 >= vxSol2
                    FinWvk_x = kxSol2
                else
                    FinWvk_x = kxSol1
                end          
            elseif vxSol1 > 0.0
                FinWvk_x = kxSol1
            elseif vxSol2 > 0.0
                FinWvk_x = kxSol2
            end
        end
        if (NewE >= 0) && (Disc >= 0.0) && (!isnan(FinWvk_x))
            #particle is  transmitted over barrier
            Wvk[1] = FinWvk_x
            Pos[1] = HBPosR + 1.0e-18

        else
            v_x_reflected = -v_x_approach # velocity is reversed in x direction
            ## solve for new kx given v_x reversal
            Wvk[1] = (v_x_reflected/(hbar/(m0*(1 + 2*alpha*Ene)))  - Wvk[2]*MIn[4] - Wvk[3]*MIn[7])/MIn[1]
            Pos[1] = HBPosR - 1.0e-18 # should not feel the field again
        end
    elseif HBCrossCondition[4] == 1  # we need to ensure that the position is not modified again
        NewE = Ene - BarHeight[Valley]
        @views MIn_p = Elip.InvM[BinID(HBPosR-1.0e-18), PhysValley, :]
        alpha_p = Mat.NonParabolic[Valley, BinID(HBPosR-1.0e-18)] # initial position non-parabolic parameter (constant in space) eV-1
        A = MIn_p[1]
        B = 2*Wvk[2]*MIn_p[4] + 2*Wvk[3]*MIn_p[7]
        C = Wvk[2]^2*MIn_p[5] + Wvk[3]^2*MIn_p[9] + 2*Wvk[3]*Wvk[2]*MIn_p[8] - 2*m0*e_c*NewE*(1 + alpha_p*NewE)/hbar^2
        Disc = B^2 - 4*A*C
        FinWvk_x = NaN # arbitrary rejection
        v_x_approach = hbar/(m0*(1 + 2*alpha*Ene)) * (Wvk[1]*MIn[1] + Wvk[2]*MIn[4] + Wvk[3]*MIn[7])  
        if Disc >= 0.0
            kxSol1 = (-B + sqrt(Disc))/(2*A)
            kxSol2 = (-B - sqrt(Disc))/(2*A)
            vxSol1 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol1*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])
            vxSol2 = hbar/(m0*(1 + 2*alpha_p*NewE)) * (kxSol2*MIn_p[1] + Wvk[2]*MIn_p[4] + Wvk[3]*MIn_p[7])  
            if (vxSol1 < 0.0) && (vxSol2 < 0.0) # require smallest for solution for kxSol
                println("Two roots found with valid velocity, smallest one selected!")
                if vxSol1 >= vxSol2
                    FinWvk_x = kxSol1
                else
                    FinWvk_x = kxSol2
                end  
            elseif vxSol1 < 0.0
                FinWvk_x = kxSol1
            elseif vxSol2 < 0.0
                FinWvk_x = kxSol2
            end
        end
        if (NewE >= 0) && (Disc >= 0.0) && (!isnan(FinWvk_x))
            #particle is  transmitted over barrier
            Wvk[1] = FinWvk_x
            Pos[1] = HBPosR - 1.0e-18
        else
            v_x_reflected = -v_x_approach # velocity is reversed in x direction
            ## solve for new kx given v_x reversal
            Wvk[1] =  (v_x_reflected/(hbar/(m0*(1 + 2*alpha*Ene)))  - Wvk[2]*MIn[4] - Wvk[3]*MIn[7])/MIn[1]
            Pos[1] = HBPosR+1.0e-18 # should not feel the field again
        end
    end

    return Wvk, Pos
end

function DriftKappa(DriftTime, initWvk,  Wvk)
    # the ith electron drifts for drift time
    AddField = ConstField
    EField_i = [AddField, 0.0, 0.0] # total electric field at intiial position
    Wvk .= initWvk .+ DriftTime.*((-e_c.*EField_i./hbar))          # wavevector after drift - first approximation
    return Wvk
end

