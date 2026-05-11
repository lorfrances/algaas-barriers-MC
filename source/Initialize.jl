function RoundToEnergyLevel(E)
    # rounds input to nearest energy level for binning
    val = round(round(E*NumEnergyLevel)/NumEnergyLevel; digits=6)
    return val
end

function InitializeMaterial()
    # intializes the material including bin dependent properties and scattering table
    ## bin-dependent initializiation
    
    EpO = zeros(NumBin)
    eps_s = zeros(NumBin)
    eps_inf = zeros(NumBin)
    DefPotA = zeros(NumValley, NumBin)
    VelSound = zeros(NumBin)
    MatDens = zeros(NumBin)
    DefPotIV = zeros(NumValley+NumValley-1, NumBin) # number of intervalley total couplings is the number of valleys plus the number of valleys which are degenerate
    EpIV = zeros(NumValley+NumValley-1, NumBin)
    Eg = zeros(NumValley, NumBin)
    m_ee = zeros(NumValley, NumBin) # DOS effective mass
    ValleySplit = zeros(NumValley+NumValley-1, NumBin)
    BarHeight = zeros(NumValley) # each valley has its own barrier height
    DopingProfile = zeros(NumBin) .+ DopDens_n
    LatConst = zeros(NumBin)
    AlContent = zeros(NumBin)
    ThermCond = zeros(NumBin)
    NonParabolic = zeros(NumValley, NumBin)
    AlloyConst = zeros(NumValley)

    # initialize tensors (flat)
    TrF = zeros(NumBin, NumPhysValley, 9)
    InvTrF = zeros(NumBin, NumPhysValley, 9)
    InvM = zeros(NumBin, NumPhysValley, 9)
    Rot = zeros(NumBin, NumPhysValley, 9)
    aVec = zeros(NumPhysValley,3)
    ########### run parameters
    LocNode = zeros(NumBin) # posiitions of nodes (center of cells)
    BinLen = zeros(NumBin) # length of bins (or cells)
    xMatB = FindxMat(HBHeight)
    println("BARRIER Al:",xMatB)

    ###############
    for i in 1:NumBin
        BinCenter = round((i*NodeLen) - 0.5*NodeLen, digits=12)  
        BinEdge = round((i*NodeLen), digits=12)
        LocNode[i] = BinCenter # position of center of bins (we call them nodes) - vector needed
        if i == 1 || i == NumBin
            BinLen[i] = NodeLen/2
        else
            BinLen[i] = NodeLen
        end

        if (HB == 1) #simulate HB - identify x_Al everywhere
            if (BinCenter >= HBPosL) && (BinCenter <= HBPosR)
                AlFraction = xMatB
            else
                AlFraction = xMatA
            end
        else
            AlFraction = xMatA
        end
        
        EpO[i], eps_s[i], eps_inf[i], DefPotA[1,i],DefPotA[2,i],DefPotA[3,i], VelSound[i], MatDens[i], DefPotIV[1,i], DefPotIV[4,i], DefPotIV[2,i], DefPotIV[3,i], DefPotIV[5,i],
        EpIV[1,i], EpIV[4,i], EpIV[2,i], EpIV[3,i], EpIV[5,i], Eg[1,i], Eg[2,i], Eg[3,i], m_ee[1,i], m_ee[2,i], m_ee[3,i], ValleySplit[1,i],
        ValleySplit[2,i], ValleySplit[3,i], ValleySplit[4,i], ValleySplit[5,i], LatConst[i], AlloyConst[1],AlloyConst[2],AlloyConst[3], AlContent[i], NonParabolic[1,i], NonParabolic[2,i], NonParabolic[3,i], ThermCond[i] = MatParVar(AlFraction) 

        if (HB == 1) && (BinCenter >= HBPosL) && (BinCenter <= HBPosR)
           DopingProfile[i] = DopDens_n*0.05
        end

        ######## assign tensors and matrices 
        for j in 1:NumPhysValley #8 physical valleys
            TrF[i, j, :], InvTrF[i, j, :], InvM[i, j, :], Rot[i, j, :], aVec[j,:] = ReturnElipMat(AlFraction, j)
        end
    end

    # define the material structure as a named tuple:
    
    Mat = (EpO=EpO, eps_s=eps_s, eps_inf=eps_inf,DefPotA=DefPotA, VelSound=VelSound,MatDens=MatDens,DefPotIV=DefPotIV,EpIV=EpIV, Eg = Eg, m_ee=m_ee, ValleySplit=ValleySplit, LatConst=LatConst, AlloyConst=AlloyConst, AlContent =AlContent, NonParabolic=NonParabolic, DopingProfile=DopingProfile, ThermCond=ThermCond)
    Elip = (TrF=TrF, InvTrF=InvTrF, InvM=InvM, Rot=Rot)
    #barrier height (if HB is on)
    if HB == 1
        for i in 1:NumValley
            if i == 1
                BarHeight[i] =  RatEcEv*(Eg[i,BinID(HBPosL + 1.0e-18)] - Eg[i, BinID(HBPosL - 1.0e-18)])
            else
                BarHeight[i] =  (Eg[i,BinID(HBPosL + 1.0e-18)] - Eg[i, BinID(HBPosL - 1.0e-18)]) - (1-RatEcEv)*(Eg[1,BinID(HBPosL + 1.0e-18)] - Eg[1, BinID(HBPosL - 1.0e-18)])
            end
        end
    end

    # equilibrium thermalized multivalley population for Ohmic contact physics
    alpha_G = 0.5*1*(Mat.m_ee[1,1]*m0)^1.5*pi^0.5*T_pConst^1.5*kB^1.5
    alpha_L = 0.5*4*(Mat.m_ee[2,1]*m0)^1.5*pi^0.5*T_pConst^1.5*kB^1.5*exp(-Mat.ValleySplit[1,1]*e_c/(kB*T_pConst))
    alpha_X = 0.5*3*(Mat.m_ee[3,1]*m0)^1.5*pi^0.5*T_pConst^1.5*kB^1.5*exp(-Mat.ValleySplit[2,1]*e_c/(kB*T_pConst))
    E_f = kB*T_pConst*log(pi^2*hbar^3*DopDens_n/(2^0.5*(alpha_G+alpha_L+alpha_X)))
    # dist
    nG = 2^0.5*(Mat.m_ee[1,1]*m0)^1.5/(pi^2*hbar^3) * (0.5*pi^0.5*T_pConst^1.5*kB^1.5*exp(E_f/(kB*T_pConst)))
    nL = 2^0.5*4*(Mat.m_ee[2,1]*m0)^1.5/(pi^2*hbar^3) * (0.5*pi^0.5*T_pConst^1.5*kB^1.5*exp((E_f-Mat.ValleySplit[1,1]*e_c)/(kB*T_pConst)))
    nX = 2^0.5*3*(Mat.m_ee[3,1]*m0)^1.5/(pi^2*hbar^3) * (0.5*pi^0.5*T_pConst^1.5*kB^1.5*exp((E_f-Mat.ValleySplit[2,1]*e_c)/(kB*T_pConst)))
    ######
    pG = nG/(nG+nL+nX)
    pL = nL/(nG+nL+nX)
    pX = nX/(nG+nL+nX)
    #####
    probValley = cumsum([pG, pL, pX])

    Contact = (probValley=probValley, aVec=aVec) # this structure is for properties of contacts



    ######## initialize the scattering using the valley-resolved scattring 
    ScattIndex = readdlm("ScatteringIndex.csv", ',', Any, '\n')
    ScattIndex = Int.(ScattIndex[2:end, 2:end])
    ScattTableG = nothing
    ScattTableL = nothing
    ScattTableX = nothing
    for i in 1:NumValley
        ScattIndexValley = ScattIndex[ScattIndex[:,3] .== i, :] #set of all scattering mechanisms corresponding to valley i.
        NumScatMechValley = size(ScattIndexValley,1)
        if i == 1
            ScattTableG = ScatteringTable(i, ScattIndexValley, NumScatMechValley, Mat)
        elseif i == 2
            ScattTableL = ScatteringTable(i, ScattIndexValley, NumScatMechValley, Mat)
        else
            ScattTableX = ScatteringTable(i, ScattIndexValley, NumScatMechValley, Mat)
        end       
    end
    ScattTable = (G = ScattTableG, L = ScattTableL, X = ScattTableX)
    ScattCount = (G = size(ScattTable.G,3), L = size(ScattTable.L,3), X = size(ScattTable.X,3))
    # take a sum for all E over all bins and valleys
    MaxRate = maximum([maximum(sum(ScattTable.G, dims=3)[:, :, 1]), maximum(sum(ScattTable.L,dims=3)[:,:,1]), maximum(sum(ScattTable.X,dims=3)[:,:,1])])

    return Mat, LocNode, BinLen, BarHeight, ScattTable, MaxRate, ScattCount, ScattIndex, Elip, Contact
end

function InitializeParticles(Mat, LocNode)
    # inializes the critcal particle-wise parameters

    ParStat = zeros(Int,NumPar)
    ParTau = zeros(NumPar)
    ParPos = zeros(NumPar,3)
    ParWvk = zeros(NumPar,3)
    ParValley = zeros(Int,NumPar)
    ParPhysValley = zeros(Int, NumPar)
    TotPar = 0

    # Assign initial particle energy and positions and free flight
    for i in 1:NumBin
        for j in 1:ParPerBin
            inx = (i-1)*ParPerBin + j
            if LocNode[i] > HBPosL && LocNode[i] < HBPosR
                ParStat[inx] = 1 #we will remove all electrons from lightly doped domain, initially
            else
                ParWvk[inx,1] = (Mat.m_ee[1,1]*m0/hbar)*randn()*sqrt(kB*T_pConst/(Mat.m_ee[1,1]*m0))
                ParWvk[inx,2] = (Mat.m_ee[1,1]*m0/hbar)*randn()*sqrt(kB*T_pConst/(Mat.m_ee[1,1]*m0))
                ParWvk[inx,3] = (Mat.m_ee[1,1]*m0/hbar)*randn()*sqrt(kB*T_pConst/(Mat.m_ee[1,1]*m0))
                Posx_i = (i-1)*NodeLen + randopen()*NodeLen
                ParPos[inx,1] = clamp(Posx_i, DomainL, LenX)
                ParTau[inx] = -TauMean*log(randopen())
                ParValley[inx] = 1
                ParPhysValley[inx] = 1
            end
            TotPar = TotPar + 1
        end
    end
    AvailPar = count(==(1), ParStat)
    println("Particles initially in contacts: ", AvailPar)

    return ParPos, ParWvk, ParTau, ParStat, ParValley, ParPhysValley
end

function FindxMat(HBHeight)
    # solves for the alloy fraction required for a given Γ barrier height (use bisection method)
    Ga_EgG = 1.4224821428571428
    Al_EgG = 3.0030361445783136

    EgG_A = Al_EgG*xMatA + Ga_EgG*(1-xMatA) - (-0.127 +1.310*xMatA)*xMatA*(1-xMatA)
    a = xMatA # must be between the alloy fraction in section A and 1.
    b = 1
    cprev = 1
    c = (a+b)*0.5
    while abs(c-cprev) > 0.00001
        BarrierHeight_a = RatEcEv*((Al_EgG*a + Ga_EgG*(1-a) - (-0.127 +1.310*a)*a*(1-a)) - EgG_A)
        BarrierHeight_c = RatEcEv*((Al_EgG*c + Ga_EgG*(1-c) - (-0.127 +1.310*c)*c*(1-c)) - EgG_A)
        f_a =  BarrierHeight_a - HBHeight
        f_c =  BarrierHeight_c - HBHeight
        if f_a*f_c<0
            b = c
        else
            a = c
        end
        cprev = c
        c = (a+b)*0.5
    end
    return c
end


function ReturnElipMat(x, v)
    # returns the matrices essential for the anisotropic treatment
    # 1. trivial gamma valley.
    T_v = zeros(9)
    T_inv_v = zeros(9)
    m_inv_v = zeros(9)
    if v == 1 # Gamma valley - trivial [000]
        m_ee = 0.15*x + 0.067*(1-x)
        R = [1 0 0 # rotation matrix to valley-aligned coordinate
            0 1 0
            0 0 1]
        t = [sqrt(1/m_ee) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee) 0
            0 0 sqrt(1/m_ee)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee 0
            0 0 1/m_ee]
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        m_inv_R = transpose(R)*m_inv*R
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 2 # L valley [111]
        m_ee_l = 1.9*(1-x) + 1.32*x
        m_ee_t = 0.0754*(1-x) + 0.15*x
        R = [1/sqrt(3) 1/sqrt(3) 1/sqrt(3) # rotation matrix to valley-aligned coordinate
            -1/sqrt(2) 1/sqrt(2) 0
            -1/sqrt(6) -1/sqrt(6) sqrt(2)/sqrt(3)]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        
        m_inv_R = transpose(R)*m_inv*R
 
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 3 # L valley [-111]
        m_ee_l = 1.9*(1-x) + 1.32*x
        m_ee_t = 0.0754*(1-x) + 0.15*x
        R = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) # rotation matrix to valley-aligned coordinate
            -1/sqrt(2) -1/sqrt(2) 0
            1/sqrt(6) -1/sqrt(6) sqrt(2)/sqrt(3)]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 4 # L valley [1-11]
        m_ee_l = 1.9*(1-x) + 1.32*x
        m_ee_t = 0.0754*(1-x) + 0.15*x
        R = [1/sqrt(3) -1/sqrt(3) 1/sqrt(3) # rotation matrix to valley-aligned coordinate
            1/sqrt(2) 1/sqrt(2) 0
            -1/sqrt(6) 1/sqrt(6) sqrt(2)/sqrt(3)]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 5 # L valley [11-1]
        m_ee_l = 1.9*(1-x) + 1.32*x
        m_ee_t = 0.0754*(1-x) + 0.15*x
        R = [1/sqrt(3) 1/sqrt(3) -1/sqrt(3) # rotation matrix to valley-aligned coordinate
            -1/sqrt(2) 1/sqrt(2) 0
            1/sqrt(6) 1/sqrt(6) sqrt(2)/sqrt(3)]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 6 # X valley [1 0 0]
        m_ee_l = 1.3*(1-x) + 0.97*x
        m_ee_t = 0.23*(1-x) + 0.22*x
        R = [1 0 0 # rotation matrix to valley-aligned coordinate
            0 1 0
            0 0 1]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 7 #  valley [0 1 0]
        m_ee_l = 1.3*(1-x) + 0.97*x
        m_ee_t = 0.23*(1-x) + 0.22*x
        R = [0 1 0 # rotation matrix to valley-aligned coordinate
            0 0 1
            1 0 0]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    if v == 8 # Gamma valley [001]
        m_ee_l = 1.3*(1-x) + 0.97*x
        m_ee_t = 0.23*(1-x) + 0.22*x
        R = [0 0 1 # rotation matrix to valley-aligned coordinate
            1 0 0
            0 1 0]
        t = [sqrt(1/m_ee_l) 0 0 # transformation matrix IN valley-alligned ccoordinate
            0 sqrt(1/m_ee_t) 0
            0 0 sqrt(1/m_ee_t)]
        T = t * R # rotated transformation matrix (rotate to valley-frame)
        m_inv = [1/m_ee_l 0 0 # effective mass inverse tensor - valley aligned
            0 1/m_ee_t 0
            0 0 1/m_ee_t]
        m_inv_R = transpose(R)*m_inv*R
        T_inv = inv(T) # inverted transformation matrix (to rotate back to cartesian)
        ############## flatten for printout
        T_v = vec(T)
        T_inv_v = vec(T_inv)
        m_inv_v = vec(m_inv_R)
        r = vec(R)
        # contact physics:
        nForward = [1 
                    0
                    0] # device direction +x (z) - unit vector
        aVectorF = transpose(T_inv) * (m_inv_R * nForward)
        aVector = aVectorF / norm(aVectorF)
    end
    return T_v, T_inv_v, m_inv_v, r, aVector
end


function MatParVar(x) 
    # material parameters as a function of the Al alloying
    # this first section is assumed to have constant properties (for now, easily changed.)

    DefPotA_G = 7.00 
    DefPotA_L = 9.20
    DefPotA_X = 9.00
    AlloyConstG = 0.5 #eV 
    AlloyConstL = 0.5 #eV
    AlloyConstX = 0.5 #eV
    IVGL = 0.8e11 #eV/m
    IVLL = 1.0e11 # eV/m
    IVGX = 1.0e11
    IVLX = 0.15e11
    IVXX = 1.0e11


    # variable properties
    EpO = RoundToEnergyLevel(0.036 + 0.014*x)
    IVEpGL = RoundToEnergyLevel(0.028 + 0.014*x)
    IVEpLL = RoundToEnergyLevel(0.029 + 0.014*x)
    IVEpGX = RoundToEnergyLevel(0.030 + 0.014*x)
    IVEpLX = RoundToEnergyLevel(0.029 + 0.014*x)
    IVEpXX = RoundToEnergyLevel(0.030 + 0.014*x)
    
    
    LatConst = 5.6611*x + 5.65235*(1-x)
    eps_s = 10.06*xMatA + 12.90*(1-xMatA)
    eps_inf = 8.16*xMatA + 10.89*(1-xMatA)
    MatDens = 3760*x + 5360*(1-x)
    C_11  = 0.100*(12.02e11*x + 11.88e11*(1-x))
    v_l_1 = (C_11/MatDens)^0.5
    VelSound = v_l_1
    Ga_EgG = 1.4224821428571428
    Ga_EgX = 1.898857142857143
    Ga_EgL = 1.7069642857142857
    Al_EgG = 3.0030361445783136
    Al_EgX = 2.164096385542169
    Al_EgL = 2.3519642857142857
    bG = -0.127 +1.310*x
    bX = 0.055
    bL = 0.0

    ########################################### note: some variables have valley dependence
    EgG = Al_EgG*x + Ga_EgG*(1-x) - bG*x*(1-x)
    EgX = Al_EgX*x + Ga_EgX*(1-x) - bX*x*(1-x)
    EgL = Al_EgL*x + Ga_EgL*(1-x) - bL*x*(1-x)

    m_ee_l_L = 1.9*(1-x) + 1.32*x
    m_ee_t_L = 0.0754*(1-x) + 0.15*x
    m_ee_l_X = 1.3*(1-x) + 0.97*x
    m_ee_t_X = 0.23*(1-x) + 0.22*x

    # DENSITY OF STATES effective mass, for scattering rates. FULL anisotropic are stored in another structure.
    m_eeG = 0.15*x + 0.067*(1-x)
    m_eeL = (m_ee_l_L*m_ee_t_L^2)^(1/3)
    m_eeX = (m_ee_l_X*m_ee_t_X^2)^(1/3)
    #

    SplitGL = RoundToEnergyLevel(EgL - EgG)   # rounding to keep the scattering code stable.
    SplitGX = RoundToEnergyLevel(EgX - EgG)
    SplitLX = RoundToEnergyLevel(EgX - EgL)
    SplitLL = 0.0
    SplitXX = 0.0

    # non-parabolicity coefficient
    alphaG = 1/EgG * (1-m_eeG)^2
    alphaL = 0.0
    alphaX = 0.0

    # thermal conductivity
    C_k  = 3.3 # thermal conductivity bowing
    k_Ga = 46
    k_Al = 80
    ThermCond = 1/(((1-x)/k_Ga)+(x/k_Al)+(((1-x)*x)/C_k))

    return EpO, eps_s, eps_inf, DefPotA_G,DefPotA_L,DefPotA_X, VelSound, MatDens, IVGL, IVLL, IVGX, IVLX, IVXX, 
    IVEpGL, IVEpLL, IVEpGX, IVEpLX, IVEpXX, EgG, EgL, EgX, m_eeG, m_eeL, m_eeX, SplitGL, SplitGX, SplitLX, SplitLL, SplitXX, LatConst, AlloyConstG, AlloyConstL, AlloyConstX, x,
    alphaG, alphaL, alphaX, ThermCond
end

function InitialiazeIO()
    # returns file handles for input-output, and writes the first lines
    FieldTimeData =  "FieldTimeData" * suff * ".dat"
    PotentialTimeData =  "PotentialTimeData" * suff * ".dat"
    CurrentDataBinAve =  "CurrentDataBinAve" * suff * ".dat"


    BinCountAve =  "BinCountAve" * suff * ".dat"
    VelDistAve =  "VelDistAve" * suff * ".dat"
    EnDistAve =  "EnDistAve" * suff * ".dat"
    ParEneDistData =  "ParEneDistData" * suff * ".dat"
    ParEneDistDataG =  "ParEneDistDataG" * suff * ".dat"
    ParEneDistDataL =  "ParEneDistDataL" * suff * ".dat"
    ParEneDistDataX =  "ParEneDistDataX" * suff * ".dat"
    ParVelDistData =  "ParVelDistData" * suff * ".dat"
    CurrentSpectra =  "CurrentSpectra" * suff * ".dat"
    CurrentCross =  "CurrentCross" * suff * ".dat"
    LValleyData = "LValleyCount" * suff * ".dat"
    GValleyData = "GValleyCount" * suff * ".dat"
    XValleyData =  "XValleyCount" * suff * ".dat"
    
    ## phonon absorption decomposition
    GPhononCountData = "GPhononCountData" * suff * ".dat"
    LPhononCountData = "LPhononCountData" * suff * ".dat"
    XPhononCountData = "XPhononCountData" * suff * ".dat"
    GLPhononCountData = "GLPhononCountData" * suff * ".dat"
    GXPhononCountData = "GXPhononCountData" * suff * ".dat"
    LGPhononCountData = "LGPhononCountData" * suff * ".dat"
    LXPhononCountData = "LXPhononCountData" * suff * ".dat"
    LLPhononCountData = "LLPhononCountData" * suff * ".dat"
    XGPhononCountData = "XGPhononCountData" * suff * ".dat"
    XLPhononCountData = "XLPhononCountData" * suff * ".dat"
    XXPhononCountData = "XXPhononCountData" * suff * ".dat"
    GLTotAbsorbData = "GLTotAbsorb" * suff * ".dat"
    GLTotEmitData = "GLTotEmit" * suff * ".dat"
    LGTotAbsorbData= "LGTotAbsorb" * suff * ".dat"
    LGTotEmitData = "LGTotEmit" * suff * ".dat"  


    # write intial file message
    open(GPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(LPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(XPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(GLPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(GXPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(LGPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(LXPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(LLPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end
    open(XGPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(XLPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(XXPhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(GLTotAbsorbData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end 

    open(GLTotEmitData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end 

    open(LGTotAbsorbData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end 

    open(LGTotEmitData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end 

    open(EnDistAve, "w") do io
        println(io,"AVERAGE PARTICLE ENERGY (eV) in bins (columns)")
    end

    open(FieldTimeData, "w") do io
        println(io,"AVERAGE ELECTRIC FIELD (V/m) in bins (columns)")
    end

    open(PotentialTimeData, "w") do io
        println(io,"AVERAGE INTERNAL POTENTIAL (V) in bins (columns)")
    end

    open(CurrentDataBinAve, "w") do io
        println(io,"AVERAGE CURRENT (A/cm^2) in bins (columns)")
    end

    open(BinCountAve, "w") do io
        println(io,"AVERAGE PARTICLE COUNT (N) in bins (columns)")
    end

    open(VelDistAve, "w") do io
        println(io,"AVERAGE VELOCITY (m/s) in bins (columns)")
    end

    open(CurrentCross,"w") do io
        println(io, "CURRENT BY PARTICLES CROSSING RIGHT BOUNDARY")
    end

    open(LValleyData,"w") do io
        println(io, "PARTICLES IN L VALLEY BY BIN")
    end

    open(XValleyData,"w") do io
        println(io, "PARTICLES in X VALLEY BY BIN")
    end

    open(GValleyData,"w") do io
        println(io, "PARTICLES IN G VALLEY BY BIN")
    end

    OutFile = (FieldTimeData =  FieldTimeData,             # compile file handles into one tuple for easy use
                PotentialTimeData =  PotentialTimeData,
                CurrentDataBinAve =  CurrentDataBinAve,
                BinCountAve = BinCountAve,
                VelDistAve = VelDistAve ,
                EnDistAve =  EnDistAve,
                ParEneDistData =  ParEneDistData,
                ParEneDistDataG =  ParEneDistDataG,
                ParEneDistDataL =  ParEneDistDataL,
                ParEneDistDataX =  ParEneDistDataX,
                ParVelDistData =  ParVelDistData, 
                CurrentSpectra = CurrentSpectra,
                CurrentCross =  CurrentCross,
                LValleyData = LValleyData,
                GValleyData = GValleyData,
                XValleyData =  XValleyData,
                GPhononCountData = GPhononCountData, 
                LPhononCountData = LPhononCountData,
                XPhononCountData = XPhononCountData,
                GLPhononCountData = GLPhononCountData, 
                GXPhononCountData = GXPhononCountData, 
                LGPhononCountData = LGPhononCountData,
                LXPhononCountData = LXPhononCountData,
                LLPhononCountData = LLPhononCountData,
                XGPhononCountData = XGPhononCountData,
                XLPhononCountData = XLPhononCountData,
                XXPhononCountData = XXPhononCountData,
                GLTotAbsorbData = GLTotAbsorbData,
                GLTotEmitData = GLTotEmitData,
                LGTotAbsorbData = LGTotAbsorbData,
                LGTotEmitData = LGTotEmitData)
    return OutFile
end

function InitializeTotals()
    # initializes the matrix/vectors that are used to store running totals for averaging
    NodeField = zeros(NumBin)
    NodePotential = zeros(NumBin)
    BinCount = zeros(NumBin)
    ValleyDist = zeros(NumValley, NumBin)
    CurrentBin = zeros(NumBin)
    VelDist = zeros(NumBin)
    EneDist = zeros(NumBin)
    ParEneDist = zeros(NumValley+1, NumEnergyLevelProp, NumBin)
    ParVelDist = zeros(NumVelLevel*2 + 1, NumBin)
    CurrentSpectra = zeros(NumEnergyLevelProp, NumBin)
    CrossFreq = 0

    Total = (NodeField=NodeField,
            NodePotential=NodePotential,
            BinCount = BinCount,
            ValleyDist = ValleyDist,
            CurrentBin = CurrentBin,
            VelDist = VelDist,
            EneDist = EneDist,
            ParEneDist = ParEneDist,
            ParVelDist = ParVelDist,
            CurrentSpectra = CurrentSpectra,
            CrossFreq = CrossFreq)
    return Total
end

function InitializeDist()
    # initializes the matrix/vectors which store the instantaneous distributions (these arre recycled for efficiency)
    ScattCountBin = zeros(TotScatt, NumBin) # store all phonon processes in this matrix
    NodeField = zeros(NumBin)
    NodePotential = zeros(NumBin)
    WeightedCharge = zeros(NumBin)
    Current = zeros(NumBin)
    ValleyDist = zeros(NumValley, NumBin)
    BinDist = zeros(NumBin)
    VelDist = zeros(NumBin)
    EneDist = zeros(NumBin)
    ParEneDist = zeros(NumValley+1,NumEnergyLevelProp, NumBin)
    ParVelDist = zeros(NumVelLevel*2+1, NumBin)
    CurrentSpectrum = zeros(NumEnergyLevelProp, NumBin)
    Dist = (ScattCountBin = ScattCountBin,
            NodeField = NodeField,
            NodePotential = NodePotential,
            Valley = ValleyDist,
            Bin = BinDist,
            Vel = VelDist,
            Ene = EneDist,
            ParEne = ParEneDist,
            ParVel = ParVelDist,
            Current = Current,
            CurrentSpectrum = CurrentSpectrum,
            )
    return Dist, WeightedCharge
end