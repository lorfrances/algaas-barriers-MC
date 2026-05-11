function Scatter(Wvk, ScattTable, ScattCount, ScattIndex, ScatterBin, Mat, Valley, PhysValley, Elip)
    # this function identifies the relevant scattering mechanism and calls functions to update their energy and momentum
    # compute electron energy (uses elliptic generalization)
    @views Tr = Elip.TrF[ScatterBin, PhysValley, :]                                         # transformation matrix for energy calc.
    # components in transformed isotropic space
    w_x = Wvk[1]*Tr[1] + Wvk[2]*Tr[4] + Wvk[3]*Tr[7]
    w_y = Wvk[1]*Tr[2] + Wvk[2]*Tr[5] + Wvk[3]*Tr[8]
    w_z = Wvk[1]*Tr[3] + Wvk[2]*Tr[6] + Wvk[3]*Tr[9]    
    alpha = Mat.NonParabolic[Valley, ScatterBin]   
    if alpha > 1e-6
        ParEn = 1/(2*alpha)*(-1 + sqrt(1 + 4*alpha*hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0*e_c)))
    else
        ParEn = hbar^2*(w_x^2+w_y^2+w_z^2)/(2*m0*e_c)
    end
    EnInx = floor(Int, ParEn/EnergyStep) 
    frac = (ParEn - EnInx*EnergyStep)/EnergyStep
    if EnInx == 0
       EnInx = 1
       frac = 0
    elseif EnInx >= NumEnergyLevel
       EnInx = NumEnergyLevel-1
       frac = 1
    end
    GammaInx = randopen()*(1/TauMean)
    SelectMech = 0 # self scattering
    GlobalMech = 0 # self scattering
    RateLo = 0.0
    RateHi = 0.0
    # pick appropriate table and number of mechanism
    if Valley == 1
        nMech = ScattCount.G
        Table = ScattTable.G
    elseif Valley == 2
        nMech = ScattCount.L
        Table = ScattTable.L
    else
        nMech = ScattCount.X
        Table = ScattTable.X
    end
    ################################################
    for j in 1:nMech
        RateLo = RateHi
        RateHi = (1-frac)*Table[ScatterBin, EnInx, j] + frac*Table[ScatterBin, EnInx+1,j] + RateHi # interpolated scattering rate
        if (GammaInx >= RateLo) && (GammaInx < RateHi)
            SelectMech = j
            break
        end
    end
    # convert the SelectMech local to valley back to global index to obtain properties of the scattering (assuming not self-scattering)
    if SelectMech != 0
        if Valley == 1
            GlobalMech = SelectMech
        elseif Valley == 2
            GlobalMech = SelectMech + ScattCount.G
        else
            GlobalMech = SelectMech + ScattCount.G + ScattCount.L
        end
        Type = ScattIndex[GlobalMech,2] # type of interaction (1=POP, 2=acoustic, 3=IV...)
        Abs = ScattIndex[GlobalMech,5]
        if Type == 1 #POP scattering
            Wvk = PolarOpticalUpdate(ParEn, w_x, w_y, w_z, ScatterBin, Valley, PhysValley, Wvk, Mat, Abs, Elip)
        elseif Type == 2 || Type == 3 # acoustic phonon scattering or acoustic phonon scattering
            Wvk = AcouPhononUpdate(Wvk, w_x, w_y, w_z, Elip, ScatterBin, PhysValley)
        elseif Type == 4 # ionized impurity scattering
            Wvk = IonizedImpurityUpdate(Wvk, w_x, w_y, w_z, ParEn, Valley, PhysValley, ScatterBin, Mat, Elip)
        elseif Type == 5 #IV phonon scalttering
            Valley_f = ScattIndex[GlobalMech,4]
            CouplingIdx = ScattIndex[GlobalMech,7]
            Wvk, Valley, PhysValley = IVPhononUpdate(Wvk, ParEn, Mat, Valley, PhysValley, Valley_f, ScatterBin, Abs, CouplingIdx, Elip)
        else
            println("ERROR: NO SCATTERING MECHANISM SELECTED (HOW?)")
        end
    end

    return Wvk, GlobalMech, Valley, PhysValley
end

function PolarOpticalUpdate(ParEn, w_x, w_y, w_z, ScatterBin, Valley, PhysValley, Wvk, Mat, Abs, Elip)
    # update particle energy
    alpha = Mat.NonParabolic[Valley, ScatterBin]/e_c
    if Abs == 1
        NewEn = (ParEn + Mat.EpO[ScatterBin])*e_c
    else
        NewEn = (ParEn - Mat.EpO[ScatterBin])*e_c
        if NewEn <= 0.0
            println("Edge case - energy flooring occured: ", ParEn*1e3)
            NewEn = E_min*e_c
        end
    end
    w_p0 = sqrt(w_x^2 + w_y^2 + w_z^2) # original wavevector norm
    w_p = sqrt(((2 * m0 * NewEn*(1 + alpha*NewEn)))/hbar^2) # new electron wavevector based on updated energy
    ge = ParEn*e_c*(1 + alpha*ParEn*e_c)
    gnew = NewEn*(1 + alpha*NewEn)

    # generate polar angle for electron scattering. Here, we work in a rotated coordinate system with the incoming electron direction along the z axis
    # we use Lundstrom's notation (Fundamentals of Carrier Transport) - alpha is polar angle in rotated system, and Beta is the azimuth in the roated system
    zeta = 2*sqrt(ge*gnew)/(ge+gnew-2*sqrt(ge*gnew))
    RandomNum = randopen()
    cos_alpha = ((zeta+1)-(2*zeta+1)^RandomNum)/zeta   # polar angle - most easily defined with cos(alpha) - favors small angles
    sin_alpha = sqrt(1-cos_alpha^2) # from Pythagoras
    beta = 2*pi*randopen()      # random azimuth (full circle)
    cos_beta = cos(beta)     
    sin_beta = sin(beta)          
    
    # now, we find the cartesian coordinates of the new wavevector WITHIN the rotated frame
    wxp = w_p*sin_alpha*cos_beta
    wyp = w_p*sin_alpha*sin_beta
    wzp = w_p*cos_alpha

    # next, we need to find the angles involved in the original rotation, when we moved from the device system to the rotated system with incoming electron direction along the z axis
    # theta and phi, again like Lundstrom
    w_xy = sqrt(w_x^2 + w_y^2)
    cos_theta = w_z/w_p0
    sin_theta = w_xy/w_p0

    if w_xy > 1e-6 # in case w_xy is zero (w ON z)
        cos_phi = w_x/w_xy
        sin_phi = w_y/w_xy
    else
        cos_phi = 1.0
        sin_phi = 0.0
    end

    # now, we have, the angles involved in rotating the coordinate system to the one aligned with invoming electron direction along the z axis. We also know the components of the new electron
    # wavevector within this rotated system. Thus, we just need to rotate back. This is acheived by rotating the aligned coordinate system about the y-axis by theta, and about the z-axis
    # by  phi. These two rotation matrices (about y, about z)are well known, and for both, we take their product. The final expression is
    # w_xyz = Z_phi x Y_theta x w_xyz_p, where Z_phi and Y_theta re the matrices.

    w_x = wxp*cos_phi*cos_theta - wyp*sin_phi + wzp*cos_phi*sin_theta
    w_y = wxp*sin_phi*cos_theta + wyp*cos_phi + wzp*sin_phi*sin_theta
    w_z= -wxp*sin_theta + wzp*cos_theta
    
    # leave the Herring-Vogt space
    @views TrInv = Elip.InvTrF[ScatterBin, PhysValley, :]
    Wvk[1] = w_x * TrInv[1] + w_y * TrInv[4] + w_z * TrInv[7]
    Wvk[2] = w_x * TrInv[2] + w_y * TrInv[5] + w_z * TrInv[8]
    Wvk[3] = w_x * TrInv[3] + w_y * TrInv[6] + w_z * TrInv[9]

    return Wvk
end

function AcouPhononUpdate(Wvk, w_x, w_y, w_z, Elip, ScatterBin, PhysValley)
    # completely randomized process
    w_p = sqrt(w_x^2 + w_y^2 + w_z^2)
    # now, we pick random polar and azimuth angles (random point on sphere), and do not need to rotate to initial wavevector direction frame
    beta = 2*pi*randopen()
    cos_alpha = 1-2*randopen()
    sin_alpha = sqrt(1 - cos_alpha^2)
    w_x = w_p*sin_alpha*cos(beta)
    w_y = w_p*sin_alpha*sin(beta)
    w_z = w_p*cos_alpha

    # leave the Herring-Vogt space
    @views TrInv = Elip.InvTrF[ScatterBin, PhysValley, :]
    Wvk[1] = w_x * TrInv[1] + w_y * TrInv[4] + w_z * TrInv[7]
    Wvk[2] = w_x * TrInv[2] + w_y * TrInv[5] + w_z * TrInv[8]
    Wvk[3] = w_x * TrInv[3] + w_y * TrInv[6] + w_z * TrInv[9]

    return Wvk
end

function IVPhononUpdate(Wvk, ParEn, Mat, Valley, PhysValley, Valley_f, ScatterBin, Abs, CouplingIdx, Elip)
    NewValley = Valley_f
    alpha = Mat.NonParabolic[NewValley, ScatterBin]/e_c
    Valley_sep = RoundToEnergyLevel(Mat.Eg[Valley_f,ScatterBin] - Mat.Eg[Valley,ScatterBin])
    if Abs == 1
        NewEn = (ParEn + Mat.EpIV[CouplingIdx, ScatterBin] -  Valley_sep)*e_c
        if NewEn <= 0.0
            NewEn = E_min*e_c
            println("Edge case - energy flooring occured: ", ParEn*1e3)
        end      
    else
        NewEn = (ParEn - Mat.EpIV[CouplingIdx, ScatterBin] -  Valley_sep)*e_c
        if NewEn <= 0.0
            NewEn = E_min*e_c
            println("Edge case - energy flooring occured: ", ParEn*1e3)
        end
    end
    
    PhysValley_new = NaN
    if Valley_f == 1
        PhysValley_new = 1
    elseif Valley_f == 2 # pick a random L valley
        PhysValley_new = floor(Int,randopen()*4) + 2
        while PhysValley_new == PhysValley
            PhysValley_new = floor(Int,randopen()*4) + 2
        end
    elseif Valley_f == 3
        PhysValley_new = floor(Int,randopen()*3) + 6 # pick a random X valley
        while PhysValley_new == PhysValley
            PhysValley_new = floor(Int,randopen()*3) + 6
        end
    end

    # using new energy, let's find an isotropic distribution for the energy state (assume isotropic scatteirng: same as acoustic modes)
    w_p = sqrt(((2 * m0 * NewEn*(1 + alpha*NewEn)))/hbar^2)
    beta = 2*pi*randopen()
    cos_alpha = 1-2*randopen()
    sin_alpha = sqrt(1 - cos_alpha^2)
    w_x = w_p*sin_alpha*cos(beta)
    w_y = w_p*sin_alpha*sin(beta)
    w_z = w_p*cos_alpha
    
    # leave the Herring-Vogt space
    @views TrInv = Elip.InvTrF[ScatterBin, PhysValley_new, :]
    Wvk[1] = w_x * TrInv[1] + w_y * TrInv[4] + w_z * TrInv[7]
    Wvk[2] = w_x * TrInv[2] + w_y * TrInv[5] + w_z * TrInv[8]
    Wvk[3] = w_x * TrInv[3] + w_y * TrInv[6] + w_z * TrInv[9]
    Valley = NewValley

    return Wvk, Valley, PhysValley_new
end

function IonizedImpurityUpdate(Wvk, w_x, w_y, w_z, ParEn, Valley, PhysValley, ScatterBin, Mat, Elip)
    alpha = Mat.NonParabolic[Valley, ScatterBin]/e_c

    # new electron wavevector same as old 
    # we deal with polar angle dependence like POP scattering
    # note w_p = w_p0 (elastic)
    w_p = sqrt(w_x^2 + w_y^2 + w_z^2)
    #
    N_i = Mat.DopingProfile[ScatterBin]
    gamma = ParEn*e_c*(1 + alpha*ParEn*e_c)
    Beta = sqrt(N_i*e_c^2/(eps0*Mat.eps_s[ScatterBin]*kB*T_pConst))
    eps_beta = hbar^2*Beta^2/(2*m0*Mat.m_ee[Valley, ScatterBin]) # use DOS effective mass...
    RandomNum = randopen()
    cos_alpha = 1 - 2*(1 - RandomNum)/(1 + 4*RandomNum*(gamma/eps_beta))
    sin_alpha = sqrt(1-cos_alpha^2)
    beta = 2*pi*randopen()
    cos_beta = cos(beta)
    sin_beta = sin(beta)
    #
    wxp = w_p*sin_alpha*cos_beta
    wyp = w_p*sin_alpha*sin_beta
    wzp = w_p*cos_alpha
    
    # rotate back
    # first, compute spherical coordinate properties for incoming electron
    w_xy = sqrt(w_x^2 + w_y^2)
    cos_theta = w_z/w_p
    sin_theta = w_xy/w_p
    if w_xy > 1e-6 # in case w_xy is zero (w ON z)
        cos_phi = w_x/w_xy
        sin_phi = w_y/w_xy
    else
        cos_phi = 1.0
        sin_phi = 0.0
    end

    # back-rotation 
    w_x = wxp*cos_phi*cos_theta - wyp*sin_phi + wzp*cos_phi*sin_theta
    w_y = wxp*sin_phi*cos_theta + wyp*cos_phi + wzp*sin_phi*sin_theta
    w_z= -wxp*sin_theta + wzp*cos_theta
    
    # leave the Herring-Vogt space
    @views TrInv = Elip.InvTrF[ScatterBin, PhysValley, :]
    Wvk[1] = w_x * TrInv[1] + w_y * TrInv[4] + w_z * TrInv[7]
    Wvk[2] = w_x * TrInv[2] + w_y * TrInv[5] + w_z * TrInv[8]
    Wvk[3] = w_x * TrInv[3] + w_y * TrInv[6] + w_z * TrInv[9]
    
    return Wvk
end