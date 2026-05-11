function ScatteringTable(Valley, ScattIndexValley, NumScatMechValley, Mat) #valley-resolved scattering rate
    ScattTable = zeros(NumBin, NumEnergyLevel, NumScatMechValley)
    for i in 1:NumScatMechValley
        Type = ScattIndexValley[i,2] # type of interaction (1=POP, 2=acoustic, 3=IV...)
        Abs = ScattIndexValley[i,5] # phonon absorption (1/0 - yes/no)
        for j in 1:NumBin
            T = T_pConst
            for k in 1:NumEnergyLevel
                E_k = k*EnergyStep*e_c # kinetic energy
                if Type == 1
                    ScattTable[j, k, i] = PolarRate(E_k, Mat, Valley, j, T, Abs)
                elseif Type == 2
                    ScattTable[j, k, i] = AcousticRate(E_k, Mat, Valley, j, T)
                elseif Type == 3
                    ScattTable[j, k, i] = AlloyRate(E_k, Mat, Valley, j, T)
                elseif Type == 4
                    ScattTable[j, k, i] = ImpurityRate(E_k, Mat, Valley, j, T)
                elseif Type == 5
                    Valley_f = ScattIndexValley[i,4]
                    CouplingIdx = ScattIndexValley[i,7]
                    Z_f = ScattIndexValley[i,8]
                    if Valley_f - Valley < 0 # downconversion
                        dir = -1
                    else
                        dir = 1
                    end
                    ScattTable[j, k, i] = IVRate(Z_f, E_k, Mat, Valley, Valley_f, j, Abs, CouplingIdx, T)
                end
            end
        end
    end
    return ScattTable
end

function PolarRate(E_k, Mat, Valley, bin, T_p, Abs) # polar optical phonon scattering rate (absorption + emission)
    Omega_pO = Mat.EpO[bin]*e_c/hbar
    alpha = Mat.NonParabolic[Valley,bin]/e_c
    if Abs == 1
        E_f = E_k + Mat.EpO[bin]*e_c
    else
        E_f = E_k - Mat.EpO[bin]*e_c
    end
    if E_f/e_c > 1e-6
        paE = 1 + alpha * E_k
        paE_f = 1 + alpha * E_f
        gamma = E_k*paE
        gamma_f = E_f*paE_f
        A = (2*paE*paE_f + alpha*(gamma + gamma_f))^2
        B = -2*alpha*gamma^0.5*gamma_f^0.5 * (4*paE*paE_f + alpha*(gamma + gamma_f))
        C = 4*paE*paE_f*(1+2*alpha*E_k)*(1+2*alpha*E_f)
        F0 = 1/C *(A * log(abs((gamma^0.5 + gamma_f^0.5)/(gamma^0.5-gamma_f^0.5))) + B)
        Const = e_c^2 * (Mat.m_ee[Valley,bin]*m0)^0.5 * Omega_pO/(4*pi*2^0.5*hbar) * ((1/(eps0*Mat.eps_inf[bin]))-(1/(eps0*Mat.eps_s[bin])))
        f_p = 1/(exp((Mat.EpO[bin]*e_c)/(kB*T_p)) - 1)
        if Abs == 1
            W = Const*((1 + 2*alpha*E_f)/(gamma^0.5))*F0*f_p
            return W
        else
            W = Const*((1 + 2*alpha*E_f)/(gamma^0.5))*F0*(f_p + 1)
            return W
        end
    else
        return 0.0
    end
end

function AcousticRate(E_k, Mat, Valley, bin, T) # aoustic phonon scattering rate
    alpha = Mat.NonParabolic[Valley, bin]/e_c
    gamma = E_k*(1 + alpha*E_k)
    F_a = ((1+alpha*E_k)^2 + 1/3*(alpha*E_k)^2)/(1+2*alpha*E_k)^2
    W_ac = (((2*Mat.m_ee[Valley,bin]*m0)^1.5 * kB*T * (e_c*Mat.DefPotA[Valley, bin])^2)/(2*pi*Mat.MatDens[bin]*Mat.VelSound[bin]^2*hbar^4))*gamma^0.5*(1 + 2*alpha*E_k)*F_a
    return W_ac
end

function AlloyRate(E_k, Mat, Valley, bin, T)
    alpha = Mat.NonParabolic[Valley, bin]/e_c
    gamma = E_k*(1 + alpha*E_k)
    V0 = (Mat.LatConst[bin]*1e-10)^3/4
    DOS = (m0*Mat.m_ee[Valley,bin])^1.5*gamma^0.5*(1 + 2*alpha*E_k)/(2^0.5*pi^2*hbar^3)
    I = (1 + 2*alpha*E_k + (4/3)*alpha^2*E_k^2)/(1 + 2*alpha*E_k)^2
    W_all = 2*pi/hbar * (Mat.AlloyConst[Valley]*e_c)^2 * Mat.AlContent[bin]*(1 - Mat.AlContent[bin]) * V0 * DOS * I  
    return W_all
end

function ImpurityRate(E_k, Mat, Valley, bin, T)
    N_i = Mat.DopingProfile[bin]
    #println(N_i)
    alpha = Mat.NonParabolic[Valley,bin]/e_c
    gamma = E_k*(1 + alpha*E_k)
    beta = sqrt(N_i*e_c^2/(eps0*Mat.eps_s[bin]*kB*T))
    eps_beta = hbar^2*beta^2/(2*m0*Mat.m_ee[Valley, bin])
    preFactor = 2^2.5 * pi * N_i * e_c^4/((4*pi*eps0*Mat.eps_s[bin])^2 * eps_beta^2 * sqrt(m0*Mat.m_ee[Valley,bin]))
    W_II = preFactor * sqrt(gamma) * (1 + 2*alpha*E_k)/(1 + 4*gamma/eps_beta)
    return W_II
end


function IVRate(Z_f, E_k, Mat, Valley, Valley_f, bin, Abs, CouplingIdx, T)
    Omega_piv = Mat.EpIV[CouplingIdx,bin]*e_c/hbar
    f_p = 1/(exp((Mat.EpIV[CouplingIdx,bin]*e_c)/(kB*T)) - 1)
    alpha_i = Mat.NonParabolic[Valley, bin]/e_c
    alpha_j = Mat.NonParabolic[Valley_f, bin]/e_c
    Valley_sep = RoundToEnergyLevel(Mat.Eg[Valley_f,bin] - Mat.Eg[Valley,bin])
    if Abs == 1
        E_f = E_k + e_c*Mat.EpIV[CouplingIdx, bin] - Valley_sep*e_c
    else
        E_f = E_k - e_c*Mat.EpIV[CouplingIdx, bin] - Valley_sep*e_c 
    end
    G_ij = 1.0   #(1 + alpha_i*E_k)*(1 + alpha_j*E_f)/((1 + 2*alpha_i*E_k)*(1 + 2*alpha_j*E_f)) # no longer assume overlap factor
    #println(G_ij)
    gamma_j = E_f*(1 + alpha_j*E_f)
    Const = Z_f * (Mat.m_ee[Valley_f, bin]*m0)^1.5 * (e_c*Mat.DefPotIV[CouplingIdx,bin])^2 / (2^0.5 * pi * Mat.MatDens[bin] * Omega_piv * hbar^3)
    if E_f/e_c > 1e-6
        if Abs == 1
            W_IV = Const*gamma_j^0.5*(1 + 2*alpha_j*E_f)*G_ij*f_p
            return W_IV
        else
            W_IV = Const*gamma_j^0.5*(1 + 2*alpha_j*E_f)*G_ij*(f_p+1)
            return W_IV
        end
    else
        return W_IV = 0.0
    end
end