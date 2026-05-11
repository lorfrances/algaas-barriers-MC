function FixedField(LocNode)
    NodeField = zeros(NumBin)
    NodePotential = zeros(NumBin)
    for i in 1:NumBin
        if (LocNode[i] >= HBPosL && LocNode[i] <= HBPosR)
            NodeField[i] = -ApplVoltage/(HBPosR-HBPosL)
            NodePotential[i] = (-ApplVoltage/(HBPosR-HBPosL))*(LocNode[i] - HBPosL)
        elseif LocNode[i] > HBPosR
            NodePotential[i] = -ApplVoltage
        end
    end
    Dist.NodeField .= NodeField
    Dist.NodePotential .= NodePotential
    return Dist
end

function ChargeDopDist(WeightedCharge, Mat, BinLen)
    ElDop = zeros(NumBin)
    ElCharge = zeros(NumBin)
    for i in 1:NumBin
        ElDop[i] = Mat.DopingProfile[i]
        ElCharge[i] = WeightedCharge[i]/BinLen[i]
    end
    return ElCharge, ElDop
end

function PoissonSolver(ElCharge, ElDop, Dist, Mat)
    NodeField = zeros(NumBin)
    b = -(e_c*(ElDop .- ElCharge)).*(NodeLen^2)./(eps0*Mat.eps_s)
    b[1] = 0.0 # zero potential on left side
    b[NumBin] = -ApplVoltage # total potential on right side
    MainDiag = 2 * ones(NumBin)
    OffDiag = -1 * ones(NumBin-1)
    A = Matrix(Tridiagonal(OffDiag, MainDiag, OffDiag))
    A[1,1] = 1
    A[1,2] = 0
    A[NumBin,NumBin] = 1
    A[NumBin,NumBin-1] = 0
    NodePotential = A\b   
    # evaluate the derivative to determine the electric field. Use a central differnec for internal points, and forward or backward differences for the edges
    for i in 2:NumBin-1
        NodeField[i] = (NodePotential[i+1] - NodePotential[i-1])/(2*NodeLen)
    end
    NodeField[1] = (NodePotential[2] - NodePotential[1])/NodeLen
    NodeField[NumBin] = (NodePotential[NumBin] - NodePotential[NumBin-1])/NodeLen
    Dist.NodeField .= NodeField
    Dist.NodePotential .= NodePotential
    return Dist
end

function PosWeightedGrid(Posx, LocNode, GridVector, Bin)
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
    WeightedVector = WeightN*GridVector[Bin] + WeightNN*GridVector[NextNode]
    return WeightedVector
end
