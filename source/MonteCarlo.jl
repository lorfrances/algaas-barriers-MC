# read all Julia files needed
include("Initialize.jl")
include("Properties.jl")
include("ScatteringRates.jl")
include("randopen.jl")
include("ParticleLoop.jl")
include("ParticleDrift.jl")
include("ParticleScatter.jl")
include("PoissonSolver.jl")
include("Contacts.jl")

# packages
using LinearAlgebra
using DelimitedFiles
using Statistics
using Printf

# global physical Constants
const kB =  1.380649e-23
const hbar = 1.054571817e-34
const m0 = 9.1093837e-31
const e_c = 1.60217663e-19
const eps0 = 8.8541878188e-12

# overall calculation conditions
# read-in the input file
input_path = ARGS[1]
InputData = readdlm(input_path, comments=true, comment_char='#')


# calculation conditionsW
const Dt = InputData[1]
const TotTime = InputData[2]
const FreqWriteOut = round.(Int, InputData[3])# property write out frequnecy
const WriteOutTime = FreqWriteOut*Dt
const NumPar = round.(Int,InputData[4])
const NumBin = round.(Int,InputData[5])
const NodeLen = InputData[6] # separation between nodes
const NumEnergyLevel = 5000
const EnergyStep = 2.0e-4
const NumEnergyLevelProp = 1000
const EnergyStepProp = 1.0e-3 # for printing data, we need to reduce file size
const NumVelLevel = 80
const MaxVel = 2e6
const PoiScheme = 2
const VelStep = MaxVel/NumVelLevel
const stepStatisticsBegin = 12000
const DomainL = 2.5e-9

## global material conditions
const NumValley = round.(Int, InputData[7])
const NumPhysValley = 8
const DopDens_n = InputData[8]
const LenX = InputData[9]
const HB = round.(Int, InputData[10])
const HBPosL = InputData[11]
const HBPosR = InputData[12]
const xMatA = InputData[13]
const HBHeight = InputData[14]
const RatEcEv = InputData[15]
const T_pConst = 300
const NumOhmic = 2
const TBR = 2.0e-9
const ApplVoltage = InputData[16]
const E_min = 1e-4 # small energy for handling edge cases (eV)
const RefBinNum = div(NumPar,NumBin) #integer
const ParCharge = DopDens_n*NodeLen/RefBinNum # sheet charge density reperesented by one elctron (electrons/m^2)

# read scattering parameters###
if HB == 0
    const suff = "_NB"   ############ suffix for ouput files (to distinguish whtehr barrier is included)
else
    const suff = ""
end
#############################

Mat, LocNode, BinLen, BarHeight, ScattTable, MaxRate, ScattCount, ScattIndex, Elip, Contact = InitializeMaterial()
println("Maximum rate: ", MaxRate)
const Nc = 2*((kB*T_pConst*Mat.m_ee[1,1]*m0*2*pi)/(hbar*2*pi)^2)^1.5
const FermiEnergy = (kB*T_pConst/e_c)*log(DopDens_n/Nc)
writedlm("ScattTable_G.dat",ScattTable.G[1,:,:])
writedlm("ScattTable_L.dat",ScattTable.L[1,:,:])
writedlm("ScattTable_X.dat",ScattTable.X[1,:,:])
const TauMean = 1/MaxRate
const ParPerBin = Int(NumPar/NumBin)
const TotScatt= ScattCount.G + ScattCount.L + ScattCount.X
#const ConstField = -0.3e6
#println(Mat.EpIV)
#println(Mat.DefPotIV)
println(BarHeight)
println(Mat.ValleySplit[:,1]," ", Mat.ValleySplit[:,200])
ParPos, ParWvk, ParTau, ParStat, ParValley, ParPhysValley = InitializeParticles(Mat, LocNode)
OutFile = InitialiazeIO()
Totals = InitializeTotals()
Dist, WeightedCharge = InitializeDist()
println("Material and scattering tables succesfully initialized!")
println("Paricle information succesfully initialized!")
println("I/O files succesfully initialized!")


########################### BEGIN SIMULATION (TIMESTEPPING)#############################################
main(ParPos, ParWvk, ParTau, ParValley, ParPhysValley, ParStat, LocNode, BinLen, Mat, Dist, Totals, OutFile, WeightedCharge, BarHeight, Elip, Contact, ScattTable, ScattCount, ScattIndex)
