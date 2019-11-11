
cd("c:\\users\\mgjo\\Box Sync\\Documents\\Research projects\\juliaprogram\\constlab")


using Tensors,PGFPlotsX,LaTeXStrings
using ForwardDiff
use DelimetedFiles

include("materials.jl");

model=2

if model==1
#model=1: linear isotropic elasticity
  nstate=1 #cannot set to 0
  npara  =2
  Emod=1000.0; nu=0.3
  mpara=[Emod,nu]
  tol=1e-3
elseif model==2  #von mises, nonlinear mixed hardening

  nstate=10
  npara  =7
  Emod=200000.0; nu=0.3; Sy=100; Hiso=Emod/2; kappainf=30; Hkin=Emod/2; alphainf=50;
  mpara=[Emod,nu,Sy,Hiso,kappainf,Hkin,alphainf]
  tol=1e-3

elseif model==3  #maxwell viscoelasticity

  nstate=1 #cannot set to 0
  npara  =4
  Emod=200000.0; nu=0.3; kappa=1*10^5; mu=1*10^5
  mpara=[Emod,nu,kappa,mu]
  tol=1e-3

end

stateold=zeros(nstate)





#Define load control
#strainctrlv=[1,0,0,0,0,0,0,0,0]
#strain control, uniaxial stress
unkn=[2,3,4,5,6,7,8,9];
notunkn=[1];
#pure strain control
#unkn=[];
#notunkn=[1,2,3,4,5,6,7,8,9];

nunkn=size(unkn,1);
nnotunkn=size(notunkn,1);
stress_unbalance=zeros(nunkn);

tol_unb=1e-3;


#Define time stepping
nsteps=300;
tend=10.; dt=tend/(nsteps-1);
timehist=0.:dt:tend;

#Strain stepping
ϵ_max =0.005; dϵ=ϵ_max/(nsteps-1)
ϵ_hist=0.:dϵ:ϵ_max;

#Stress stepping
σ_max =500.; dσ=(σ_max/(nsteps-1))
σ_hist=0.:dσ:σ_max;
σ_hist=σ_hist.*0.


aa=zeros(nsteps,13)
for i=1:nsteps
   aa[i,:]=[timehist[i],ϵ_hist[i],0.,0.,0.,0.,0.,σ_hist[i],0.,0.,0.,0.,0.]
end
writedlm("loadhist.txt",aa)

ϵ_out = zeros(nsteps) #Float64[]
σ_out = zeros(nsteps) #Float64[]
σ_out2 = zeros(nsteps) #Float64[]



#
ϵ_out,σ_out,σ_out2=time_hist(nsteps,tol_unb,model,tol,σ_hist,ϵ_hist,nstate,stateold,dt,npara,mpara,unkn,notunkn)


println("stress=",σ_out)




#μ = Emod

#=@pgf Axis({thick, no_marks, legend_pos={north_west}, xmajorgrids, ymajorgrids, ylabel=L"\bullet / E", xlabel=L"\epsilon_{11}",width="12cm", height="8cm"},
    Plot(Coordinates(ϵ_hist, σ_out./Emod)),
    LegendEntry(L"\sigma_{11}"),
    Plot(Coordinates(ϵ_out, σ_out2./Emod)),
    LegendEntry(L"\sigma_{22}")
)
=#

@pgf Axis({thick, no_marks, legend_pos={north_west}, xmajorgrids, ymajorgrids, ylabel=L"{\rm [MPa]}", xlabel=L"\epsilon_{11}",width="12cm", height="8cm"},
    Plot(Coordinates(ϵ_hist, σ_out)),
    LegendEntry(L"\sigma_{11}"),
    Plot(Coordinates(ϵ_out, σ_out2)),
    LegendEntry(L"\sigma_{22}")
)
