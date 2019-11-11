


#######
### time stepping function

function time_hist(nsteps,tol_unb,model,tol,σ_hist,ϵ_hist,nstate,stateold,dt,npara,mpara,unkn,notunkn)


ϵ_out=zeros(nsteps); σ_out=zeros(nsteps); σ_out2=zeros(nsteps);
strainv=zeros(9); stressv=zeros(9); dstrain=zeros(9); strainold=zeros(9); stressold=zeros(9);
stress_out=zeros(9); strain_out=zeros(9); statenew=zeros(nstate); lconv=1
#println("ppp",notunkn)
#println("qqq",unkn)
nunkn=size(unkn,1)

kalle=zeros(9)

for n=1:nsteps


  stressv[unkn]=σ_hist[n]*ones(nunkn);
  strainv[notunkn]=ϵ_hist[n]*ones(nnotunkn);


  stress_out,strain_out,statenew,lconv=solve_stress_iteration(tol_unb,model,tol,strainold,stressold,nstate,stateold,strainv,dt,npara,mpara,stressv,unkn)

  #push!(ϵ_out, c)
  ϵ_out[n]=strain_out[1]
  #push!(σ_out, stress_out[1])
  σ_out[n]=stress_out[1]
  #push!(σ_out2,stress_out[2])
  σ_out2[n]=stress_out[2]

  stressv,strainv,stressold,strainold,stateold=update_variables(stress_out,strain_out,statenew)



end


#end time stepping function

return ϵ_out,σ_out,σ_out2
end


function update_variables(stress_out,strain_out,statenew)
  stressv=stress_out*1.; strainv=strain_out*1.;
  stressold=stress_out*1.; strainold=strain_out*1.;
  stateold=statenew*1.

  return stressv,strainv,stressold,strainold,stateold

end


####### solve if some stress components are prescribed

function solve_stress_iteration(tol_unb,model,tol,strainold,stressold,nstate,stateold,strainv,dt,npara,mpara,stressv,unkn)



#(tol_unb::Float64,model::Int32,tol::Float64,strainold::Float64,stressold::Float64,stateold::Float64,strainv::Float64,dt::Float64,mpara::Float64,stressv::Float64,unkn::Int32)

  res_unb=1e10; niter=0
  stress_out=zeros(9); strain_out=zeros(9); statenew=zeros(nstate); lconv=1;
  #Array{Float64}(undef, 2)


  #
  while res_unb>tol_unb

     dstrain=strainv-strainold


     stress_out, statenew, stiffv, lconv = constitutive_driver(model,tol,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)


     if (size(unkn,1)>0)

         stress_unbalance=stress_out[unkn]-stressv[unkn]
         res_unb=norm(stress_unbalance)
         niter=niter+1
         println("unb=  ",res_unb,"   niter=",niter) #,"   dstrain=",dstrain)
         #println("s11= ",stress_out[1])

         #Newton update
         strainv[unkn]=strainv[unkn]-inv(stiffv[unkn,unkn])*stress_unbalance

     else
         res_unb=0.e0
     end


    #
  end #while

    strain_out=strainv*1.

    return stress_out,strain_out,statenew,lconv


end

#############################



struct hookeisopar
    #constant, not updated for every element
    E::Float64
    ν::Float64
    C::Matrix{Float64}
end

#=function Hookeiso(E::Float64, ν::Float64, plainstrain::Bool)
        μ=E/(2*(1+ν))
        λ=E*ν/((1+ν)*(1-2*ν))
        if plainstrain
                C=[     2μ+λ    λ       0;
                        λ       2μ+λ  0;
                        0       0       μ]
        else
                C=2μ/(2μ+λ)*[   2(μ+λ)  λ       0;
                                λ       2(μ+λ)  0;
                                0       0       (2μ+λ)/2]
        end
        Hookeiso(E, ν, C)
end =#




#not used for 2D case
function hookeiso_stiff(E::AbstractFloat, ν::AbstractFloat)::Tensor{4,3,Float64}
        δ(i,j) = i==j ? 1.0 : 0.0
        μ=E/(2*(1+ν))
        λ=E*ν/((1+ν)*(1-2*ν))
        #C=Tensor{4,3,Float64}((i,j,k,l)->μ*(δ(i,l)*δ(j,k)+δ(i,k)*δ(j,l)) + λ*δ(i,j)*δ(k,l))
        C=Tensor{4,3,Float64}((i,j,k,l)->2*μ*(δ(i,k)*δ(j,l)) + λ*δ(i,j)*δ(k,l))
        return C
end

function hookeiso_stiff_GK(G::AbstractFloat, K::AbstractFloat)::Tensor{4,3,Float64}
        δ(i,j) = i==j ? 1.0 : 0.0
        μ=G
        λ=K-2*G/3
        C=Tensor{4,3,Float64}((i,j,k,l)->2*μ*(δ(i,k)*δ(j,l)) + λ*δ(i,j)*δ(k,l))
        return C
end


###############

function constitutive_driver(model,tol,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)

stressout=zeros(9); statenew=zeros(nstate); stiffv=zeros(9,9); lconv=1;

if model==1

  stiffv=tovoigt(hookeiso_stiff(mpara[1],mpara[2]))
  stress_out=stressold+stiffv*dstrain
  statenew=stateold*1.
elseif model==2

  stiffv=tovoigt(hookeiso_stiff(mpara[1],mpara[2]))
  stress_out=stressold+stiffv*dstrain
  #println(stress_out,stressold,"****")
  statenew=stateold*1.
  xinit=zeros(20)
  xinit[1:9]=stressold; xinit[10:19]=stateold; xinit[20]=0.
  stress_out,statenew,stiffv=local_problem(xinit,model,tol,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)

else model==3

  stiffv=tovoigt(hookeiso_stiff(mpara[1],mpara[2]))
  stress_out=stressold+stiffv*dstrain
  #println(stress_out,stressold,"****")
  statenew=stateold*1.
  xinit=zeros(9)
  xinit[1:9]=stressold;
  stress_out,statenew,stiffv=local_problem_ve(xinit,model,tol,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)

end



return stress_out, statenew, stiffv, lconv

end
################# solve local problem
###  Here my plasticity code was .../ M Ekh
############ For local problem
#
function extract_statevariables_model_plast_prototype(stateold)

     kappaold=stateold[1]; alphaold=fromvoigt(Tensor{2,3},stateold[2:10]);
     return kappaold,alphaold
end

function collect_statevariables_model_plast_prototype(kappa,alpha)
     state=zeros(eltype(kappa),10)
     state[1]=kappa
     state[2:10]=tovoigt(alpha)
     return state
end

function extract_unknowns_model_plast_prototype(x)

     stress=fromvoigt(Tensor{2,3},x[1:9]); kappa=x[10]; alpha=fromvoigt(Tensor{2,3},x[11:19]); dlambda=x[20]
     return stress,kappa,alpha,dlambda
end
#




################# solve local problem viscoelasticity
function local_problem_ve(xinit,model,tol,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)


nunkn=size(xinit,1)
iter=0; x=zeros(nunkn)
x=xinit
maxitr=50
stress_out=zeros(9)
statenew=zeros(nstate)

R=zeros(eltype(x), nunkn)
#dR=zeros(eltype(x), [20,20])

#elastic trial step
C=hookeiso_stiff(mpara[1],mpara[2])


Jacob=zeros(nunkn,nunkn)
while true; iter+=1

    R=local_unbalance_ve(x,model,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)

    if norm(R)< tol
      Jacob=ForwardDiff.jacobian(y->local_unbalance_ve(y,model,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara),x)
      break
    end
    if iter > maxitr
                error("Local problem did not converge aftter $maxitr iterations.")
    end
    dR=ForwardDiff.jacobian(y->local_unbalance_ve(y,model,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara),x)
#   println("***  ",ForwardDiff.jacobian(y->local_unbalance(y,model,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara),x))
    #Δx = -dR\R
    #x += Δx
    #println(dR)
    x=x-dR\R

end

invJacob=inv(Jacob)
ats=tovoigt(dcontract(fromvoigt(Tensor{4,3},invJacob[1:9,1:9]),C))


stress=fromvoigt(Tensor{2,3},x[1:9]);
statenew=stateold



stress_out=tovoigt(stress)



return stress_out,statenew,ats
end





############ Local problem
function local_unbalance_ve(x,model,strainold,stressold,nstate,stateold,dstrain,dt,npara,mpara)
     nunkn=size(x,1)
     Rr=zeros(eltype(x),nunkn)
     if model==3

      stress=fromvoigt(Tensor{2,3},x[1:9])

      C=hookeiso_stiff(mpara[1],mpara[2])
      Dv=hookeiso_stiff_GK(1/(2*mpara[3]),1/(2*mpara[4]))
      stresstrial=fromvoigt(Tensor{2,3},stressold)+dcontract(C,fromvoigt(Tensor{2,3},dstrain))
      Rstress=stress-stresstrial+dt*dcontract(dcontract(C,Dv),stress)
      Rr[1:9]=tovoigt(Rstress);

     else

       println("model not implemented ")

     end

     return Rr
end





