module BandStructure
#############################################################################

using Distributed

import DelimitedFiles; const DlmF = DelimitedFiles
import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
import Arpack

import Algebra,Utils



#===========================================================================#
#
# Decide between dense and sparse diagonalization
#
#---------------------------------------------------------------------------#

function get_eigen(H,evect=false;tol=1e-8,nr_bands=nothing,sigma=tol/10)


  if isnothing(nr_bands)

    !evect && return (k) -> (LA.eigvals(LA.Hermitian(Array(H(k)))),nothing)
 
    
    return (k)-> LA.eigen(LA.Hermitian(Array(H(k)))) |> e->(e.values,e.vectors)
 
  end

  sigma = isnothing(sigma) ? tol/10 : sigma

  Areigs(k) = Arpack.eigs(H(k),nev=nr_bands,sigma=sigma,tol=tol,ritzvec=evect)


	evect && return (k) -> Areigs(k) |> function (e)


													order = sortperm(real(e[1]))

#													println("\n",real(e[1])[order],"\n")
														
													return (e[1][order], e[2][:,order])


												end

	return (k) -> (sort(Areigs(k)[1],by=real),nothing)


end





#===========================================================================#
#
# Evals and operators on evects
#
#---------------------------------------------------------------------------#

function get_nrbands(nr_bands,H0,lim=0.1)

  nH = size(H0,1)

  !isnothing(nr_bands) && nr_bands < nH*lim && return nr_bands

  return nothing

end



function Diagonalize(H, kPoints, filename=nothing; 
										 kLabels=nothing,
										 kTicks=[0],
										 filemethod="new", parallel=false, operators=[[],[]],
										 tol=1e-8,  nr_bands=nothing, sigma=tol/10)

  k0 = kPoints[1,:]

  eig = get_eigen(H, !isempty(operators[1]); 
									tol=tol, nr_bands=get_nrbands(nr_bands,H(k0)), sigma=sigma)

  #k_label, k, [En, Op1, Op2, ...]

	psi0 = eig(k0)[2] # [1] is the eigenvalue; wfs are on columns


  bounds = cumsum([[0,1];[size(Op(psi0[:,1:1],k=k0),2) 
																							for Op in operators[2]]])


  result = vcat(
    (parallel ? pmap : map)(eachrow(kPoints)) do k

      E,P = eig(k)

      return hcat(E,[Op(P,k=k) for Op in operators[2]]...)

    end...)




  new_item, outdict = Utils.Write_NamesVals(filename, filemethod,
																		["Energy";operators[1]], result, bounds)





	if size(kPoints,1)==1
	
		return new_item("kLabels", 
										reshape(range(0,1,length=size(result,1)),:,1), outdict)

	end	

	kLabels = Utils.Assign_Value(kLabels, range(0,1,length=size(kPoints,1)))

	kLabels = reshape(repeat(kLabels,
													 inner=div(size(result,1), length(kLabels))),:,1)


	# kPoints should also be repeated!


	for (name,value) in zip( ["kLabels",  "kTicks"], #,"kPoints",], 
														[kLabels, kTicks] )#, kPoints]	)

		outdict = new_item(name, value, outdict)

	end

  return outdict


end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Read_Bands(filename,colors=[])

  return Utils.Read_NamesVals(filename,[["Energy","kLabels"];colors])

end

#===========================================================================#
#
# Time evolution of a state. psi(t) = e^iHt psi0
#
#---------------------------------------------------------------------------#

function TimeEvolution_expH(psi0,es,vs,ts,filename,operators=[[],[]];tol=1e-7)

  V = reshape(vs'*psi0,1,:).*vs

  states = V * exp.(1im*es.*reshape(ts,1,:))

  states[:,1] = psi0  

  weights(P;kwargs...) =  abs.(transpose(P)*conj(vs))



  Write, = Utils.Write_NamesVals(filename,tol=tol)

  Write("Time",ts)
  Write("Eigenvals",es)

  for (name,oper) in zip([operators[1];"Coef"],[operators[2];weights])

    Write(name,oper(states))

  end


end




#===========================================================================#
#
# Analyzes the energy spectrum find a reasonable energy step
#
#---------------------------------------------------------------------------#


function Analyze_EnSplitting(En,kPoints)


  En = reshape(En,:,size(kPoints,1))

#  En[band_index,k_index]
# Algebra.OuterDist(E,E) |> d -> d[LA.triu(trues(size(d)),1)] |> d -> sort(d[d .> tol])


  bandsep = abs.(diff(En,dims=1))[:] #all band splittings

#  mesh = sort(vcat(meshsplit.(eachrow(En))...)) # all splittings En(k)->En(k')



  if size(kPoints,1) ==1

    closestk = 1

  else

    closestk = hcat(sortperm.(eachcol(Algebra.OuterDist(kPoints,kPoints)))...)[2,:]
  end

  

  p=2 

  D(i) = minimum(abs.(En[i:i,:] .- En[i-p.<=axes(En,1).<=i+p,closestk]),dims=1)
     #  splittings of level i|En[i,k] - En[closest i',closest k]| for all ks

  meshsplit = hcat(D.(axes(En,1))...)[:]# |> x -> sum(x)/length(x)

#  sort(x)[div(length(x),10)]

  return [sum(x)/length(x) for x in [meshsplit,bandsep]]
# DOS: delta should be smaller than the typical band separation (bandsep), but larger than the energy splitting due to the finite k-mesh (meshsplit)

end


#===========================================================================#
#
# Read the DOS
#
#---------------------------------------------------------------------------#


function Read_DOS(filename;weights="Lorentzian",LDOS=false)#::String)

  En,kPts = DlmF.readdlm.(filename.(["Energy","kPoints"]))

  scales = vcat(minimum(En),maximum(En),Analyze_EnSplitting(En,kPts)...)

  LDOS = LDOS && DlmF.readdlm(filename("LDOS"))

  return scales,fDOS(En,LDOS,weights)

end



function fDOS(En,LDOS=false;weights="Lorentzian",get_weights=false)

  weights = Algebra.get_Distribution(weights)


  function f1(E,delta)

    W = weights(E,En,delta)

    return sum(W,dims=2)

  end


  function f2(E,delta)

    W = weights(E,En,delta)
#    W = weights(E,transpose(En),delta)



    return W*LDOS./(sum(W,dims=2) .+ 1e-20)

  end

  function f3(E,delta)

    W = weights(E,En,delta)

    return sum(W,dims=2),W

  end


  function f4(E,delta)

    W = weights(E,En,delta)
#    W = weights(E,transpose(En),delta)

    return W*LDOS./(sum(W,dims=2) .+ 1e-20),W

  end


  if !get_weights

    !isa(LDOS,AbstractArray) && return f1

    return f1,f2

  else

    !isa(LDOS,AbstractArray) && return f3

    return f3,f4

  end

end













#===========================================================================#
#
# Wilson loop operators
#
#---------------------------------------------------------------------------#


#function WLO(H,kPoints,filename=nothing;kLabels=axes(kPoints,1),parallel=false,tol=1e-8,operators=[[],[]],nr_bands=nothing,sigma=tol/10)

function occupied(E,P)

  return P[:,sortperm(E) .<= div(size(E,1),2)]

end

function WLO_(psi,kPoints)

  Psi = map(psi,eachrow(kPoints[1:end-1,:]))

  return Psi[1]'*mapreduce(P->P*P',*,Psi[2:end])*Psi[1]

end

function WLO(H,kPoints;tol=1e-8,nr_bands=nothing,sigma=tol/10)

  eig = get_eigen(H,true;tol=tol,nr_bands=nr_bands,sigma=sigma)

  return WLO_(k->occupied(eig(k)...),kPoints)

end


function WLO_Spectrum(W,kPoints,filename=nothing;kLabels=axes(kPoints,1),parallel=false)

  evals(k) = angle.(LA.eigvals(W(k)))/2pi

  E = vcat((parallel ? pmap : map)(evals,eachrow(kPoints))...)


  function write_(name,matrix)

    
    if !isnothing(filename)
      open(filename(name),"w") do fout

        DlmF.writedlm(fout,ndims(matrix) == 1 ? reshape(matrix,:,1) : matrix)

      end
    end 

    return matrix
  end


  outdict = Dict()

  outdict["Energy"] = write_("Energy",E)

  labels = repeat(kLabels,inner=div(size(E,1),length(kLabels)))

  outdict["kLabels"] = write_("kLabels",labels)

  outdict["kPoints"] = write_("kPoints",kPoints)

  return outdict




end





function WLO2(H,W,kPoints;tol=1e-8,nr_bands=nothing,sigma=tol/10)

  eigH = get_eigen(H,true;tol=tol,nr_bands=nr_bands,sigma=sigma)

  eigWLO(k) = LA.eigen(W(k)) |> e-> (angle.(e.values),e.vectors)

  w(k) = occupied(eigH(k)...)*occupied(eigWLO(k)...)

  return WLO_(w,kPoints)

end




#===========================================================================#
#
# 
#
#---------------------------------------------------------------------------#























#############################################################################

end


