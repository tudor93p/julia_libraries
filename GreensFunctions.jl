module GreensFunctions
#############################################################################

import LinearAlgebra; const LA = LinearAlgebra
##import SparseArrays; const SpA = SparseArrays

import Utils,Operators,Graph,LayeredSystem


#===========================================================================#
#
# Green's function using brute-force inversion
#
#---------------------------------------------------------------------------#



""" gives the GF at energy E of an isolated system S with Hamiltonian H """

GF(E::Number,H::AbstractMatrix) = inv(Array(E*LA.I - H))


""" the system S is now coupled to another system S'
 									such that the self-energy is SE			"""

GF(E::Number,H::AbstractMatrix,SE::AbstractMatrix) = GF(E,H+SE)


"""
	- the system S is now coupled to n = length(args) 
				other systems {S_1, S_2, ..., S_n} with known GFs g_n

	- args is a series of tuples. Element i cooresponds to system S_i:
				args[i] = (	H( S_i -> S ), g_i )  ---> as required by SelfEn
"""

GF(E::Number,H::AbstractMatrix,args...) = GF(E,H+SelfEn(args...))


""" the Hamiltonian of S is null """
 
GF(E::Number,args...) = GF(E, SelfEn(args...))


""" the GF of isolated S is already known, g = inv(E-H)
							with or without coupling, the total GF is, respectively:"""

GF(g::AbstractMatrix,args...) = GF(0.0,-inv(g),SelfEn(args...))

GF(g::AbstractMatrix) = g


#===========================================================================#
#
#	Self energy due to coupling to a system with known GF
#
#---------------------------------------------------------------------------#

	"""
	Self energy associated with the coupling to another system

	A := present system; B:= target system with GF g
	
	U = HoppingMatrix(B -> A)

	"""

SelfEn(U::AbstractMatrix,g::AbstractMatrix) =  U'*g*U

SelfEn((U,g)) = SelfEn(U,g)
			# if one tuple is given
			
SelfEn(args...) = sum(SelfEn.(args))
			#	if more than one tuple is give




#===========================================================================#
#
#	Helper functions for dealing with the physical system stored in the graph
#
#---------------------------------------------------------------------------#


function GraphLayeredSystem_Utils(g)

	node = Graph.node_by_prop!(g,:name)
	
	f() = ()
	
	f(n::Int64,args...) = (n,f(args...)...)
	
	f(T::String,I::Int64,args...) = (node(T,I),f(args...)...)
	
	
	
	get_prop(p) = Graph.get_prop(g,p)
	get_prop(p,ns...) = Graph.get_prop(g,ns...,p)
	


	islead(n) = occursin("Lead",get_prop(:type,n))

	islayer(n) = get_prop(:type,n)=="Layer"


	
	
	function H(args...)

		ns = f(args...)
		
		(length(ns) == 1 || ns[1]==ns[2]) && return get_prop(:H,ns[1])
	
		Graph.has_prop(g,ns...,:H) && return get_prop(:H,ns...)
	
		rns = reverse(ns)
	
		Graph.has_prop(g,rns...,:H) && return get_prop(:H,rns...)'

		return zeros(size(get_prop(:H,ns[1]),1),size(get_prop(:H,ns[2]),1))

	end
	
	
	function setEnergy_LeadGF(g_,Energy)
	
		for ll in get_prop(:LeadLabels)
			
			gfs = get_prop(:GFf,f(ll,1)[1])(Energy)	
		
			for uc in 1:get_prop(:UCsLeads)

				Graph.set_prop!(g_,f(ll,uc)[1],:GF,gfs[min(uc,end)])
			end
		end
	
		return n -> Graph.get_prop(g_,n,:GF)
	end
	
	
	
	function left(n)
	
	  ns = Graph.in_neighbors(g,n)
	
		return isempty(ns) ? nothing : ns[1]
	
	end
	
	function right(n)
	
	  ns = Graph.out_neighbors(g,n)
	
		return isempty(ns) ? nothing : ns[1]
	
	end
	
	
	function next(n::Int64,dir::String) 
	
		dir == "right" && return right(n)
	
		dir == "left" && return left(n)
	
		error("dir should be left or right")
	
	
	end
	


	function meets_layer(n,dir::String="both")::Bool
	
		isnothing(n) && return false

		(dir == "both" || islayer(n)) && return true
	
		return meets_layer(next(n,dir),dir)
	
	end
	
	
	
	
	
	function bring_closer(n,m) 
	
		isempty(Graph.FindPath(g,n,m)) &&	return false,left(n),right(m)
	
		return true,right(n),left(m)
	
	end
	
	function lead_extends(n,dir)
	
		return islead(n) && (dir == "both" || !meets_layer(n,dir))
	
	end



	return 	setEnergy_LeadGF,
					islead,
					meets_layer,
					lead_extends,
					H,
					next,
					bring_closer,
					f


end	


#===========================================================================#
#
# 	Green's function using decimation method for at most two leads.
#
#		Implementation as explained by Lewenkopf+Mucciolo in
#
#		"The recursive Green’s function method for graphene"
#
#				J Comput Electron (2013) 12: 203–231
#
#---------------------------------------------------------------------------#

function GF_Decimation(HoppMatr,NrLayers;LeftLead=nothing,RightLead=nothing;translate=nothing)

	g = LayeredSystem.LayeredSystem_toGraph(HoppMatr,NrLayers;LeftLead=LeftLead,RightLead=RightLead)

	return  Energy -> GF_Decimation_fromGraph(Energy,g,translate=translate)

end



function GF_Decimation_fromGraph(Energy,g,translate=nothing)

	setEnergy_LeadGF,
	islead,
	meets_layer,
	lead_extends,
	H,
	next,
	bring_closer,
	node	= GraphLayeredSystem_Utils(g)



	SemiInfLeadGF = setEnergy_LeadGF(g,Energy)


	# ----- the dictionary where the values are stored -------- #
	
	Gd = Dict{Tuple{Int64,Int64,String},Array{Complex{Float64},2}}()



	function G(n,m,dir) 
	
		n==m && islead(n) && !meets_layer(n,dir) && return SemiInfLeadGF(n)

							#	this should not be stored, already stored in the graph!

		n==m && return get!(Gd,(n,m,dir)) do 
																			Gnn(n,dir) end

		return get!(Gd,(n,m,dir)) do 
																Gnm(n,m,dir) end

	end


	# ----- methods to calculate the GF ------- #	

	function Gnn(n,dir)
	
		lead_extends(n,dir) && return GF(SemiInfLeadGF(n),
																			coupling_toSystem(n,dir)...)
	
	 	return GF( Energy, H(n), coupling(n,dir)... )
	
	end

	function Gnm(n,m,dir)

		isnleft, n1, m1  = bring_closer(n,m)
	
		for (d,test) in zip(["left","right"],[isnleft,!isnleft])
	
			if dir in ["both",d] 
	
				test && return G(n,m1,d)*H(m1,m)*G(m,m,dir)
	
				return G(n,n,dir)*H(n,n1)*G(n1,m,d)
	
			end
		end

	end



	# ---- coupling a layer left or right ----- #

	function coupling(src::Int64,dir::String)
	
		return filter(!isempty,map(["right","left"]) do d
	
										dst = (dir in [d,"both"]) ? next(src,d) : nothing
	
										return isnothing(dst) ? () : (H(dst,src), G(dst,dst,d))
	
									end)
	end
	
	
	# ------- coupling for a lead ---------- #

	
	function coupling_toSystem(src,d)
	
		for dir in ["right","left"]
	
			meets_layer(src,dir) && return coupling(src,dir)
	
		end
	
	end
	

	return if	isnothing(translate)
		
						function out1(name1::String,index1::Int64,
												 	name2::String=name1,index2::Int64=index1;
																							dir::String="both")
				
							return G(node(name1,index1,name2,index2)...,dir)
				
						end

				
				else


					function out2(name1::String,index1::Int64,
											 	name2::String=name1,index2::Int64=index1;
																										dir::String="both")
						
						n1i1n2i2,slice = translate(name1,index1,name2,index2) 

						return G(node(n1i1n2i2...)...,dir)[slice...]

					end


				end

end




#function GF_Decimation(Energy,H_intracell,H_intercell,NrLayers,gRightLeft,include_leads=2)
#				"""
#	Energy -> value of the energy at which the GF is evaluated
#
#	Inhomogeneous hopping:
#
#		H_intracell(n) -> intra-layer Hamiltonian of layer n
#
#		H_intercell(n,m) -> inter-layer hopping between layer n and layer m
#
#	NrLayers -> upper bound of n,m  (lower bound is 1)
#
#	gRightLeft = gR,gL -> functions which give the GF of the semi-infinite leads
#
#				"""
#
#	
#	function coupling(source::Int64,dest::Int64)
#
#		dest > source && return (H_intercell(dest,source),GR(dest,dest))
#
#		dest < source && return (H_intercell(dest,source),GL(dest,dest))
#
#
#	end
#
#	function G_lead_layer(lead,m,getGF)
#
#		lead > NrLayers && return GR(lead,m+1)*H_intercell(m+1,m)*getGF(m,m)
#
#		lead < 1 && return GL(lead,m-1)*H_intercell(m-1,m)*getGF(m,m)
#
#
#	end
#
#
#
#	gRs,gLs = [g(Energy) for g in gRightLeft]
#	
#
#	# methods to calculate #
#	
#	
#	function GL_(n,m)
#
#		n==m && return	GF(Energy,H_intracell(n),coupling(n,n-1))
#
#		n==0<m && return	G_lead_layer(n,m,GL)
#
#		return nothing
#
#	end
#
#	function GR_(n,m)  
#
#		n==m && return GF(Energy,H_intracell(m),coupling(m,m+1))
#
#		n==NrLayers+1>m && return G_lead_layer(n,m,GR)
#
#		return nothing
#	end
#
#
#
#	# combine sweeps into actual GFs 
#
#	function G_(n,m=n)
#
#		n==m>NrLayers && return GF(GR(n,n),coupling(n,n-1))
#
#		n==m<1 && return GF(GL(n,n),coupling(n,n+1))
#
#		n==m && return GF(Energy,H_intracell(n),coupling(n,n-1),coupling(n,n+1))
#
#		(( n<1 ) | (n>NrLayers )) && return G_lead_layer(n,m,G)
#
#		abs(n-m)==1 && return G(n)*H_intercell(n,m)*(m>n ? GR : GL)(m,m)
#
#		n in [0,NrLayers+1] && return G_lead_layer(n,m,G)
#
#	end
#
#
#	# calculate *and* store #
#
#	GLd,GRd,Gd = Dict(),Dict(),Dict()
#
#	GL(n,m) = n==m<1 ? gLs[min(end,1-n)] : get!(GLd,(n,m)) do  
#																						GL_(n,m) end
#
#	GR(n,m) = n==m>NrLayers ? gRs[min(end,n-NrLayers)] : get!(GRd,(n,m)) do 
#																						GR_(n,m) end
#
#
#	G(n,m=n) = get!(Gd,(n,m)) do 
#																G_(n,m) end
#
#
#
#	return GL,GR,G
#
#end
#

#	for (n,m) in zip(1:NrLayers+1+include_leads,NrLayers:-1:0-include_leads)

#	  GLd[(n,n)] = GF(Energy,h(n)+U(n,n-1)*GL(n-1,n-1)*U(n-1,n))

#		GLd[(n,n)] = GF(Energy,h(n),(U(n-1,n),GL(n-1,n-1)))

#	  GLd[(0,n)] = GL(0,n-1)*U(n-1,n)*GL(n,n)

#	  GRd[(m,m)] = GF(Energy, h(m) + U(m,m+1)*GR(m+1,m+1)*U(m+1,m))

#		GRd[(m,m)] = GF(Energy,h(m),(U(m+1,m),GR(m+1,m+1)))

#	  GRd[(NrLayers+1,m)] = GR(NrLayers+1,m+1)*U(m+1,m)*GR(m,m)

#	end

#


#===========================================================================#
#
# Local DOS from a Green's function given as a matrix
#
#---------------------------------------------------------------------------#

function LDOS(Gr;Op=[1],kwargs...)

  return Operators.Trace_Orbitals(-1/pi*imag(LA.diag(Gr));Op=Op,kwargs...)

end


#===========================================================================#
#
# Local DOS from a Green's function given as a function of layers
#
#---------------------------------------------------------------------------#

function LDOS_Decimation(GD,NrLayers,indsLayer;Op=[1],LeftLead=nothing,RightLead=nothing)

# corresponding to vcat(Atoms,LeftLead...,RightLead...)


	LNames,LSizes = begin

		local Ls = [LeftLead,RightLead]

		local I = findall(!isnothing,Ls)

		["LeftLead","RightLead"][I],[size.(Ls[i],1) for i in I]

	end

	cumLSizes = Utils.recursive_cumsum(LSizes)

	indsLayers = indsLayer.(1:NrLayers)


	ldos = zeros(mapreduce(length,+,indsLayers) + cumLSizes[end][end])


	for (L,inds) in enumerate(indsLayers)

		ldos[inds] = LDOS(GD("Layer",L),Op=Op,nr_at=length(inds))

	end


	for (LN,LS,cLS) in zip(LNames,LSizes,cumLSizes)
																			 
		for (uc,(nr_at,cls)) in enumerate(zip(LS,cLS))

			ldos[cls-nr_at+1:cls] = LDOS(GD(LN,uc),Op=Op,nr_at=nr_at)

		end
	end

	return ldos
	
end


#===========================================================================#
#
# Reconstruct the full GF matrix 
#
#---------------------------------------------------------------------------#


#function fullGF_fromGD(GD,NrLayers;iLayer,nr_orb,nr_at=nothing,include_leads=false)
#
#  iterLayers = include_leads ? collect(0:NrLayers+1) : collect(1:NrLayers)
#
#  if isnothing(nr_at)
#    nr_at = sum(length.(iLayer.(iterLayers)))
#  end
#
#  out = SpA.spzeros(Complex{Float64},nr_at*nr_orb,nr_at*nr_orb)
#
#  i(n) = TBmodel.Hamilt_indices_all(1:nr_orb,iLayer(n),nr_orb,flat=true) 
#
#  for n in iterLayers
#
#    for (a,b) in [(n,n),(n-1,n),(n,n+1)]
#      
#      if (a in iterLayers) && (b in iterLayers)
#
#        out[i(a),i(b)] = GD(a,b)
#
#      end
#
#    end
#
#  end
#
#  return out 
#
#end
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#








#===========================================================================#
#
#		Bulk, left-surface and right-surface GF, as presented in 
#
#				M. P. Lopez Sancho, J. M. Lopez Sancho, J. Rubio
#
#					J. Phys. F: Met. Phys. 15, 851-858 (1985) 
#
#			"Highly convergent schemes for the calculation of bulk and
#																surface Green functions"
#
#---------------------------------------------------------------------------#


function GF_SanchoRubio(Energy,H_intracell,H_intercell;
												target="bulk+-",max_iter=50,tol=1e-8)

				# The system is assumed homogeneous
				# in general, H(cell_i,cell_j) and H(cell_i)


	# ------ understanding target ------ #
	
	output_keys = ["bulk","+","-"]

	input_keys = [["bulk"],["+","plus","positive"],["-","minus","negative"]]

	desired_keys = [any(occursin.(k,[string(target)])) for k in input_keys]

	if !any(desired_keys) 
	  println("Target not understood. Returning all.")
		desired_keys = [true,true,true]
	else
		output_keys = output_keys[desired_keys]
	end


	# ------ error calculation and convergence criterion ------ #

  Errors = zeros(max_iter)

	function updateErrors(iter,alpha,beta)

    Errors[iter] = LA.norm(alpha) + LA.norm(beta)

  end

  function converged(iter)

		if iter>=3 && all(Errors[iter-3+1:iter] .< tol)

			return true

		elseif iter==max_iter 

			println("\nSancho-Rubio algorithm did not converge after ",max_iter," iterations.\n",Errors)

			return true
			
		end

		return false

	end


	# ------ performing the iterations of the algorithm ----- # 


  function SR_iter(alpha=Array(H_intercell),
									 beta=alpha',
									 epsBulk=Array(H_intracell),
									 epsSurf=epsBulk,
									 epsDualSurf=epsBulk,
									 i=1)

		updateErrors(i,alpha,beta)

		converged(i) && return GF.(Energy,[epsBulk,epsSurf,epsDualSurf][desired_keys])
			
    gBulk = GF(Energy,epsBulk)
																		# auxiliary vars
    agb,bga = alpha*gBulk*beta, beta*gBulk*alpha

    return SR_iter( alpha*gBulk*alpha,
										beta*gBulk*beta,
										epsBulk + agb + bga,
										epsSurf + agb,
										epsDualSurf + bga,
										i+1
									)
  end 

	# -------- returning the desired target GF --------- #


	length(output_keys)==1 && return SR_iter()[1]
  
	return Dict(zip(output_keys, SR_iter()))



end








#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





























































#############################################################################
end
