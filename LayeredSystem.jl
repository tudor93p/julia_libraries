module LayeredSystem
#############################################################################
#
#		Contains the methods needed to implement the algorithm presented in 
#
#					Lima, Dusko, Lewenkopf - 	PRB 97, 165405 (2018)
#
#			"Efficient method for computing the electronic transport properties 
#								of a multiterminal system "
#
#		The two-terminal counterpart is presented in
#
#				Lewenkopf, Mucciolo - J Comput Electron 12, 203–231 (2013)
#
#						"The recursive Green’s function method for graphene"
#
#############################################################################


import TBmodel,Utils,Graph


#===========================================================================#
#
# Helper function, prepares the lead dictionary for further operations
# 			:label and :coupling must be additionally provided
#
#---------------------------------------------------------------------------#


function PrepareLead(pyLead,BridgeAtoms,HoppMatr,LeadGF)
""" 
	- pyLead = python Lattice object; Lead aligned and attached to Atoms

	- BridgeAtoms::Matrix{Float64} is non-empty only for unregular lead-atom
				iterfaces. Contains additional atoms to ensure maximum contact

	- HoppMatr::Function(atoms1,atoms2) gives the Hamiltonian

	- LeadGF::Function(E) lead GF, computed elsewhere

"""

	LeadAtoms = pyLead.PosAtoms()

	LeadIntra = HoppMatr(LeadAtoms)

	LeadInter = HoppMatr(LeadAtoms,LeadAtoms .+ pyLead.LattVect) 


	isnothing(BridgeAtoms) && return Dict(
																				
															:head => [LeadAtoms], 
															
															:intracell => [LeadIntra],

															:intercell => [LeadInter],

															:GF => E->[LeadGF(E)]

																				)

	BridgeIntra = HoppMatr(BridgeAtoms) 

	BridgeToLead = HoppMatr(BridgeAtoms,LeadAtoms)


	return Dict(

		:head => [BridgeAtoms,LeadAtoms],

		:intracell => [BridgeIntra,LeadIntra],

		:intercell => [BridgeToLead,LeadInter],

		:GF => 	E-> (g->[RGF(E,BridgeIntra,(BridgeToLead',g)),g])(LeadGF(E))

						)
end





#===========================================================================#
#
#	Groups the atoms into several layers, such that
#					- all atoms connected to leads are in layer 1
#					- for n>1, layer n contains all neighbors of atoms in layer (n-1)
#
#---------------------------------------------------------------------------#

function Distribute_Atoms(Atoms,isbond,LeadContacts)
			"""	
	Atoms::Array{Float64,2} -> positions of atoms in the scattering region

	isbond(a,b)::Function -> says whether there's a bond between a and b
										(possibly much faster than computing the Hamiltonian)

	LeadContacts

			"""


	out = Dict{Int64,Int64}()

			# -------- find layers iteratively ----------- # 

	function f(layer,candidates)

		isempty(candidates) && return layer-1

		get!.([out],candidates,layer)

		return f(layer+1, filter(setdiff(axes(Atoms,1),keys(out))) do i

							 					ai = Atoms[i,:]

									  		return any(c->isbond(Atoms[c,:],ai),candidates)

											end)
	end


	N = f(1,vcat(LeadContacts...) |> lc -> isempty(lc) ? [1] : lc)


	# -------- sanity check and return -------------------- #

	length(out) == size(Atoms,1) || error("There are unallocated atoms.")


	LayerOfAtom, IndsAtomsOfLayer = Utils.FindPartners(out,sortfirst=true)




	return Dict(

			:NrLayers=> N,

			:LayerOfAtom => LayerOfAtom,

			:IndsAtomsOfLayer => IndsAtomsOfLayer,

			:AtomsOfLayer => L->Atoms[IndsAtomsOfLayer(L),:],

						)


end


function LayerSlicer(;LayerOfAtom,IndsAtomsOfLayer,Nr_Orbitals,kwargs...)
	
	return function (name::String,index::Int64)

		name == "Layer" && return (name,index),(Colon(),)
	
		layer = LayerOfAtom(index)
	
		atoms = IndsAtomsOfLayer(layer)
	
		return ("Layer",layer),(if length(atoms)==1 Colon()
#														elseif isnothing(Nr_Orbitals) nothing
														else 
	
			TBmodel.Hamilt_indices(	vcat(1:Nr_Orbitals...),
															findfirst(isequal(index),atoms),
															Nr_Orbitals)
														end,)

	end
		
end
#===========================================================================#
#
#	Sanity check for the atom <--> layer relationship:
#			any lead should be couple to a single & terminal layer
#
#---------------------------------------------------------------------------#


function Check_AtomToLayer(LeadContacts=[];kwargs...)
	#NrLayers=nothing,LayerOfAtom=nothing,IndsAtomsOfLayer=nothing,AtomsOfLayer=nothing)
	
	for key in [:NrLayers,:LayerOfAtom,:IndsAtomsOfLayer,:AtomsOfLayer]
		
		haskey(kwargs,key) || return false

	end

	N,LayerOfAtom = kwargs[:NrLayers], kwargs[:LayerOfAtom]

	for LC in LeadContacts

		contacts = LayerOfAtom.(LC)

		!any(boundary -> all(contacts.==boundary),[1,N]) && return false

	end

	return true

end




#==========================================================================#
#
#	Combine several (disjoint) leads into a single one 
#
#---------------------------------------------------------------------------#

function Combine_Leads(leads,atoms,label)

	isempty(leads) && return nothing

	coupling(l) = l[:coupling](l[:head][1],atoms)


	if length(leads)==1

		return Dict(
			:label =>	label,
			
			:coupling => coupling(leads[1]),

			(k=>leads[1][k] for k in [:head, :intracell,:intercell,:GF])...
			
						), size.(leads[1][:intracell],1)

	end


					#	helper function
			
	f(k) = [l[k] for l in leads]

	f(k,j) = [item[min(j,end)] for item in f(k)]

#	f(k,j,E) = [item(E) for item in f(k,j)]

	nr_ucs = maximum(length.(f(:intracell)))


	NewLead = Dict(

		:label =>	label,

		:head => map(1:nr_ucs) do j vcat(f(:head,j)...) end,

		:coupling => vcat(coupling.(leads)...),

		:intracell => map(1:nr_ucs) do j Utils.BlkDiag(f(:intracell,j)) end,

		:intercell => map(1:nr_ucs) do j Utils.BlkDiag(f(:intercell,j)) end,

		:GF => function (E)

							gfs_at_E = [l[:GF](E) for l in leads]	# one single evaluation!
							
							return map(1:nr_ucs) do j 		

								Utils.BlkDiag(q[min(j,end)] for q in gfs_at_E)

							end

						end
		
					)

	subsizes = map(1:nr_ucs) do j size.(f(:intracell,j),1) end

	return NewLead,subsizes

end






#===========================================================================#
#
# Groups the user-given leads into two "virtual" leads: LeftLead & RightLead
#
#---------------------------------------------------------------------------#

function Distribute_Leads(Leads, LeadContacts; NrLayers, LayerOfAtom, AtomsOfLayer, kwargs...)

	isempty(Leads) && return Dict()

	lead_distrib = Dict(("LeftLead",1)=>[],
											("RightLead",NrLayers)=>[])

	
	for (iL,LC) in enumerate(LeadContacts)

		isempty(LC) && error("Lead does not couple to any atoms")
		
		layers = LayerOfAtom.(LC)

		for (side,n) in keys(lead_distrib)

			if all(layers .== n) 
				
				push!(lead_distrib[(side,n)],iL)
			end
		end	
	end

	if sum(length.(values(lead_distrib))) != length(Leads)

		error("A lead is not coupled to terminal layers.")

	end




	VirtLeads,LeadSizes = Dict(), Dict()

	for ((side,n),i) in filter!(p->!isempty(p.second),lead_distrib)

		VirtLeads[Symbol(side)],LeadSizes[side] = Combine_Leads(Leads[i],
																														AtomsOfLayer(n),
																														side)
	end

	
	SideOfLead, LeadsOfSide_ = Utils.FindPartners([string(Leads[i][:label])=>string(side) for ((side,n),I) in lead_distrib for i in I ], sortfirst=string)


	LeadsOfSide(s) = LeadsOfSide_(string(s))


	function slicer(name::String,index::Int64)

		name in ["LeftLead","RightLead"] && return (name,index),(Colon(),)

		side = SideOfLead(name)

		allleads = LeadsOfSide(side)

		length(allleads)==1 && return (side,index),(Colon(),)

	 	lead = findfirst(isequal(name),allleads)
	
		boundaries = cumsum([0;LeadSizes[side][min(end,index)]])
	
		return (side,index),(boundaries[lead]+1:boundaries[lead+1],)

	end

	
	return VirtLeads,Dict(	:SideOfLead => SideOfLead, 
													:LeadsOfSide => LeadsOfSide,
													:LeadSlicer => slicer
											)
end




#===========================================================================#
#
#		Brings together the (layer <--> atom) and (real lead <--> virtual lead)
#					relations into a single translator function;
#					
#			↪	 also gives the sector corresponding to the desired atom/real lead
#
#---------------------------------------------------------------------------#

function LeadLayerSlicer(;LeadSlicer=nothing,kwargs...)

	layer_slicer = LayerSlicer(;kwargs...)
														 

	slicer(name::Symbol,index::Int64,args...) = slicer(string(name),index,args...)
	slicer(name::Symbol,index::Int64) = slicer(string(name),index)
	
		

	function slicer(name::String,index::Int64,args...)

		ni1,slice1 = slicer(name,index)

		ni2,slice2 = slicer(args...)


		return (ni1...,ni2...),(slice1...,slice2...)
	
	end


	function slicer(name::String,index::Int64)

		name in ["Atom","Layer"] && return layer_slicer(name,index)
	
		return LeadSlicer(name,index)

	end

	return slicer
end










#===========================================================================#
#
#		Verify that the lead dictionary has all necessary information
#
#---------------------------------------------------------------------------#

#function CheckLead(Lead=nothing,label=nothing)
#	
#	isnothing(Lead) && return nothing
#
#	if !all(haskey.([Lead],[:intracell,:intercell,:GF,:coupling,:head]))
#		error("Incomplete lead information")
#	end
#
#	get!(Lead,:label,label)::String
#	
#	Lead[:coupling]::AbstractArray
#
#	[a::AbstractArray for k in [:intracell,:intercell] for a in Lead[k]]
#
#	typeof(Lead[:GF]) <: Function || error("Provide lead GF as a function!")
#
#  return Lead
#
#end



#===========================================================================#
#
#		Map a layered TB system to a linear graph (directed path, graph)
#				
#				Store on the graph only the unique Hamiltonians (no H.C.)
#
#---------------------------------------------------------------------------#


function LayeredSystem_toGraph(HoppMatr,NrLayers;LeftLead=nothing,RightLead=nothing)
					"""			 
	HoppMatr(unit_cell_n,unit_cell_m) = Hamiltonian between the two layers

	NrLayers -> upper bound of n,m  (lower bound is 1)

	If any leads are specified, they must contain
		-	an identifier :label => label::String
		-	coupling matrix :coupling => U::AbstractArray
		-	hopping matrix :intercell => H::AbstractArray
		-	the GF :GF => g::Function

				"""


#	LeftLead  = CheckLead(LeftLead,"LeftLead")
#	RightLead = CheckLead(RightLead,"RightLead")

	Leads = filter(!isnothing,[LeftLead,RightLead])


	NrLeadUCs = if isempty(Leads) 0 else

									maximum([length(L[:intracell]) for L in Leads])+1

							end

#	NrLeadUCs = 5

  g = Graph.MetaDiPath(NrLayers)

	Graph.set_props!(g,Dict(:NrLayers=>NrLayers,
													:UCsLeads=>NrLeadUCs,
													:LeadLabels=>[Lead[:label] for Lead in Leads]
											 ))

	for i in 1:NrLayers

		Graph.set_props!(g,i,Dict(:type=>"Layer",
														 	:name=>("Layer",i),
														 	:H=>HoppMatr(i))
														)
	
		i>1 && Graph.set_prop!(g, i-1, i, :H, HoppMatr(i-1,i))


	end


	for Lead in Leads,i in 1:NrLeadUCs

			Graph.add_vertex!(g,Dict(vcat(
							:type=>"VirtualLead",
							:name=>(Lead[:label],i),
							:H => Lead[:intracell][min(i,end)],
							(i>1 ? [] : [:GFf => Lead[:GF],])... 
					# the GF for all lead unit cells will be stored in the lead head. 
												 )))
	end

  node = Graph.node_by_prop!(g,:name)


	
	if !isnothing(RightLead)

		Graph.add_edge!(g,NrLayers,
											node(RightLead[:label],1),
											:H, RightLead[:coupling]')

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g,node(RightLead[:label],i),
													node(RightLead[:label],i+1),
													:H, RightLead[:intercell][min(i,end)]
											)
			end # intercell is given in the u direction, as here: i -> i+1 
	end

	if !isnothing(LeftLead)

		Graph.add_edge!(g,node(LeftLead[:label],1),
											1,
											:H, LeftLead[:coupling])

		for i in 1:NrLeadUCs-1

				Graph.add_edge!(g,node(LeftLead[:label],i+1),
													node(LeftLead[:label],i),
													:H, LeftLead[:intercell][min(i,end)]'
											)
		end # intercell is given in the u direction: opposite to i+1->i 
	end


	return g

end














































#############################################################################
end
