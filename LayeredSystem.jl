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

function DistribAtoms_toLayers(Atoms,isbond,LeadContacts)

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


	return N,Utils.FindPartners(out,sortfirst=true)



end


#===========================================================================#
#
#	Sanity check for the atom <--> layer relationship:
#			any lead should be couple to a single & terminal layer
#
#---------------------------------------------------------------------------#


function Check_AtomToLayer(N,LayerAtomRel,LeadContacts=[])
													 
	LayerOfAtom = LayerAtomRel[1]
	
	for LC in LeadContacts
	
		if !(all(LayerOfAtom.(LC).==1) || all(LayerOfAtom.(LC).==N))
	
			return false
	
		end
	end

	return true

end




#===========================================================================#
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

function Distribute_Leads(Leads, LeadContacts,
																	Atoms, N, LayerAtomRelations)
	
	isempty(Leads) && return Dict(),nothing,Dict()

	LayerOfAtom, IndsAtomsOfLayer = LayerAtomRelations

	AtomsOfLayer(l) = Atoms[IndsAtomsOfLayer(l),:]

	LeadFamilies = Dict(
											
		map(enumerate(zip(Leads,LeadContacts))) do (iL,(L,LC))
	
			isempty(LC) && error("Lead does not couple to any atoms")
		
			layers = LayerOfAtom.(LC)
	
			all(layers .== 1) && return L[:label]=>(iL,"LeftLead",1)
			
			all(layers .== N) && return L[:label]=>(iL,"RightLead",N)
		
			error("Lead is not coupled to terminal layers.")
		
		end)
	

	
	VirtualLeads = Dict(map(["LeftLead","RightLead"]) do side 

		leads = sort([k for (k,v) in LeadFamilies if v[2]==side],by=string)

		return Symbol(side) => if isempty(leads) nothing else

			Combine_Leads([Leads[LeadFamilies[l][1]] for l in leads],
										AtomsOfLayer(LeadFamilies[leads[1]][3]),
										side)

				# Tuple(Lead,SubSizes)
														end
	end)

	return Dict(k=>v[1] for (k,v) in VirtualLeads if !isnothing(v)),
					Utils.FindPartners([(String(k),v[2]) for (k,v) in LeadFamilies],
																													sortfirst=string),
					Dict(k=>v[2] for (k,v) in VirtualLeads if !isnothing(v))
				
end




#===========================================================================#
#
#		Brings together the (layer <--> atom) and (real lead <--> virtual lead)
#					relations into a single translator function;
#					
#			↪	 also gives the sector corresponding to the desired atom/real lead
#
#---------------------------------------------------------------------------#

function AtomLead_toLayerVirtLead(LayerAtomRel=nothing,
																	NrOrbitals=nothing,
																	VirtRealLeadRel=nothing,
																	LeadSizes=nothing)
									

	(isnothing(LayerAtomRel) | isnothing(NrOrbitals)) && return nothing
																	
	function translate(name::String,index::Int64,args...)

		ni1,slice1 = translate(name,index)
	
		ni2,slice2 = translate(args...)
	
		return (ni1...,ni2...),(slice1...,slice2...)
	
	end


	function translate(name::String,index::Int64)

		if name in ["Layer","LeftLead","RightLead"]
	
			return (name,index),(Colon(),)
	
		elseif name == "Atom"
	
			LayerOfAtom, IndsAtomsOfLayer = LayerAtomRel

			layer = LayerOfAtom(index)
	
			atoms = IndsAtomsOfLayer(layer)
	
			slice = if length(atoms)==1 Colon() else
			
				TBmodel.Hamilt_indices(1:NrOrbitals,findfirst(index .== atoms),NrOrbitals)
	
							end
	
			return ("Layer",layer),(slice,)
	
		
		elseif !isnothing(VirtRealLeadRel) & !isnothing(LeadSizes)
	
			SideOfLead, LeadsOfSide = VirtRealLeadRel

			side = SideOfLead(name)
			
			allleads = LeadsOfSide(side)
		
			slice = if length(allleads)==1 Colon() else
	
		    lead = findfirst(name .== allleads)
		
				boundaries = cumsum([0;LeadSizes[side][min(end,index)]])
	
				boundaries[lead]+1:boundaries[lead+1]
	
							end
	
			return (side,index),(slice,)

		else

			error("Wrong input")

		end

	end

	return translate
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

									maximum([length(L[:intracell]) for L in Leads])

							end



  g = Graph.MetaDiPath(NrLayers)

	Graph.set_props!(g,Dict(:NrLayers=>NrLayers,
													:UCsLeads=>NrLeadUCs,
													:LeadLabels=>[Lead[:label] for Lead in Leads]
											 ))

	for i in 1:NrLayers

		Graph.set_props!(g,i,Dict(#:type=>"Layer",
														 	:name=>("Layer",i),
														 	:H=>HoppMatr(i))
														)
	
		i>1 && Graph.set_prop!(g, i-1, i, :H, HoppMatr(i-1,i))


	end


	for Lead in Leads,i in 1:NrLeadUCs

			Graph.add_vertex!(g,Dict(vcat(
#							:type=>"VirtualLead",
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
													:H, RightLead[:intercell][i]
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
													:H, LeftLead[:intercell][i]'
											)
		end # intercell is given in the u direction: opposite to i+1->i 
	end


	return g

end














































#############################################################################
end
