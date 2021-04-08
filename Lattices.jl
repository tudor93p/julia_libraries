module Lattices
#############################################################################

using OrderedCollections:OrderedDict

import LinearAlgebra; const LA=LinearAlgebra
import Utils, Algebra




export Lattice



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

const TOLERANCE = 1e-5


# Latt : python object
# Properties:
#		Latt.LattVec
# Methods which modify: 
# 	Latt.Add_atoms(sublattices=nothing)
# 	Latt.Shift_Atoms(...)
# 	Latt.Reduce_Dim(...)
# Methods which don't:
#		Latt.Distances(...)
# 	Latt.Supercell(...)
# 	Latt.PosAtoms(sublattices)
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#




to_myMatrix(x::Real)::AbstractMatrix{Float64} = fill(x, 1, 1)
to_myMatrix(x::Nothing)::AbstractMatrix{Float64} = fill(0.0, 1,1)


function to_myMatrix(x::AbstractMatrix)::AbstractMatrix{Float64} 

	x

end 

function to_myMatrix(x::Utils.List)::AbstractMatrix{Float64}

	length(x)==0 && return zeros(Float64, 0, 0)

	if all(Utils.isList, x)
	
		out = zeros(Float64, length(x[1]), length(x))
	
		for (j,xj) in enumerate(x) # each column j is a vector

			out[:,j] = collect(xj)

		end 

		return out 

	end 

	all(isreal, x) && return reshape(collect(x),:,1)


end 


function to_myODict(x::Nothing)::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict()

end


function to_myODict(x::AbstractDict)::OrderedDict{Any, AbstractMatrix{Float64}}

	if x isa OrderedDict 

		valtype(x)==AbstractMatrix{Float64} && return x

		return OrderedDict(k=>to_myMatrix(v) for (k,v) in x)

	end 

	return OrderedDict(map(sort(collect(keys(x)), by=string)) do k

									 k => to_myMatrix(x[k])

								 end)

end 


function to_myODict(x::Utils.List)::OrderedDict{Any, AbstractMatrix{Float64}}

	isempty(x) && return OrderedDict()

	all(isa.(x,Pair)) && return OrderedDict(k=>to_myMatrix(v) for (k,v) in x)

	return OrderedDict("A"=>to_myMatrix(x))

end

function to_myODict(x::AbstractMatrix)::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict("A"=>to_myMatrix(x))

end



function to_myODict(x::Pair)::OrderedDict{Any, AbstractMatrix{Float64}}

	OrderedDict(x.first => to_myMatrix(x.second))

end


to_myODict(x...)::OrderedDict{Any, AbstractMatrix{Float64}} = to_myODict(x)






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

												

function CombsOfVecs(A::AbstractMatrix{<:T}, coeff::AbstractMatrix{<:V}; dim=2)::AbstractMatrix{promote_type(T,V)} where T where V

	dim==2 && return A*coeff

	dim==1 && return coeff*A

	error()

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


mutable struct Lattice

	LattVec::AbstractMatrix{Float64}
	
	Sublattices::OrderedDict{Any, AbstractMatrix{Float64}}

	Vacancies::OrderedDict{Any, AbstractMatrix{Float64}}

	function Lattice(LV, SL=nothing, VA=nothing; mode::Symbol=:cartesian)
									 
		lv, sl, va = to_myMatrix(LV), to_myODict(SL), to_myODict(VA) 
		
		for k in [keys(sl);keys(va)], odict in [sl,va]

			if haskey(odict, k)

				s = size(odict[k],1)

				if mode==:fractional

					size(lv,2)!=s && error("Dimensionality mismatch. $s coefficient(s) was(were) provided for sublattice '$k', but ",size(lv,2)," are necessary")

					odict[k] = CombsOfVecs(lv, odict[k])
				
				elseif mode==:cartesian

					size(lv,1)!=s && error("Dimensionality mismatch. The positions for sublattice '$k' are $s-dimensional, but the lattice vectors are ",size(lv,1),"-dimensional")

				else 

					error("Mode $mode not defined")

				end 


			end 

		end 

		return new(lv,sl,va)

	end 
	
end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function CombsOfVecs(latt::Lattice, coeff::AbstractMatrix{<:Real}; kwargs...)::AbstractMatrix{Float64}

	CombsOfVecs(latt.LattVec, coeff; kwargs...)

end



function CombsOfVecs((A,c)::Tuple{<:Union{Lattice, AbstractMatrix{<:Real}},
																	AbstractMatrix{<:Real}}; 
										 kwargs...)::AbstractMatrix{Float64}

	CombsOfVecs(A, c; kwargs...)

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function copy(latt::Lattice)

	Lattice(latt.LattVec,
					latt.Sublattices,
					latt.Vacancies)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function sublatt_labels(latt::Lattice; kind=:Sublattices)

	kind!=:Both && return collect(keys(getproperty(latt, kind)))

	return union(sublatt_labels(latt; kind=:Sublattices),
							 sublatt_labels(latt; kind=:Vacancies))

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function sublattices_contain(latt::Lattice, inp::Nothing=nothing; kwargs...)
	
	sublattices_contain(latt,""; kwargs...)

end 


function sublattices_contain(latt::Lattice, inp::Real; kwargs...)

	sublattices_contain(latt, string(Int(trunc(inp))); kwargs...)

end 


function sublattices_contain(latt::Lattice, inp::AbstractString=""; kwargs...)

	sublattices_contain(latt, [inp]; kwargs...)

end


function sublattices_contain(latt::Lattice, inp::Utils.List; kwargs...)
	
	# find all sublattices which contain items of inp 

	all_sl = sublatt_labels(latt; kwargs...)

	sl = filter(s->any(i->occursin(i,s), inp), all_sl) 

	return isempty(sl) ? all_sl : sl 

	#isempty(sl) && error("Wrong sublattice labels!")

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function EmptyPos(latt::Lattice)::AbstractMatrix{Float64}

	zeros(Float64, VecDim(latt), 0)

end




function PosAtoms(latt::Lattice;
									labels_contain=nothing,
									label=nothing,
									f=(sl,P)->P,
									kind=:Sublattces,
									kwargs...
									)

	SL = if !isnothing(label)
		
					in(label, sublatt_labels(latt, kind=kind)) ? [label] : []
					
					else 
		
						sublattices_contain(latt, labels_contain)

				end 




	for sl in SL 

		isempty(latt.Sublattices[sl]) && continue

		return hcat([f(s, latt.Sublattices[s]) for s in SL]...)

	end 


	return EmptyPos(latt)
	

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function VecDim(latt::Lattice)::Int

	size(latt.LattVec,1)>0 && return size(latt.LattVec,1) 

	for prop in [:Sublattices, :Vacancies], (k,v) in getproperty(latt, prop)
			
			size(v,1)>0 && return size(v,1)
	
	end 

	return 0

end 


function LattDim(latt::Lattice)::Int

	size(latt.LattVec,2)

end 


function LattVec(latt::Lattice)::AbstractMatrix{Float64}

	latt.LattVec

end 

LattVec(latt::Nothing)::Nothing = nothing


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function CombsOfVecs10(A_or_latt::T, stopstart...; dim=2) where T<:Union{Lattice,AbstractMatrix{<:Real}}

	D = if T<:Lattice 

					LattDim(A_or_latt) 
					
			elseif T<:AbstractMatrix

					size(A_or_latt, dim)

			else 

				error("Type '$T' not supported")

			end 

	return CombsOfVecs(A_or_latt,
										 Utils.vectors_of_integers(D, stopstart...; dim=dim),
										 dim=dim)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function UnitCells(latt::Lattice, stopstart...)

	CombsOfVecs10(latt, stopstart...; dim=2)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function OuterDist(args::Vararg{Any,2}; kwargs...)::AbstractArray

	Algebra.OuterDist(map(args) do X

		isa(X,Lattice) ? PosAtoms(X; kwargs...) : X

	end...; dim=2)

end 


function FlatOuterDist(args::Vararg{Any,2}; kwargs...)::AbstractArray

	Algebra.FlatOuterDist(map(args) do X

		isa(X,Lattice) ? PosAtoms(X; kwargs...) : X

	end...; dim=2)

end 







function Distances(latt::Lattice; nr_uc=2, nr_neighbors=1, with_zero=true, kwargs...)

	TOLERANCE>1e-1 && error("The tolerance is too large!")

	AtomsUC = PosAtoms(latt; kwargs...)
	
	Ds = FlatOuterDist(AtomsUC[:, 1:1].-AtomsUC, UnitCells(latt, nr_uc))

	Utils.Unique!(Ds; tol=TOLERANCE, sorted=true)

	!with_zero && deleteat!(Ds, Ds.<TOLERANCE)

	return nr_neighbors=="all" ? Ds : Ds[1:min(nr_neighbors+1,end)]
  
end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function Correct_LattCoeff((latt,n)::Tuple{Lattice, T})::AbstractMatrix{Int} where T

	Correct_LattCoeff(latt, n)

end
														

function Correct_LattCoeff(latt::Lattice, n::T)::AbstractMatrix{Int} where T

	@show n T  

	if T<:Float64 || Utils.isList(T,Float64) || T<:AbstractMatrix{<:Float64} ||(
												Utils.isList(n,Real) && any(e->isa(e,Float64),n))

		n1 = T<:AbstractMatrix ? n : vcat(n...)

		if isapprox(trunc.(n1), n1, atol=TOLERANCE) 

			return Correct_LattCoeff(latt, Int64.(trunc.(n1)))

		end

		error("n is not an integer/an array of integers!")


	elseif T<:Int64
	
		return Correct_LattCoeff(latt, LA.Diagonal(repeat([n], LattDim(latt))))


	elseif Utils.isList(n,Real) && all(e->isa(e,Int),n) 

		if length(n)==LattDim(latt)

			return Correct_LattCoeff(latt, LA.Diagonal(vcat(n)))

		elseif length(n)==1

			return Correct_LattCoeff(latt, LA.Diagonal(repeat([n[1]], LattDim(latt))))

		end 

		error("Incorrect length of the vector 'n' provided")


	elseif T<:AbstractMatrix{Int}
		
		all(size(n) .== LattDim(latt)) && return n 
			
		error("Incorrect size of the matrix 'n' provided")

	end 

	error("Type '$T' of input 'n' not understood")

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#===========================================================================#
#
# Generate vertices of a d-dimensional body based on d vectors
#
#---------------------------------------------------------------------------#

function BodyVertices_fromVectors(v::AbstractMatrix{T}; dim=2)::AbstractMatrix{T} where T

	CombsOfVecs10(v, 1, 0; dim=dim)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function is_left(P0, P1, P2)

	(P1[1] - P0[1]) * (P2[2] - P0[2]) - (P2[1] - P0[1]) * (P1[2] - P0[2])

end 


function prepare_polygon_vertices(V::AbstractMatrix; 
																	order_vertices=false, dim=2)::AbstractMatrix

	if order_vertices

		return prepare_polygon_vertices(
								sortslices(V, dims=dim, by=R->atan(R[2],R[1]));
								order_vertices=false, dim=dim)

	end 

	return cat(selectdim(V, [2,1][dim], 1:2),
					selectdim(selectdim(V, [2,1][dim], 1:2), dim, 1:1),
					dims=dim)

end 

function PointInPolygon_wn(V::AbstractMatrix; kwargs...)

	V = prepare_polygon_vertices(V; kwargs...) 
															 
	return P -> PointInPolygon_wn(P, V; kwargs..., prepare_vertices=false)

end 



function PointInPolygon_wn(P::Utils.List, V::AbstractMatrix; 
													 dim=2, prepare_vertices=true, kwargs...)

  wn = 0   # the winding number counter

  						# repeat the first vertex at the end

	if prepare_vertices

		V = prepare_polygon_vertices(V; dim=dim, kwargs...)

	end 

	v(i) = selectdim(V, dim, i)

	v(i,j) = v(i)[j]



  # loop through all edges of the polygon

	for i in 1:size(V, dim)-1  # edge from v(i) to v(i+1) 

		if v(i,2) <= P[2] # start y <= P[1]

			if v(i+1,2)>P[2] && is_left(v(i), v(i+1),P)>0 
									# an upward crossing & P left of edge
			
				wn += 1           # have a valid up intersect

			end 

			 
		else # start y > P[2] (no test needed)
			 
			if v(i+1,2)<=P[2] && is_left(v(i), v(i+1), P)<0

						
						# a downward crossing &  P right of edge

				wn -= 1           # have a valid down intersect


			end

		end 

	end 

  return wn != 0


end 




#===========================================================================#
#
# Obtain the ucs in big UC by checking which ucs are inside a polygon
#
#---------------------------------------------------------------------------#

function ucsUC_Polygon(v, mM, polygon; dim=2, asfunction=false)

	is_inside = PointInPolygon_wn(polygon, order_vertices=true, dim=dim)


	uc_candidates = Utils.vectors_of_integers(2, 	selectdim(mM, dim, 2).+1,
																								selectdim(mM, dim, 1).-1;
																								dim=dim)


	out(cand) = selectdim(cand, dim, 
												mapslices(is_inside, cand, dims=[2,1][dim])[:]
												) #|> collect


	!asfunction && return out(uc_candidates)


	shifts = sortslices(CombsOfVecs10(polygon, 1; dim=dim),
											dims=dim, by=LA.norm)

	return (
					
		size(shifts,dim),
					
		function shiftout(i)

			s = selectdim(shifts, dim, i:i)
	
			v = TOLERANCE * s / (LA.norm(s) + TOLERANCE/100)
	
			return out(v.+uc_candidates)

		end )



end 



#===========================================================================#
#
# Calculates the positons of the small unit cells inside the large Unit Cell
#		defined by the integer vectors in N
#
#---------------------------------------------------------------------------#

ucs_in_UC(N::Real; kwargs...) = ucs_in_UC(hcat(N); kwargs)

ucs_in_UC(N::Utils.List; kw...) = ucs_in_UC(LA.Diagonal(vcat(N...)); kw...) 

function ucs_in_UC(N::AbstractMatrix{T}; dim=2, method2D_cells=:Polygon, kwargs...) where T

	if T<:Int 

		# Nothing to do, N is ready.

	elseif T<:Float64 || any(n->isa(n,Float64), N)

		tN = trunc.(N) 

		isapprox(tN, N, atol=TOLERANCE) || error("N should contain integers")
		
		return ucs_in_UC(convert(AbstractMatrix{Int}, tN); kwargs...)

	else 

		error("Type '$T' not supported")

	end 


	correct_nr = Int(round(abs(LA.det(N))))
	# this should be the size of the new unit cell -- like | a1.(a2 x a3) |


	polygon = BodyVertices_fromVectors(N; dim=dim) 

	mM = cat(cat.(extrema(polygon, dims=dim)..., dims=[2,1][dim])..., dims=dim)


#  if len(N) == 1:
#
#    iter_ucs = [np.arange(*maxmin[::-1]).reshape(-1,1)]
#

	n, iter_ucs = if size(N,2)==2 && method2D_cells == :Polygon
									
									ucsUC_Polygon(N, mM, polygon; dim=dim, asfunction=true)

								else 

									error("Not implemented")

								end 

#    elif method2D_cells == "Spiral":
#      iter_ucs = ucsUC_Spiral(N,maxmin,correct_nr)
#
#    else:
#      raise ValueError("The second argument, 'method2D_cells', must be either 'Polygon' or 'Spiral'.")
#
#
#    
#  else:
#
#    iter_ucs = method2(N,maxmin)
#
#
	for i=1:n 

		ucs = iter_ucs(i)

		size(ucs, dim)==correct_nr && return ucs
		
		# if the correct number of uc-s in the large UC is found

	end 

	error("Could not construct big unit cell correctly.")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

#function CollectAtoms(lattices::Utils.List, latt_vac::Symbol, label)::AbstractMatrix{Float64}
#
#	hcat(map(lattices) do L 
#
#		L::Lattice
#
#		!hasproperty(L, latt_vac) && error("'$latt_vac' not found")
#
#		return get(getproperty(L, latt_vac), k, zeros(Float64, VecDim(L), 0))
#
#	end...) 
#
#end 

function parse_input_RsNs( latt::Union{Nothing,Lattice}=nothing;
													 Rs::Union{Nothing,AbstractVecOrMat{<:Real}}=nothing,
													 Ns=nothing,
													 A::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
													 StopStart=nothing,
													 dim=2,
													 kwargs...
													)

	!isnothing(Rs) && return Rs 

	A = Utils.Assign_Value(A, LattVec(latt))

	isnothing(A) && error("Cannot proceed with unknown vectors")


#nr_vect = size(A,dim)
#
#dim==2: size(A,dim)==size(N,[2,1][dim])
#
#dim==1: size(N,dim)==size(A,[2,1][dim])


#size([N,A][dim],dim) = size([A,N][dim],[2,1][dim])


#@show Ns 
#Correct_LattCoeff(latt, Ns)














	!isnothing(Ns) && return CombsOfVecs(A, 
																			 if ndims(Ns)==2 Ns else
[fill,reshape][ndims(Ns)+1](Ns, ([1,size(A,dim)][dim], [size(A,dim),1][dim]))
																			end,
																			 dim=dim)

	!isnothing(StopStart) && return CombsOfVecs10(A, StopStart...; dim=dim)

	error("Not enough non-nothing kwargs")

end

function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ucs::AbstractMatrix{<:Real};
											 kwargs...)::AbstractMatrix{Float64}

	Algebra.FlatOuterSum(atoms, ucs; dim=get(kwargs,:dim,2))


end 


function Atoms_ManyUCs(atoms::AbstractMatrix{<:Real},
											 ns::AbstractMatrix{<:Real},
											 A::AbstractMatrix{<:Real};
											 kwargs...)::AbstractMatrix{Float64}

	Atoms_ManyUCs(atoms, parse_input_RsNs(;Ns=ns, A=A, dim=2, kwargs...))

end 



function Atoms_ManyUCs(latt::Lattice;
#											 UCs::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
#											 Ns::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
#											 A::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
#											 StopStart=nothing,
											 kwargs...)::AbstractMatrix{Float64}

	Atoms_ManyUCs(PosAtoms(latt; dim=2, kwargs...),
								parse_input_RsNs(latt; dim=2, kwargs...)
								)
	

end







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Superlattice(latt::Lattice, n; kwargs...)
	
	Superlattice([latt], [n]; kwargs...)

end


function Superlattice(Components::Utils.List, Ns::Utils.List; Labels=nothing, kwargs...)

	Ns = Correct_LattCoeff.(zip(Components,Ns))

	Supervectors = CombsOfVecs(Components[1], Ns[1])

	for S in CombsOfVecs.(zip(Components[2:end], Ns[2:end]))

		isapprox(Supervectors, S, atol=TOLERANCE) && continue

		error("The 'Ns' provided are not compatible; one should recover the same supervectors for each pair (N,Component).")
	
	end



	Ns = ucs_in_UC.(Ns; kwargs...)

	Labels = 	if isnothing(Labels)
	
							repeat([""], length(Components))

						else 

							Utils.Assign_Value.(Labels, "")

						end 

	empty_labels = isempty.(Labels)


	return Lattice(Supervectors, map([:Sublattices, :Vacancies]) do p

		Utils.flatmap(unique(vcat(sublatt_labels.(Components; kind=p)...))) do k 

			all_atoms = map(zip(Components,Ns)) do (latt,ns)

				Atoms_ManyUCs(latt; Ns=ns, A=Supervectors, label=k, kind=p)
			
			end 

			good_atoms = .!isempty.(all_atoms)

			together_labels = good_atoms .&   empty_labels
			separate_labels = good_atoms .& .!empty_labels
			

			return vcat(if !any(together_labels) 
									
										[] 
									
									else 

										 k=>hcat(all_atoms[together_labels]...)

									end,

									map(findall(separate_labels)) do i 

										string(Labels[i], k) => all_atoms[i]

									end)

		end  # Utils.flatmap ends 

	end...) # map over [:Sublattices, :Vacancies] ends

end




#===========================================================================#
#
# 			some safety measureas when adding atoms 
#
#
#---------------------------------------------------------------------------#

function Initialize_PosAtoms(latt::Lattice, pos_atoms...)

	Initialize_PosAtoms(VecDim(latt), pos_atoms...)

end 

function Initialize_PosAtoms(vect_dim::Int, pos_atoms=nothing)
	
#    if type(pos_atoms)==type(None) or np.size(pos_atoms)==0:
#      print("\n *** No position given. The atom is set at the origin. ***\n")

	# if nothing is given, take the origin

	pos_atoms = to_myMatrix(if isnothing(pos_atoms) || isempty(pos_atoms) 

														zeros(Float64, vect_dim, 1)

													else 

														pos_atoms

													end)

	dim,nr = size(pos_atoms)

	if dim>vect_dim

		error("Atoms cannot have more coordinates than the dimension of the lattice vectors!")

	elseif dim==vect_dim

		return pos_atoms

	else 

	# in case additional coordinates are needed, use zeross

		return vcat(pos_atoms, zeros(Float64, vect_dim-dim, nr))

	end 

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function parse_input_sublattices(latt::Lattice, inp::Nothing, args...)

	good = !in(sublatt_labels(latt; kind=:Both))

	for c in string.('A':'Z')

		good(c) && return parse_input_sublattices(latt, [c], args...)

	end 

	error()

end 

function parse_input_sublattices(latt::Lattice, inp::T, args...) where T<:Union{AbstractString, Char, Symbol, Real}

	parse_input_sublattices(latt, [inp], args...)

end 


function parse_input_sublattices(latt::Lattice, inp::Utils.List, nr_at, method)

	repeat(map(method, inp), outer=max(1,Int(ceil(nr_at/length(inp)))))[1:nr_at]
		# repeat the list if necessary

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function AddAtoms!(latt::Lattice, positions=[], labels="A"; kind=:Sublattices)

	positions = Initialize_PosAtoms(latt, positions)

	new_labels = parse_input_sublattices(latt, labels, size(positions,2), string)

	D = getproperty(latt, kind)


	for (k,i) in zip(Utils.Unique(new_labels, inds=:all)...)

		isempty(i) && continue

		D[k] = hcat(get(D, k, EmptyPos(latt)), view(positions, :, i))

	end 

	return latt 

end 







#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function ShiftAtoms!(latt::Lattice; n=nothing, r=nothing, kind=:Both)

	R = parse_input_RsNs(latt; Ns=n, Rs=r, dim=2)

	for K in (kind==:Both ? [:Sublattices, :Vacancies] : [kind])

		for l in sublatt_labels(latt, kind=K)
			
			getproperty(latt, K)[l] .+= R 

		end 

	end 

	return latt 

end


function Shift_Atoms(latt::Lattice; n=nothing, r=nothing, kind=:Both)

	R = parse_input_RsNs(latt; Ns=n, Rs=r, dim=2)

	return Lattice(

					LattVec(latt),

					map([:Sublattices, :Vacancies]) do K 
							
						D = getproperty(latt, K)

						kind in [:Both,K] || return copy(D)

						return [l=>D[l].+R for l in sublatt_labels(latt, kind=K)]

					end...)

end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Reduce_Dim(latt::Lattice; dim=nothing, complement=false)

	if complement && (isnothing(dim) || length(dim)>1)
		
		error("The complement can be computed only if one single dimension is given.")

	end 


	keep_dims = if complement 

					[dim[1]]

							else 

					isnothing(dim) ? [] : filter(!in(vcat(dim...)),1:LattDim(latt))

						end

	return Lattice(LattVec(latt)[:,keep_dims],
								 copy(latt.Sublattices),
								 copy(latt.Vacancies)
								 )

end 

   
function Reduce_Dim!(latt::Lattice; dim=nothing, complement=false)

	if complement && (isnothing(dim) || length(dim)>1)
		
		error("The complement can be computed only if one single dimension is given.")

	end 


	keep_dims = if complement 

					[dim[1]]

							else 

					isnothing(dim) ? [] : filter(!in(vcat(dim...)),1:LattDim(latt))

						end

	latt.LattVec = latt.LattVec[:,keep_dims]

	return latt 

end 










#############################################################################
end
