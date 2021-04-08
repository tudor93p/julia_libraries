module Utils
#############################################################################



import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
import Dates, Combinatorics


import DelimitedFiles; const DlmF = DelimitedFiles
import FileIO#,JLD

using OrderedCollections:OrderedDict



const List = Union{AbstractVector, AbstractSet, Tuple}

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function DistributeBallsToBoxes(balls::Int, boxes::Int)

	map(Combinatorics.combinations(1:(balls+boxes-1), boxes-1)) do d

		diff(vcat(0, d, balls+boxes)) .- 1

	end

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function ReplaceByNeighbor!(bad_element::Function, A::AbstractArray, start=CartesianIndex(first.(axes(A))))
	
	cart_ind = findnext(bad_element, A, CartesianIndex(start)) 

	isnothing(cart_ind) && return A


	for k in 1:maximum([sum(abs, cart_ind.I .- C) for C in [1, size(A)]])
					# neighbors 

		for J in sort(DistributeBallsToBoxes(k, ndims(A)), by=LA.norm)

			for S in Base.product(fill([1,-1],ndims(A))...)	# Signs

				new_cart_ind = CartesianIndex(cart_ind.I .+ Tuple(S.*J))

				all(1 .<= new_cart_ind.I .<= size(A)) || continue 

				a = A[new_cart_ind]

				bad_element(a) && continue

				A[cart_ind] = a 

				return ReplaceByNeighbor!(bad_element, A, cart_ind)

			end 

		end 

	end 
				
	error("The whole array has 'unwanted' elements")

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function multiply_elementwise(A, B, dim=nothing)


	any(x -> isa(x,Number), (A,B)) && return A*B

	sA,sB = filter.(!isequal(1), size.((A,B)))

	inds = filter(!isnothing,  indexin(vcat(sA...), vcat(sB...)))

	if length(inds)!=minimum(length,(sA,sB)) || any(diff(inds).<=0)

		error("Shape mismatch")
	
	end 

	ndims(A)==ndims(B) && return A.*B


	if issubset(sA, sB) # A is smaller than B => expand A

		new_size_A,inds = vcat(size(B)...), vcat(1:ndims(B)...)

		for i in inds 

			!isnothing(dim) && i==dim && continue

			nr = count(new_size_A[i].==sA)

			nr==0 && setindex!(new_size_A, 1, i)

			nr>1 && error("Ambiguous result. Too many occurences found")

		end

		return reshape(A, new_size_A...).*B
		

	elseif issubset(sB, sA)
		
		return multiply_elementwise(B, A)

	end 

	error()

end 


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Rescale(A, mM0, mM1=A)

  m0,M0 = extrema(mM0)

	length(A)==1 && return A-A .+ m0

	m,M = extrema(mM1)
	
	return (A .- m)*(M0-m0)/(M-m) .+ m0

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#
is_float(S) = any(Ti -> S<:Ti, [Float64, Complex{<:Float64}])

is_exact(S) = any(Ti -> S<:Ti, [Integer, Rational, AbstractString, AbstractChar, Complex{<:Integer}, Complex{<:Rational}])



function Unique!(V::AbstractArray{T}; sorted=false, kwargs...) where T

	i = Unique(V; kwargs..., sorted=false, inds=:first)[2]

	deleteat!(V, filter(!in(i), axes(V,1)))

	sorted && sort!(V)

end



function Unique(V::AbstractArray{T}; 
								tol=1e-8, inds=nothing,
								sorted=false,
								check_type=true) where T

#	v isa AbstractArray &&
#	isList(v) || error("Type ",typeof(v)," not supported.")

	if !check_type || is_exact(T) || all(is_exact.(typeof.(V)))

		U = (sorted ? sort : identity)(unique(V))
	
		isnothing(inds) && return U

		i_f = findfirst(inds.==[:first,:all,:last])

		isnothing(i_f) && error("Cannot get indices with '$inds'")

		f = [findfirst, findall, findlast][i_f]

#		typeof(inds)<:Function || error("Cannot get indices with ", typeof(inds))
	
		return U, [f(isequal(u), V) for u in U]

	end 


	function get_newV()

		is_float(T) && return round.(V, digits= Int(ceil(-log10(tol))))

		fs, newT = typeof.(V) |> t -> (findall(is_float, t), promote_type(t...))

		W = AbstractArray{newT}(V)
	
		W[fs] = round.(W[fs], digits=Int(ceil(-log10(tol))))

		return W

	end

	I = Unique(get_newV(); tol=tol, inds=Assign_Value(inds, :first),
						 						 sorted=sorted, check_type=false)[2]

	isnothing(inds) && return V[I]

	return V[first.(I)], I


end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function flatmap(f,itr)

	vcat(map(f,itr)...)

end

function flatmapif(f, pred, itr)

	vcat(mapif(f, pred, itr)...)

end 

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function mapif(f::Function,pred::Bool,itr)

	return pred ? map(f,itr) : []

end

function mapif(f::Function,pred::Function,itr)

	filter(pred,map(f,itr))

end

function zipmap(args...)

	zip(map(args...)...)

end 


invmap(arg, fs...) = invmap(arg, fs)

function invmap(args, fs)



	map(f->f(args...), fs)

end 




#===========================================================================#
#
#	Dense matrix from lists of indices and values
#				(each row in inds represents a catesian index)
#
#---------------------------------------------------------------------------#

function Array_from_ListIndsVals(inds::AbstractMatrix, vals)

	inds, vals = collect(inds), collect(vals)

	isList(vals) || error("Wrong type")

	v0 = first(vals)

	A = zeros(typeof(first(v0)), size(v0)..., maximum.(eachcol(inds))...)

	for (I,V) in zip(eachrow(inds), vals)
		
		A[fill(:, ndims(V))..., I...] = V 

	end 

	return A

end

#===========================================================================#
#
#	Cumulative sum for a list of lists
#
#---------------------------------------------------------------------------#

function recursive_end(array::T) where T


	T <: Number && return array

	isempty(array) && return 0

	T <: AbstractArray{<:Number} && return array[end]

	i = findlast(!isempty,array)

	isnothing(i) && return 0

	return recursive_end(array[i])


end



function recursive_cumsum(array::T,start=0) where T

	T <: Number && return array + start

	nonempty_arrayinds = findall(!isempty,array)


	isempty(nonempty_arrayinds) && return array


	T <: AbstractArray{<:Number} && return cumsum(array) .+ start


	out = copy(array)

	for (order_nonempty,array_index) in enumerate(nonempty_arrayinds)

		out[array_index] = recursive_cumsum(array[array_index],

																				if order_nonempty==1 start else

									recursive_end(out[nonempty_arrayinds[order_nonempty-1]]) 

																				end)
	end

	return	out

end







#===========================================================================#
#
#	Apply function on a given axis. Interate the rest 
#
#---------------------------------------------------------------------------#



function ApplyF_OnAxis(f::Function, y::AbstractArray, dim::Int64=1)

	inds = setindex!(collect(Any, axes(y)), [:], dim)

	i1 = first.(inds)

	y_new1 = f(y[i1...])

	y_new = zeros(eltype(y_new1),
								setindex!(collect(size(y)), length(y_new1), dim)...)

	y_new[i1...] = y_new1

	for i in Base.Iterators.drop(Base.product(inds...),1)

		y_new[i...] = f(y[i...])
		
	end

	return y_new

end



#===========================================================================#
#
#	Identify consecutive integers in a list
#
#---------------------------------------------------------------------------#

#[0,0,0,0,0,0,3,1,1,1] => [1:6,7:7,8:10]

function IdentifySectors(list,sectors=[],start=0; tol=1e-8)

	for i in axes(list,1)

		if isapprox(list[i],list[1],atol=tol)		# within the sector
		
			i==length(list) && return vcat(sectors,[1+start:i+start])
																			# sector ended abruptly

		else  # sector ended at i-1

			return IdentifySectors(list[i:end],
														 vcat(sectors,[1+start:i-1+start]),
														 start+i-1;
														 tol=tol)

		end
	end

end

function IdentifyRanges(list)

	D2 = (diff(diff(list)).==0)

	!any(D2) && return list

	sectors =	map(filter(s->all(D2[s]),IdentifySectors(D2))) do s

							return minimum(s):maximum(s)+2
						end


	conflict(S) = (j -> maximum(S[j-1])>=minimum(S[j]))

	while true 

		i = findfirst(conflict(sectors),2:length(sectors))

		isnothing(i) && break

		s1,s2 = sectors[i:i+1]

		if length(s1) < length(s2)
			
			sectors[i] = minimum(s1):minimum(s2)-1

		else

			sectors[i+1] = maximum(s1)+1:maximum(s2)

		end

		sectors = filter(s->length(s)>2,sectors)

	end

#	if isempty(sectors)
#
#		d = filter(di->di>0,Algebra.FlatOuterDiff(q,q)[:])
#		
#		u = unique(d)
#		
#		c = zeros(Int64,length(u))
#		
#		for n in d c[indexin(n,u)[1]] +=1  end
#
#	end

	sector(i) = findfirst(s->i in s,sectors) 

	vcat(map(IdentifySectors(map(sector,axes(list,1)))) do S

		!in(S,sectors) && return list[S]

		(a,b) = extrema(S)
		
		step = list[a+1]-list[a]

		return [step==1 ? (list[a]:list[b]) : (list[a]:step:list[b])]

	end...)


end




#===========================================================================#
#
# check if arg is a vector or a tuple 
#
#---------------------------------------------------------------------------#




function isTuple(arg::T, Ti=Any) where T

	(T<:Type ? arg : T) <: Tuple{Vararg{<:Ti}}

end



function isList(arg::T, Ti=Any) where T

	isTuple(arg, Ti) && return true 

	for S in [AbstractVector, AbstractSet]

		(T<:Type ? arg : T) <: S{<:Ti} && return true 

	end 

	return false 

end

#===========================================================================#
#
#	Block diagonal matrix 
#
#---------------------------------------------------------------------------#

function BlkDiag(arg::T) where T
	
	isList(T,AbstractMatrix) && return cat(arg...,dims=(1,2))

	T<:Base.Generator && return BlkDiag(collect(arg))

	T<:AbstractMatrix{<:Number} && return arg

	T<:AbstractVector{<:Number} && return LA.diagm(0=>arg)

	T<:Number && return hcat(arg)

	isList(T,Any) && return BlkDiag(map(BlkDiag,arg))

	error("Input not understood. Type: ",T)

end

BlkDiag(args...) = BlkDiag(args)





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function flat(list_of_lists...; keep=nothing)

	i = findfirst(isList, list_of_lists)

	if isnothing(i) 
		
		isnothing(keep) && return vcat(list_of_lists...)
		
		return filter(keep, vcat(list_of_lists...))

	end

	return flat(list_of_lists[1:i-1]..., 
							list_of_lists[i]..., 
							list_of_lists[i+1:end]...;
							keep=keep)

end 

#===========================================================================#
#
# Given a relationship (item of type1) => (item of type2),
# Builds two dictionaries and returns two functions:
#				 - first gives the type-2 partner of a type1-item
#				 														(as the initial dict would do)
#				 - the second gives the type-1 partner of a type-2 item
#																		(as an inverse dict would do)
#
#
#---------------------------------------------------------------------------#


function FindPartners(pairs12;sortfirst=nothing)

	pairs12 = Dict(pairs12)

	isempty(pairs12) && return nothing,nothing

	pairs21 =	begin

							local K = vcat(keys(pairs12)...)
					
							isa(sortfirst,Bool) && sortfirst && sort!(K)

							isa(sortfirst,Function) && sort!(K, by=sortfirst)

							local V = unique(vcat(values(pairs12)...))

							Dict(v=>filter(k->in(v,vcat(pairs12[k])),K) for v in V)

						end



	return (key1 -> get(pairs12, key1, nothing), 
					key2 -> get(pairs21, key2, nothing))


#	return (key,label="") -> label==label1 ? pairs12[key] : pairs21[key]
end


#===========================================================================#
#
#	Quadrant of angle
#
#---------------------------------------------------------------------------#


function Quadrant(A::AbstractMatrix)

	size(A,2) == 2 && return Quadrant.(eachrow(A))

	size(A,1) == 2 && return Quadrant.(eachcol(A))

	error("Input not understood")

end

function Quadrant(xy::AbstractVector)

	length(xy) == 2 && return Quadrant(atan(xy[2],xy[1]))

	error("Input not understood")

end


function Quadrant(theta::Number)

	theta > pi && return Quadrant(theta-2pi)
	theta <-pi && return Quadrant(theta+2pi)

	for (i,(m,M)) in enumerate(([0,pi/2],[pi/2,pi],[-pi,-pi/2],[-pi/2,0]))
	
		m<=theta<=M && return i

	end


end

#===========================================================================#
#
# Unit matrix (built-in expression not intuitive)
#
#---------------------------------------------------------------------------#

function UnitMatrix(d::Int64,dtype=Float64)

  return Matrix{dtype}(LA.I,d,d)

end

function UnitMatrix(A::AbstractMatrix)

	one(A)

end

function UnitMatrix(n::Nothing=nothing)

	nothing

end

#===========================================================================#
#
# put "@time" if conditon==true
#
#---------------------------------------------------------------------------#

#macro timeif(cond,funwithargs)
#
#  @show funwithargs
#
#  cond && return @time eval(funwithargs)
#
#
#
#  return eval(funwithargs)
#
#end
#
#macro printif(cond,args...)
#
#  cond && println(join(args," "))
#
#end

#function timeif(cond)
#
#  cond==false && return f(expr,args...) = eval(expr)
#
#  return function f(expr,args...)
#
#    println(join(args," "))
#  
#    return @time eval(expr)
#
#  end
#
#
#end






#===========================================================================#
#
# Make a big sparse array with the elements of given array at certain indices
#
#---------------------------------------------------------------------------#

function Bigger_Matrix(new_Is,new_Js,size_Is,size_Js,matrix)

	# By converting matrix to sparse
    I,J,V = SpA.findnz(SpA.sparse(matrix))

	# Withouth converting to sparse -- much slower!
#    IJ = findall(!iszero,matrix)

#    (I,J),V = vcat.(Tuple.(IJ)...), getindex.([matrix],IJ)

    return SpA.sparse(new_Is[I],new_Js[J],V,size_Is,size_Js)

end

#===========================================================================#
#
# convert from string to Symbol
#
#---------------------------------------------------------------------------#

function DictKey_Symbol(d)

  return Dict([Symbol(k)=>v for (k,v) in pairs(d)])

end


#===========================================================================#
#
# Distribute list of parameters to several processes and print progress
#
#---------------------------------------------------------------------------#

function Distribute_Work(allparamcombs,do_work;arg_pos=1,kwargs0...)


  nr_scripts, idproc = get_arg.(1,arg_pos .+ (0:1) ,Int64)

  nr_scripts = min(length(allparamcombs),nr_scripts)

  if idproc > nr_scripts
    println(string("\nI am ",idproc," and I am not doing any jobs (nr.jobs < nr.processes).\n"))
    return

  end
  
  doparams, which_, njobs = distribute_list(allparamcombs,nr_scripts,idproc)
  
  println(string("\nI am ",idproc,"/",nr_scripts," and I am doing jobs ",which_," out of ",njobs,".\n"))


  function print_progress(ip,t1)

		println(string("\nI am ",idproc,"/",nr_scripts," and I completed ",ip,"/",length(which_)," jobs (last one: ",which_[ip],"/",which_," in ", Int(round(Dates.value(Dates.now()-t1)/1000)),"s)."))

#    println(strout)

#    open("/home/pahomit/progress.txt","a") do fout
#
#      DlmF.writedlm(fout,[strout])
#
#    end

    return Dates.now()

  end



  for (ip,p) in enumerate(doparams)
  
    time1 = Dates.now()
  
    do_work(p;kwargs0...)

    print_progress(ip,time1)
  
  
  end

#      return enumerate(doparams),print_progress,args_for_print
    




#  return ([args_for_print(ip),p][2] for (ip,p) in enumerate(doparams))

	return 


end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Slice_LikeArraysInList(arrays, X, dims=1)

	i,f = fill(Colon(), dims-1), fill(Colon(), ndims(X)-dims)

	ms = cumsum([0;length.(arrays)]) |> B -> [B[j-1]+1:B[j] for j=2:length(B)]

	return [X[i..., m, f...] for m in ms]

#	isnothing(out_vector) && return inds

#	return [out_vector[i] for i in inds]

end

#===========================================================================#
#
# Generate a list of n distinct random items from a list
#
#---------------------------------------------------------------------------#

function Random_Items(list,n;distinct=true)

  distinct && n > size(list,1) && error("Cannot have so many distinct items")

  new_list = []

  while length(new_list) < n
  
    new_item = rand(axes(list,1))
 
    (distinct & new_item in new_list) || push!(new_list,new_item)
  
  end

  return list[new_list]

end

#===========================================================================#
#
# Write results in different files
#
#---------------------------------------------------------------------------#

function Write_NamesVals(filename,  storemethod, names, result, bounds; tol=1e-7, filemethod="new")


	Write!, outdict = Write_NamesVals(filename, storemethod;
																		filemethod=filemethod,
																		tol=tol)


  for (i,name) in enumerate(names)

 		Write!(name, result[:,bounds[i]+1:bounds[i+1]], outdict)

  end

  return Write!, outdict

#  outdict = new_item("kLabels",reshape(repeat(kLabels,inner=div(size(result,1),length(kLabels))),:,1),outdict)

#  outdict = new_item("kPoints",kPoints,outdict)

#  return outdict


end

#===========================================================================#
#
# Write 
#
#---------------------------------------------------------------------------#

function has_significant_imagpart(matrix; atol=1e-8, rtol=1e-5, mute=false)

  re_matrix = real(matrix)

	inds =  findall(abs.(re_matrix) .> atol)

	if !isempty(inds)

		rel_im_matrix = imag(Array(matrix)[inds])./Array(re_matrix)[inds]

		inds2 = findall(abs.(rel_im_matrix) .> rtol)
	
		if !isempty(inds2) 
		
			

			!mute && println("Number of elements with significant imaginary part: ",
																length(matrix[inds][inds2])
											)

			!mute && println("Maximum relative imaginary part: ",
																maximum(abs.(rel_im_matrix[inds2]))
											)

			return true 
		end
	end

	return false

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Extension_Storemethod(storemethod)

	storemethod in ["dat","jld"] && return ".$storemethod"

	error("The 'storemethod' $storemethod is not supported")

end


function isLegendFile(storemethod)

	storemethod == "dat" && return x->occursin("_Legend",x)

	storemethod == "jld" && return x->false

	error("The 'storemethod' $storemethod is not supported")

end 

function Name_fromLegend(storemethod)

	storemethod == "dat" && return x->split(x,"_Legend")[1]
	
	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 

function LegendFile_fromName(storemethod)

	storemethod == "dat" && return x->x*"_Legend"

	storemethod == "jld" && return identity

	error("The 'storemethod' $storemethod is not supported")

end 




function Write_NamesVals(filename, storemethod="jld"; 
												 tol=1e-7, filemethod="new")


	function writable(matrix::AbstractArray{<:Number}, name)

		if has_significant_imagpart(matrix; atol=tol)
			
			println("Observable '$name'")

			!isnothing(filename) && println(filename(name))

			println()

			@warn "Observable '$name' not real!"

		end

		return real(matrix)

	end

	writable(matrix::AbstractArray{<:AbstractString},name) = matrix

	writable(matrix::AbstractArray{<:Char},name) = matrix

#	writable(matrix::Char, name) = matrix

	writable(x::Number, name) = writable(vcat(x), name) 

	function writable(D::AbstractDict, name) 
		
		Dict(k=>writable(v,name) for (k,v) in pairs(D))

	end 



	ext = Extension_Storemethod(storemethod)


  function write_(name, data)

    if !isnothing(filename) 

			fn = filename(name)*ext

			if storemethod=="dat"

	      open(fn, filemethod=="append" ? "a" : "w") do fout

  	       DlmF.writedlm(fout, data)

    	  end

			elseif storemethod=="jld"

#				JLD.save(FileIO.File(FileIO.format"JLD",fn), name, data)
				FileIO.save(FileIO.File(FileIO.format"JLD",fn), name, data)

			end 

    end
    
    return data 
		
  end



	function new_item!(name::AbstractString, val, outdict::AbstractDict=Dict())

#		@show name typeof(val) keys(outdict)

		if !isnothing(name) & !isnothing(val)

			outdict[name] = write_(name, writable(val,name))

		end 

    return outdict

  end

	function new_item!((name, val)::T, outdict::AbstractDict=Dict()) where T<:Tuple{AbstractString,Any}

		new_item!(name, val, outdict)

	end 	

  return new_item!,Dict()

end




function Delete_NamesVals(filename, names, storemethod)

	ext = Extension_Storemethod(storemethod)

	get_legend = LegendFile_fromName(storemethod)

	Ns = isa(names,AbstractString) ? [names] : names

#	println("\nDelete ",join(Ns,", "),"\n")

	for fn in filename.(Ns)

		rm.(filter(isfile, unique([fn,get_legend(fn)]).*ext))

	end 

end






function Read_NamesVals(filename, names, storemethod)
	
	names = isa(names,AbstractString) ? [names] : names

	filenames = unique(filename.(names).*Extension_Storemethod(storemethod))


	if storemethod=="jld"

		FNs = filter(isfile, filenames)

		isempty(FNs) && return Dict()

		return merge(map(FNs) do fn 

				FileIO.load(FileIO.File(FileIO.format"JLD",fn))

		end...) 

	elseif storemethod=="dat"

	  outdict = Dict()
	
		for (n,fn) in zip(names,filenames)

			if isfile(fn)

				outdict[n] =DlmF.readdlm(fn)

			end 
	
	  end
	
	  return outdict

	end 

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function is_dict_or_JLDAW(D)

	D isa AbstractDict && return true 

	all(in(propertynames(D)), [:keys, :values]) && return true

	return false 

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function FoundFiles_NamesVals(filename, names, storemethod)
	
	fn = filename.(vcat(names)).*Extension_Storemethod(storemethod)

	return all(isfile, unique(fn))

end




#===========================================================================#
#
# Prepare dict for writing by flattening its values 
# 										(if array, use cartesian inds for the new keys)
#
#---------------------------------------------------------------------------#

function flattenDictEntries(D::AbstractDict; onlykeys=false, 
																							onlyvals=nothing,
																							)

	any(occursin.("_",keys(D))) && error("The keys already contain the character '_'. The result will not be read correctly.")

	T = typeof( first(D).second )


	onlykeys && !isnothing(onlyvals) && error("The kwargs 'onlykeys' and 'onlyvals' are mutually exclusive")

#	getkeys(D) = [D[k] for k in Keys]

	if T<:Number 
		
		onlykeys && return collect(keys(D))
		
		isnothing(onlyvals) && return D

		onlyvals isa Bool && return collect(values(D))

		onlyvals isa AbstractString && return D[onlyvals]

		isList(onlyvals, AbstractString) && return [D[k] for k in onlyvals]

		error("'onlyvals=$onlyvals' not supported")

	end 
	
	T<:AbstractArray || error("$T not supported")

	key(k,I) = string(k, "_", join(Tuple(I),"-")) 

	iter = ((k,I,V[I]) for (k,V) in pairs(D) for I in CartesianIndices(V))


	onlykeys && return [key(k,I) for (k,I,v) in iter]

	onlyvals isa Bool && return [v for (k,I,v) in iter]

	if onlyvals isa AbstractString 
	
		for (k,I,v) in iter 

			key(k,I)==onlyvals && return v

		end

	end



	out = Dict(key(k,I)=>v for (k,I,v) in iter)

	isnothing(onlyvals) && return out

	isList(onlyvals, AbstractString) && return [out[k] for k in onlyvals]


	error("'onlyvals=$onlyvals' not supported")


end								



#function restoreDictEntries(d::AbstractDict)
#
#	!occursin("-",first(keys(d))) && return d 
#
#	K = collect(keys(d))
#
#	sK = hcat(split.(K,"-")...)
#
#	return Dict(map(unique(sK[1,:])) do label 
#
#		select = label.==sK[1,:]
#		
#		return label => Array_from_ListIndsVals(
#												transpose(parse.(Int, sK[2:end, select])),
#												[d[k] for k in K[select]])
#	end)
#
#end 

function restoreDictEntries(matrix::AbstractMatrix, legend)
											
	K = string.(legend[:])

	!all(occursin.("_",K)) && return Dict(zip(K,eachcol(matrix)))

	sK = hcat(map(legend) do L
					
						left,right = split(L,"_")

						return vcat(left,split(right,"-"))

					 end...)

#	column: [label, inds...]


	return Dict(map(unique(sK[1,:])) do label 

		select = label.==sK[1,:] 	# all elements corresp to label
		
		return label => Array_from_ListIndsVals(
												transpose(parse.(Int, sK[2:end, select])),
												eachcol(matrix[:,select])
																						)
	end)

end 







#===========================================================================#
#
#	Write physical observable. Possibly with a legend file
#
#---------------------------------------------------------------------------#


function Write_PhysObs(filename, storemethod; tol=1e-7)

	function assert_type(v)

		isa(v, Number) && return 1 
		
		isa(v, AbstractDict) && return 3 

		isList(v, Number) && return 2 


		error(typeof(v)," not supported")

	end 


	concat(f,vals) = vcat([reshape(f(v),1,:) for v in vals]...)
# above: preallocate and avoid vcat ?

	function get_legend(vals)

		assert_type(vals[1]) in [1,2] && return nothing
		
		return sort(flattenDictEntries(vals[1], onlykeys=true))

	end 

	function get_matrix(vals)

		assert_type(vals[1]) in [1,2] && return concat(vcat, vals)

		return concat(v->flattenDictEntries(v; onlyvals=get_legend(vals)), vals)


	end 

	function get_data(vals)

		assert_type(vals[1]) in [1,2] && return concat(vcat, vals)

		K = collect(keys(vals[1]))

		(T,S,I0) = vals[1][K[1]] |> a -> (eltype(a), size(a), fill(:, ndims(a)))

		data = Dict(k=>zeros(T, (length(vals), S...)) for k in K) 

		for (i,val) in enumerate(vals), (k,v) in pairs(val) 

			data[k][i, I0...] = v

		end

		return data

	end



	Write!, = Write_NamesVals(filename, storemethod; tol=tol)

	legend_file = LegendFile_fromName(storemethod)



	function new_item!((obs, vals), outdict::AbstractDict=Dict())

		new_item!(obs, vals, outdict)

	end



	function new_item!(obs, vals, outdict=Dict())
		
	
#		legend,matrix,data = legend_matrix_data(vals)
	
		if storemethod=="jld"

			#return 
			Write!(obs, get_data(vals), outdict)

		elseif storemethod=="dat"

			Write!(obs, get_matrix(vals)) 
								# ! careful ! this object is written, not returned 

			Write!(legend_file(obs), get_legend(vals), outdict)

			outdict[obs] = get_data(vals)

		end

		return outdict

	end
	
	
	return new_item!,Dict()

end




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

FoundFiles_PhysObs = FoundFiles_NamesVals


function Read_PhysObs(filename, names, storemethod)

	if storemethod=="jld"

		return Read_NamesVals(filename, names, storemethod)

	elseif storemethod=="dat"

		legf = LegendFile_fromName(storemethod)

		obsf = Name_fromLegend(storemethod)

		out = Read_NamesVals(filename, [names; legf.(vcat(names))], storemethod)
													
		for legend in filter(isLegendFile(storemethod), keys(out))

			obs = obsf(legend)

			!haskey(out, obs) && error("Legend exists, but no obs?")

			out[obs] = restoreDictEntries(pop!(out, obs), pop!(out, legend))
	
		end
	
		return out 

	end 

end 





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





function ChangeStoreMethod_PhysObs(filename, names, source, dest)


	if source!=dest 

		Write!, = Utils.Write_NamesVals(filename,  dest)

#		@show names vcat(names) source isLegendFile(source).(vcat(names))


		for name in filter(!isLegendFile(source), vcat(names))

#			println("\nMoving '$name' from '$source' to '$dest'\n")

#			@show 1,name 
#			M = Read_PhysObs(filename, name, source)[name]
#			@show 2,name
#			M isa AbstractDict && @show keys(M) size.(values(M))
#			M isa AbstractArray && @show size(M)

#			@show 3,name 
#			@show typeof(M)

#			Write!(name, M)
			Write!(name, Read_PhysObs(filename, name, source)[name])

			Delete_NamesVals(filename, name, source)

		end 

	end

	return FoundFiles_PhysObs(filename, names, dest)

#	return true # if succesful, files exist 

end 



#===========================================================================#
#
# Make directory tree if it doesn't exist
#
#---------------------------------------------------------------------------#

# function exsits already! 'mkpath'


#function MakeDir(filepath)
#
##  a = split(strip(path,'/'),'/')
#  a = split(filepath,'/')
#
#  for i in 1:length(a)-1
#
#    f = join(a[1:i],"/")
#
#    !isdir(f) && mkdir(f)
#
#  end
#
#  return filepath
#
#end
#


#===========================================================================#
#
# Add item to dictionary if "key" does not exist, 
#		otherwise return existing value
#
#---------------------------------------------------------------------------#

#function Update_Dict(dict,key,f,args)
#
#  if !(key in keys(dict))
#  
#    dict[key] = f(args)
#
#  end
#
#  return dict[key]
#
#end


#===========================================================================#
#
# Read variables from the arguments of the command
#
#---------------------------------------------------------------------------#

function get_arg(default,i,typ=String)

  if length(ARGS) >= i 

    typ!=String && return parse(typ,ARGS[i])

    return ARGS[i]

  end

  return default

end





#===========================================================================#
#
# "inverse map" --> applies the same function to one argument 
#
#---------------------------------------------------------------------------#


function imap(fs,args...)




end 



#===========================================================================#
#
# Executes a julia function in parallel
#
#---------------------------------------------------------------------------#




function execute_julia_function_(file_function;args=[],Nr_Procs=8,libs=[])


  lib_local  = "/media/tudor/Tudor/Work/scripts/julia_libraries/"
  lib_remote = "/net/horon/scratch/pahomit/apps/julia-1.2.0/lib/mylib/"

  lib_root = gethostname() == "tudor-HP" ? lib_local : lib_remote

		# make up some name for the script to be launched
  frun = string("auxrun_",join(file_function),rand(40000:50000),".jl")


  open(frun,"w") do f
    filename,funname = file_function

    text = string("include(\"",filename,".jl\")\n",funname,"()")
		# inclunde the file and run the target function

    write(f,text)
 #   println(text)
  end

  libs = "-L " *lib_root.*libs.*".jl"

  cmd = [["julia -p",Nr_Procs];libs;[frun];args]

  cmd = join(string.(cmd),' ')
		# give the input and output files and run in parallel

#  println(cmd)

  run(`sh -c $cmd`)

  rm(frun)

end



function execute_julia_function(args,nr_results,read_write_funs,file_function;Nr_Procs=8,libs=[])

  !isdir("./aux") && mkdir("./aux")

  write_args, read_result = read_write_funs



		# write the given arguments with the particular method
  f_args = write_args(args...)


		# make up some filenames for the results
  f_result = "./aux/res_"*join(file_function).*string.(rand(20000:30000) .+ collect(1:nr_results) )

		# make a separate script and run the function in parallel 
  execute_julia_function_(file_function,args=[f_args;f_result],Nr_Procs=Nr_Procs,libs=libs)
		# read the results with the particular method
  return read_result(f_result)

end



#===========================================================================#
#
# Distributes a list to n ~equal parts
#
#---------------------------------------------------------------------------#

function distribute_list(list,n,i)

  nitems = size(list,1)

  inds  = cumsum([0;div(nitems,n) .+ (1:n.<=nitems%n)])


  which_ = inds[i]+1:inds[i+1]


  return list[which_], which_, nitems


end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function DictRandVals(params::T) where T<:AbstractDict

	T(k=>(typeof(v) <: AbstractVector) ? rand(v) : v 
												for (k,v) in pairs(params))

end

function DictRandVals(params::T) where T
	
	isList(T,AbstractDict) && return DictRandVals.(params)

	error("Type not understood: $T")

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#



function DictFirstVals(params::T) where T<:AbstractDict

	T(k=>(typeof(v) <: AbstractVector) ? v[1] : v 
							for (k,v) in pairs(params))

end

function DictFirstVals(params::T) where T
	
	isList(T,AbstractDict) && return DictFirstVals.(params)

	error("Type not understood: $T")

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function AllValCombs(params::T;
										 constraint=nothing,sortby=nothing) where T<:AbstractDict

	vals(v) = (typeof(v) <: AbstractVector) ? v : [v]

	K = vcat(keys(params)...)

	dicts = [T(zip(K,vs)) 
					 for vs in Base.product([vals(params[k]) for k in K]...)][:]

  !isnothing(constraint) && filter!(constraint,dicts)

  !isnothing(sortby) && sort!(dicts,by=sortby)

  return dicts
end


function AllValCombs(params::AbstractVector{<:OrderedDict};kwargs...)

	collect(Base.product(AllValCombs.(params;kwargs...)...))[:]

end

#===========================================================================#
#
# Construct a NamedTuple from a list of keys and a dict/NT/array
#
#---------------------------------------------------------------------------#

function NT(data,Keys=nothing)


  if isa(data,AbstractDict) | isa(data,NamedTuple)

		Keys = Assign_Value(Keys,vcat(keys(data)...))

    data = [data[k] for k in Keys]
  end

#  data = try [data[k] for k in Keys]
#
#	 catch error
#
#           !any(isa.(error,[MethodError,ArgumentError])) && error("Please provide a list-type container (Array/...) or an associative collection (Dict/NamedTuple) with the right keys.")
#
#	   data
#
#	 end

  return NamedTuple{Tuple(Keys)}(data)

end


#===========================================================================#
#
# NamedTuple to Arrays (python compatibility)
#
#---------------------------------------------------------------------------#

function NT_toLists(NT)

  return [string(k) for k in keys(NT)], [NT[k] for k in keys(NT)]


end


#===========================================================================#
#
# For the common keys, the non-nothing values in 'nt' replace those in NT
#
#---------------------------------------------------------------------------#

function NT_ReplaceItems(NT,nt=NamedTuple{}())

  vals = Dict([k=>(k in keys(nt) ? nt[k] : NT[k]) for k in keys(NT)])

  K = Tuple([k for k in keys(NT) if !any(isnothing.(vals[k]))])

  return NamedTuple{K}([vals[k] for k in K])

end




function NT_KeepItems(NT,keep_keys=[])

  K = Tuple(intersect(keys(NT),keep_keys))

  return NamedTuple{K}([NT[k] for k in K])

end



function NT_AllCombs(params::NamedTuple;constraint=nothing,sortby=nothing)

  NTs = map(Base.product(values(params)...)) do vs

    return NamedTuple{keys(params)}(vs)
  end[:] 


  !isnothing(constraint) && filter!(constraint,NTs)

  !isnothing(sortby) && sort!(NTs,by=sortby)


  return NTs

#  vals = map(Base.product([Base.product(p...) for p in params]...)) do V
#
# 
#           [NamedTuple{ks}(vs) for (vs,ks) in zip(V,keys.(params))]
#
#         end
#

#  if isa(iter_last,Symbol)
#
#    iter_last = [iter_last]
#
#  elseif isa(iter_last,AbstractArray)
#
    length(iter_last) == 0 && return NT_AllCombs_(params...)
# 
#  else
#
#    error("'iter_last' must be 'Symbol' or Array of 'Symbol'")
# 
#  end
#
#
#  length(intersect(keys.(params)...,iter_last)) > 0 && error()
#
#
#  aux = NamedTuple{Tuple(iter_last)}([[nothing] for k in iter_last])
#
#  iter1 = [NT_ReplaceItems(p,aux) for p in params]
#		# replace desired item with nothing
#
#  iter2 = NamedTuple{Tuple(iter_last)}(
#		map(k->[p[k] for p in params if k in keys(p)][1],iter_last))
#
#
#  return [[[NT_ReplaceItems(Pi,p) for Pi in P]
#		for (p,) in NT_AllCombs(iter2)]
#			for P in NT_AllCombs(iter1...) ]
#

end



#function NT_AllCombs_(params...)
#
#
#  vals = map(Base.product([Base.product(p...) for p in params]...)) do V
#
# 
#           [NamedTuple{ks}(vs) for (vs,ks) in zip(V,keys.(params))]
#
#         end
#
#
# # vals = [[vcat(vi...) for vi in v] for v in vals]
##
#
##  inds = Base.product([Base.product(axes.(p,1)...) for p in params]...)
#
##  inds = [[vcat(ii...) for ii in i] for i in inds]
#
##  println.(inds)
#
#
##  if !isnothing(names) & !isnothing(files_exist)
#
#
##
##    names = .Assign_Value(names,string.(["Parameter "],1:size(params,1)))
##   
##   
##  
## 
###  state(i) = files_exist(value(i))
##
##
##
##
##
##
## 
##  function nr_occur(u,c)
##  
##    check = files_exist.( values( inds[c .== u,:] ))
##   
##    s,n = sum(check),length(check)
##  
##    return string(s,"/",n," (",round(100*s/n,digits=1),"%)")
##  
##  end
##  
##
##  
##  for (name,param,column) in zip(names,params,eachcol(inds))
##  
##    for u in unique(column)
##  
##      println("'",name,"' = ",round(param[u],digits=3), ": ",nr_occur(u,column))
##  
##    end
##  
##    println()
##  
##  end
##
##  end
##
##  return values(inds[.!files_exist.(values(inds)),:])
#
#  return vals[:]
#end 
#
#
#


#===========================================================================#
#
# Combine two lists of parameters (named tuples, dicts)
#
#---------------------------------------------------------------------------#

function Combine_NamedTuples(input_param,default_param)

  issubset(keys(input_param),keys(default_param)) || error("The parameters must be among "*join(map(string,keys(default_param)),", ")*".")

  return merge(default_param,input_param)

end




#===========================================================================#
#
# Pad a float left and right
#
#---------------------------------------------------------------------------#

function lrpad(x,N=[2,2];p='.')

	isa(N,Number) && return lrpad(x,[N,N];p=p)

  return rstrip(join([f(s[1:min(length(s),n)],n,'0') for (s,f,n) in zip(split(string(Float64(x)),'.'),[lpad,rpad],N)],p),p)

end


function nr2string(nr::T,digits=2) where T


	if T<:Union{Tuple,AbstractArray}# && any(t->typeof(nr[1])<:t,types)

		return map(n->nr2string(n,digits),nr) |> join
	
	end

	Union{Number,AbstractString} |> t-> T<:t || error("Please provide $t")



	isempty(digits) && return string(nr)

	digits isa Number && return nr2string(nr,[digits,digits])

  if nr isa AbstractString

    nr_ = tryparse(Float64,nr)
 
    return isnothing(nr_) ? nr : nr2string(nr_,digits)
 
  end


  left,right = split(string(round(Float64(abs(nr)),digits=digits[2])),'.')

  left = repeat("m",nr<0)*lpad(left,digits[1],'0')

  digits[2]==0 && return left

	right = rpad(right[1:min(end,digits[2])],digits[2],'0')

  return join([left,right],'p')

end

#===========================================================================#
#
#  Returns a path of n points which connects the inputed points
#
#---------------------------------------------------------------------------#

function PathConnect(points,n;end_point=true,bounds=[0,1],fdist=identity)

  size(points,1) == 1 && return points,[0]

  dist = LA.normalize(diff(points,dims=1) |>eachrow .|>LA.norm .|>fdist,1)

  n -= end_point


  ns_ = max.(1,Int64.(round.(dist*n)))

  while sum(ns_) > n 
    ns_[argmax(ns_)] -= 1
  end

  while sum(ns_) < n
    ns_[argmin(ns_)] += 1
  end
  

  path = vcat(map(axes(points,1)[1:end-1]) do i

    ep = end_point && i==axes(points,1)[end-1]

    t = range(0, 1, length = ns_[i] + 1 )[1:end-1+ep]

    return hcat(1 .- t , t)*points[i:i+1,:]

  end...)



  xticks = bounds[1] .+ cumsum([0;dist]).*diff(bounds)
#  xticks = cumsum([1:ns_])
 
  return path,xticks

end



#===========================================================================#
#
# return F(args) = f1(args) + f2(args) + ...
#
#---------------------------------------------------------------------------#

function Sum_functions(functions...)::Function
 

#  functions = [f for f in functions if isa(f,Function)]

#  return (args...) -> mapreduce(F -> F(args...),+, Functions)
  return (args...) -> reduce(+,map(f -> f(args...), functions))

end


#===========================================================================#
#
# Make array out of list of lists, like np.array(...)
#
#---------------------------------------------------------------------------#

function subarray(A,js=[])

  (isa(A,Number) | (length(js)==0)) && return A

  return subarray(A[js[1],axes(A)[2:ndims(A)]...],js[2:length(js)])

end


function get_shape(A)

  shape = Int64[]

  while true
  
    isa(A,Number) && return Tuple(shape)

    push!(shape,size(A,1))

    A = subarray(A,[1])
   
  end

end

	# much faster than my function
#function nparray(L)::Array
#
#  np = PyCall.pyimport("numpy")
#
#  return np.array(L)
#
#end

#function ListToArray(L,dtype=Float64)
#
#  A = zeros(dtype,get_shape(L))
#
#  for i in CartesianIndices(A)
#
#    A[i] = subarray(L,Tuple(i))
#
#  end  
#
#  return A
#
#end

#===========================================================================#
#
# vectors of integers, each component is a range
#
#---------------------------------------------------------------------------#

function vectors_of_integers(D::Int64, stop, start=-stop; dim=1, sortby=nothing)

	dim2 = [2,1][dim]


	boundaries = zipmap([Assign_Value(start,-stop), stop]) do x
		
		isa(x,Int) && return  fill(x,D)
	
		Utils.isList(x,Int) || error("Type ",typeof(x)," not supported")

		length(x)>=D && return vcat(x...)[1:D]

		length(x)==1 && return fill(x[1], D)

		error("Wrong dimensions")
	
	end 


	Ls = [1;[f-i+1 for (i,f) in boundaries];1]


	out = similar(Array{Int64}, [prod(Ls), D][[dim,dim2]]...)

	
	for (i,(a,b)) in enumerate(boundaries)

		setindex!(out,
							repeat(a:b, inner=prod(Ls[i+2:end]), outer=prod(Ls[1:i])),
							[:,i][dim], [:,i][dim2]
							)


	end

	!isa(sortby, Function) && return out

	return sortslices(out, dims=dim, by=sortby)





end



#===========================================================================#
#
#
#---------------------------------------------------------------------------#


#function flat_indsvals(i::Int64,j::Int64,v::Matrix{Complex{Float64}},d0::Int64=1) 
#
##  d0==1 && return [i],[j],[v[1]]
#
#  small_is, small_js, vs = SpA.findnz(SpA.sparse(v))
#
#  return hcat((i-1)*d0 .+ small_is, (j-1)*d0 .+ small_js, vs)
#
#
#end 




#===========================================================================#
#
# Check if we have the same number or the same array
#
#---------------------------------------------------------------------------#


function fSame(atol)

  function same(a::Number,b::Number=0.0,atol=atol)::Bool

    return isapprox(a,b,atol=atol)

  end



  function same(a::AbstractArray,b::AbstractArray,atol=atol)::Bool

    return isapprox(a,b,atol=atol)

  end



  function same(a::AbstractArray,b::Number=0.0,atol=atol)::Bool

    return isapprox(LA.norm(a),b,atol=atol)
  
  end



  function same(a::Number,b::AbstractArray,atol=atol)::Bool

    return isapprox(LA.norm(b),a,atol=atol)
  
  end


  return same

end



#function Same(a,b,order::Real=6)::Bool
#
#  return LA.norm(a.-b) < 1.0/10.0^order
#
#end


#===========================================================================#
#
# Transforms a n-dimensional array intro a list of (n-1)-dimensional arrays
#
#---------------------------------------------------------------------------#


function ArrayToList(a)

  return [a[i,axes(a)[2:ndims(a)]...] for i in axes(a,1)]

end


#===========================================================================#
#
# Assigns value if variable is nothing
#
#---------------------------------------------------------------------------#


function Assign_Value(input_value, std_value)

  isnothing(input_value) && !isnothing(std_value) && return std_value
  
  return input_value

end


function Assign_Value(input_value, get_std_value, args...; kwargs...)

	!isnothing(input_value) && return input_value
	
	!isnothing(get_std_value) && return get_std_value(args...; kwargs...)

  return nothing

end



#############################################################################

end
