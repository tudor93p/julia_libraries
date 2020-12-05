module Utils

import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
import Dates

import DelimitedFiles; const DlmF = DelimitedFiles

using OrderedCollections:OrderedDict




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Rescale(A, mM0)

  m0,M0 = extrema(mM0)

	m,M = extrema(A)

	return (A .- m)*(M0-m0)/(M-m) .+ m0

end

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function flatmap(f,itr)

	vcat(map(f,itr)...)

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
#	Block diagonal matrix 
#
#---------------------------------------------------------------------------#

function isList(arg::T,Ti=Any) where T

	if T <: Type 
		
		return any(S -> arg<:S, [AbstractVector{<:Ti},Tuple{Vararg{<:Ti}}])

	else

		return isList(typeof(arg),Ti)

	end
end

#isList(a::Any,Ti=Any) = isList(typeof(a),Ti)

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

							typeof(sortfirst)<:Function && sort!(K,by=sortfirst)

							local V = unique(values(pairs12))

							Dict(v=>filter(k->pairs12[k]==v,K) for v in V)
						end

	return key1 -> pairs12[key1], key2->pairs21[key2]


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

    println(string("\nI am ",idproc,"/",nr_scripts," and I completed ",ip,"/",length(which_)," jobs (last one: ",which_[ip],"/",which_," in ", round(Dates.value(Dates.now()-t1)/1000),"s)."))

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

function Write_NamesVals(filename, filemethod, names, result, bounds; tol=1e-7)

  new_item,outdict = Write_NamesVals(filename, filemethod, tol=tol)

  for (i,name) in enumerate(names)

    outdict = new_item(name, result[:,bounds[i]+1:bounds[i+1]], outdict)

  end

  return new_item, outdict

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




function Write_NamesVals(filename=nothing, filemethod="new"; tol=1e-7)

	function writable(matrix::AbstractArray{<:Number}, name)

		if has_significant_imagpart(matrix; atol=tol)
			
			println("Observable '$name'")

			!isnothing(filename) && println(filename(name))

			println()

			error("Observable '$name' not real!")

		end

		return real(matrix)

	end


	writable(matrix::AbstractArray{<:String},name) = matrix

	writable(matrix::AbstractArray{<:Char},name) = matrix

	writable(x::Number, name) = writable(vcat(x), name) 

  function write_(name,matrix)

    if !isnothing(filename) 

      open(filename(name), filemethod=="append" ? "a" : "w") do fout

         DlmF.writedlm(fout, matrix)

      end
    end
    
    return matrix
  end




  function new_item(name,val,outdict=Dict())

		outdict[name] = write_(name,writable(val,name))

    return outdict

  end


  return new_item,Dict()

end




function Read_NamesVals(filename,names)

  outdict = Dict()

  for name in (isa(names,String) ? [names] : names)

    if isfile(filename(name))

      outdict[name] = DlmF.readdlm(filename(name))

    end

  end

  return outdict

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

	return T(k=>(typeof(v) <: AbstractVector) ? rand(v) : v 
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

	return T(k=>(typeof(v) <: AbstractVector) ? v[1] : v 
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
##    names = Utils.Assign_Value(names,string.(["Parameter "],1:size(params,1)))
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

	types = [Number,String]


	if ((T<:Tuple) | (T<:AbstractArray))# && any(t->typeof(nr[1])<:t,types)

		return map(n->nr2string(n,digits),nr) |> join
	
	end

	any(t->T<:t,types) || error("Please provide $types")


	isempty(digits) && return string(nr)

	digits isa Number && return nr2string(nr,[digits,digits])

  if nr isa String

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

function vectors_of_integers(dim::Int64,stop,start=-stop)


#  if dim == 0
#    return 
#  end

  start = vcat(start); stop=vcat(stop)


  if length(start) != dim
    start = repeat(start,dim)
  end

  if length(stop) != dim
    stop  = repeat(stop,dim)
  end



  iter = Iterators.product([a:b for (a,b) in zip(start,stop)]...)

  return collect(transpose(hcat([collect(x) for x in vcat(iter...)]...)))




#  function f(stop_::Array{Int64,1},start_::Array{Int64,1})
#
#    if length(stop_)==1 & length(start_)==1
#      return collect(start_[1]:stop_[1])
#    end
#
#
#    current_level = f(stop_[1:1],start_[1:1])
#    lower_level = f(stop_[2:length(stop_)],start_[2:length(start_)])
#
#    inner = [size(lower_level)[1],1]
#    outer = [size(current_level)[1],1]
#
#    current_level = repeat(current_level,inner=inner)
#    lower_level = repeat(lower_level,outer=outer)
#
#    return hcat(current_level,lower_level)
#
#  end
#
#  return f(stop,start)


end
#



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


function Assign_Value(input_value,std_value)

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
