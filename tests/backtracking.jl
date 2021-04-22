include("../Utils.jl")

using LinearAlgebra
import Random

#A = rand(1000,3000)
##A = rand(3,2,15,6,7)
###
#slice_dim  =2 #slice operation 
##f_dim = (2,5) # f operation 
#f_dim = 1
##
##f(a) = dropdims(sum(a,dims=f_dim),dims=f_dim[1])
#
#X = rand(566,size(A,1))
#
#function f(a) 
#
#	#s = [size(a,i) for i in 1:ndims(a) if !in(i,f_dim)]
#
#
#	return sum(X*a,dims=f_dim)
#	s1,s2 = [size(a,i) for i in f_dim]
#
##	@show s1 s2 
#
#	x = selectdim(a, f_dim[1], [max(1,div(s1,i)) for i in 1:div(s2,2)])
#
##	@show size(x)
#
#	y= selectdim(a,f_dim[2], 1:div(s2,2))
#
##	@show size(y)
#	
##	z = repeat(sqrt.(x) - abs2.(y),3)[1:s2]
#
#	z = sum(sqrt.(x),dims=f_dim) .+ selectdim(abs2.(y),f_dim[2],[3])
#
#
##	@show s1 s2 size(z)
#
#	return sum(a, dims=f_dim)
#
#	return z
#
#
#
#end 
#						
##@show size(A) size(f(A))
##
##
#
#inds = [rand(axes(A,slice_dim),i) for i in [4,5,3,2,5]]
#
##aux(A,f,f_dim,slice_dim,inds)
#
#
#Utils.ApplyF_IndsSets_Check(f, A, slice_dim, inds, f_dim, keepdims=true)
#









L = ["a","b","c","d"]
n=3 
distict = distinct=false 


sol = Utils.Backtracking((L, n, distict), possible_extensions, promising_candidate, accept_sol, output)[1]
@time sol = Utils.Backtracking((L,n, distinct), possible_extensions, promising_candidate, accept_sol, output)[1]

	
	println([L[sol[i]] for i=1:n],"\n")


@show Utils.Random_Items(L,n,distinct=distinct)
@time Utils.Random_Items(L,n,distinct=distinct)


























nothing
