module Algebra

#using Distributed 

import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
#import PyCall


import Utils

gethostname()=="tudor-HP" && import Dierckx,FFTW

#export #PauliMatrices,OuterSum




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

function fft(x::AbstractMatrix, w::Number=0; dim::Int=1, addup=true)
	
	E = exp.(-1im*w*(axes(x,dim).-1))

	Y = Utils.multiply_elementwise(E, x, dim)

	!addup && return Y
	
	return dropdims(sum(Y, dims=dim),dims=dim)

end 


function fft(x::AbstractVector, w::Union{Number,AbstractVector}=2pi*(axes(x,1).-1)/length(x); addup=true, kwargs...)

	E = exp.(-1im*OuterBinary(w, axes(x,1).-1, *))

	addup && return E*x

	return E .* reshape(x,1,:)

end 



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Interp1D(x, y, k::Int)

	if k==0

		get_ind = Interp1D(vcat(x...), 1:length(x), 1)

		out(X::Number) = y[Int(round(get_ind(X)))]

		out(X::AbstractArray) = y[Int.(round.(get_ind(X)))]
		
		return out 

	end 


	return Dierckx.Spline1D(vcat(x...), vcat(y...), k=k)

end 

function Interp1D(x, y, k::Int, X)

	Interp1D(x,y,k)(X)

end



#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function Mean(A::AbstractArray, dims::Union{Tuple,Int})

	dims_ = unique(vcat(dims...))

	i = [d in dims_ ? 1 : Colon() for d in 1:ndims(A)]

	return sum(A, dims=dims)[i...]/prod(size(A)[dims_])

end


function Mean(A::AbstractArray, dims::Nothing=nothing)

	sum(A)/length(A)

end




#===========================================================================#
#
# Linear independence of arrays
#
#---------------------------------------------------------------------------#


function Nr_LinearIndependentArrays(array_list)

	LA.rank(Array(hcat(reshape.(array_list,:)...)))

end

function LinearIndependentArrays(array_list)
	
	maxrank = Nr_LinearIndependentArrays(array_list)
	
	indep(inds) = Nr_LinearIndependentArrays(array_list[inds])==length(inds)


	function add_new(new_ind=2, inds=[1]) 

		(length(inds)==maxrank || new_ind==length(array_list)+1) && return inds 

		for i in inds

			!indep([i, new_ind]) && return add_new(new_ind+1, inds)
			
		end

		indep([inds;new_ind]) && return add_new(new_ind+1, [inds;new_ind])

		return add_new(new_ind+1, inds)

	end

	return matrix_list[add_new()]

end

#===========================================================================#
#
# Modified Gram-Schmidt
#
#---------------------------------------------------------------------------#

function GramSchmidt(V::T, iter_axis=1; 
										 						normalize=true, tol=1e-10) where T

	inds(i) = insert!(repeat(Any[:], ndims(V)-1), iter_axis, i)

#	inds(i) = [j==iter_axis ? i : Colon() for j in 1:ndims(V)]


	iter = if T<:AbstractArray{<:Number}

							A->eachslice(A; dims=iter_axis)

					elseif T<:AbstractVector{<:AbstractArray{<:Number}}

							identity

					end 


	U = copy(V)*0.0

	for (j,Vj) in enumerate(iter(V))
		
		u = copy(Vj)

		for i in 1:j-1

			u -= (Ui->LA.dot(Ui,u)*Ui)(U[inds(i)...])

		end 

		norm = LA.norm(u)

		norm<tol && error("Invalid input matrix for Gram-Schmidt!")

		normalize && setindex!(U, u/norm, inds(j)...)

	end


  return U

end





#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


function Commutator(A,B)

	A*B-B*A

end




#===========================================================================#
#
# Orthonormal eigenvectors of a Hermitian matrix
#			grouped according to their energy (degenerate eigenvectors together)
#
#---------------------------------------------------------------------------#

function eigenvectors_degenerate_sectors(H; tol=1e-10)

	eigen = LA.eigen(LA.Hermitian(Array(H)))


	return [GramSchmidt(eigen.vectors[:,s],2) 
									for s in Utils.IdentifySectors(eigen.values; tol=tol)]

end

#===========================================================================#
#
#	Simultaneously diagonalize a set of commuting Hermitian matrices
#
#---------------------------------------------------------------------------#
	
project(Ms,U) = [U'*M*U for M in Ms]

function SimultDiagonaliz_CommutingHMatrices(Hs; tol=1e-10)

	all(length.(Hs).==1) && return Utils.UnitMatrix(1)

	projectors = eigenvectors_degenerate_sectors(Hs[1], tol=tol)


	length(Hs)==1 && return projectors



	return vcat(map(projectors) do p 
	
		np = SimultDiagonaliz_CommutingHMatrices(project(Hs[2:end], p); tol=tol) 

		return [p*new_p for new_p in np]

	end...)


end




#===========================================================================#
#
# Structure factors for a set of matrices
#
#---------------------------------------------------------------------------#

function structure_factors(Ls, test=true)

	Utils.isList(Ls, AbstractMatrix) || error("List of matrices needed")

	F = zeros(Complex{Float64}, length(Ls), length(Ls), length(Ls))

	S = inv([LA.dot(L1,L2) for L1 in Ls, L2 in Ls])

	for (a,La) in enumerate(Ls)
		
		for (b,Lb) in enumerate(Ls[1:a-1])

			Cab = Commutator(La, Lb)

			F[a,b,:] = S*[-1im*LA.dot(L, Cab) for L in Ls]

			F[b,a,:] = -F[a,b,:]

		end

	end


	if test 

		e = map(eachcol(rand(axes(Ls,1),2,10))) do (a,b)
			
			LA.norm(Commutator(Ls[a],Ls[b]) .- im*mapreduce(prod,+,zip(F[a,b,:],Ls)))
		
		end

		maximum(e) > 1e-10 && println("####### Wrong structure factors")

	end

	return F

end






#===========================================================================#
#
# Computes the Killing form of a set of matrices 
# 													(or from the structure factors directly)
#
#---------------------------------------------------------------------------#

function Killing_form(arg::T) where T
	
	F = eachslice(if T<:AbstractArray{<:Number,3}

												arg
									
								elseif Utils.isList(T, AbstractMatrix) 

												structure_factors(arg)
								
								else 
												error("Wrong input for Killing form")
								end,

								dims=1)

	return [sum(Fa.*transpose(Fb)) for Fa in F, Fb in F]

end



#===========================================================================#
#
#	The center of a Lie algebra
#
#---------------------------------------------------------------------------#

function Center_LieAlgebra(Ls, K=Killing_form(Ls))

	Utils.isList(Ls, AbstractMatrix) || error("Wrong input") 

	return [mapreduce(prod,+,zip(l,Ls)) for l in eachcol(LA.nullspace(K))]

end


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function ConvoluteVectorPacket(weights, values, centers, 
															 delta, vectors::AbstractMatrix;
															 dim=1,
															 get_weights=false, normalize=true, 
															 keepdims=true, kwargs...) 

	W = getCombinedDistrib(weights, values, centers, delta;
												 normalize=normalize)

	out = dim |> function multipl(dim)
	
					v = sum(
									
							if dim==1 
						
						 W*view(vectors, axes(W,2),:)

							elseif dim==2 
								
								view(vectors, :, axes(W,2))*transpose(W)

							end,

							dims=dim)

							return keepdims ? v : dropdims(v, dims=dim)
					
					end 

	return get_weights ? (out,W) : out

end



#===========================================================================#
#
# Find root using bisection method
#
#---------------------------------------------------------------------------#

function FindRoot_Bisection(f,u10,u20=nothing;val10=nothing,val20=nothing,a=0,b=a+1e-8,move_right=x->x+1,move_left=x->x-1)

  function bisection(u1,u2)


 #   println("\nbisection")

 #   println("bisection ",u1," ",u2)

    um = (u1+u2)/2

    val = f(um)

    if val < a 
      return bisection(um,u2)

    elseif val > b
      return bisection(u1,um)

    else 
      return um

    end

  end



  uvalu(u) = (u,f(u))
 
# --- Find the interval, if not provided (properly) --- #

  function interval(u1,u2,val1,val2)

    (val1 < a) && (val2 > b) && return bisection(u1,u2)

#      println("bisection")
#      return bisection(u1,u2)
#    end

    if val1 > a

      if val1 > b
        (u2,val2) = (u1,val1)
      end

      (u1,val1) = uvalu(move_left(u1))
    end

    if val2 < b 

      if val2 < a
        (u1,val1) = (u2,val2)
      end

      (u2,val2) = uvalu(move_right(u2))
    end


#    println("interval ",u1," ",u2)

    return interval(u1,u2,val1,val2)

  end



  if isnothing(val10)

    val10 = f(u10)

  end




  if (isnothing(u20)  || isapprox(u10,u20,atol=1e-12))

    return interval(u10,u10,val10,val10)

  

  end

  return interval(u10,u20,val10,isnothing(val20) ? f(u20) : val20)

  


end





#===========================================================================#
#
# Partially invert matrix
#
#---------------------------------------------------------------------------#


function Cofactor(A,row,col)

  return  (-1)^(col + row) * LA.det(A[axes(A,1).!=row,axes(A,2).!=col])

end



function Partial_Inverse(A,elements;compute="some")

  if compute == "some" 

    a,b = 1e-10,1e+10

    detA = LA.det(A)

    val = abs(detA)



    if a < val < b # && return [Cofactor(A,j,i) for (i,j) in elements]/detA

      return Dict([(i,j)=>Cofactor(A,j,i)/detA for (i,j) in elements])

    end


    Q1 = FindRoot_Bisection(u->abs(LA.det(A*u)),1,a=a,b=b,val10=val,
			move_right=x->x*2, move_left=x->x/2)

    Q2 = Q1^(1/(1-1/size(A,1)))

#    return A*Q2 |> B->[Cofactor(B,j,i) for (i,j) in elements]/LA.det(Q1*A)

    detA, B = LA.det(Q1*A), Q2*A

    return Dict([(i,j) => Cofactor(B,j,i)/detA  for (i,j) in elements])


#    return Partial_Inverse(A,elements,compute="full")


  elseif compute == "full" 

    return inv(Array(A)) |> B -> Dict([i=>B[i...] for i in elements])


  else
    
   error("'compute' must be 'full' or 'some'")

  end

end




#function Different_ijs(A,B=A)

#  A .!= reshape(B,1,A)
#
# 
#  for (i,a) in enumerate(A)
#
#    B[axes(B,1) .!= i]
# 
#
#  (j,
#
#  i,j = axes(a,1),axes(a,1)
#
#  i,j = repeat(i,inner=length(j)),repeat(j,outer=length(i))
#
#  w = (i.!=j)
# 
#
#

  
#  return [[y,x] for (x,y) in Base.product(B,A) if x!==y]




#end

#===========================================================================#
#
#  
#
#---------------------------------------------------------------------------#


#===========================================================================#
#
# simulate localized wavefunctions 
#
#---------------------------------------------------------------------------#

#function Localized_WF(atoms,weight_orb,centers=atoms,localiz_shape="Gaussian",spread=1)
#
#  if ndims(centers) == 1 
#    centers = reshape(centers,1,:)
#  end
#
#  nr_orb = size(weight_orb(atoms[1,:]),1)
#
#  psi0 = zeros(Complex{Float64},size(atoms,1)*nr_orb,size(centers,1))
#
#  weights = get_Distribution(localiz_shape)
#
##  dist = OuterDist(atoms,centers)
#
#
#  W2 = weights.(eachcol(atoms),eachcol(centers),[spread])
#
#  W2 = Normalize_Columns((*).(W2...))
##  W2[i,j] : weight on atom i for the wf with center j
#
#
#  for (iatom,(atom,w2)) in enumerate(zip(eachrow(atoms),eachrow(W2)))
#
#    w = Normalize(weight_orb(atom)).*reshape(w2,1,:)
#
#    psi0[TBmodel.Hamilt_indices(1:nr_orb,iatom,nr_orb),:] = w
#
#  end
#
#  return psi0
#
#end

#===========================================================================#
#
# project a single wavefunction to certain eigenfunctions around E0 
#
#---------------------------------------------------------------------------#

function Project_WF_toE0(eigenen,eigenwf,newE0,newE0wid;projection_method="impose_weights",energy_shape="Gaussian")


  function change_coef()

    weights = get_Distribution(energy_shape)(eigenen,newE0,newE0wid)

    if projection_method == "impose_weights"

      return coef -> (coef./abs.(coef)) .* weights
    
    elseif projection_method == "partial_spectrum"
   
      weights2 = Complex{Float64}.(weights.>maximum(weights,dims=1)/100)

      return coef -> coef .* weights2
    
    end
    
  end


  newcoef = change_coef()

  return wf -> eigenwf*Normalize_Columns(newcoef(eigenwf'*reshape(wf,:,1)))


end





function diff_vc(v::Tv, c::Tc) where {Tv,Tc}

	all(T->T<:AbstractArray,(Tv,Tc)) && return reshape(c,1,:) .- reshape(v,:,1)

	return c .- v

	error()

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

struct myDistrib

	F::Function 

	function myDistrib(f::Function; normalize=false) 
	
		hasmethod(f, (Real, Real)) || error("Wrong input function")
	
		g(x::Real, w::Real)::Real = f(x,w)
	
		g(v::Real, c::Real, w::Real)::Real = f(c-v, w)
	
		g(x::AbstractArray, w::Real)::Array = f.(x,w)
	
		g(v::AbstractArray, c::Real, w::Real)::Array = f.(c.-v, w)
	
		g(v::Real, c::AbstractArray, w::Real)::Array = f.(c.-v, w)
	
		g(v::AbstractVector, c::AbstractVector, w::Real)::Matrix = f.(reshape(c,1,:) .- reshape(v,:,1), w)
	
	
		function g(w::Real)::Function 
			
			h(x::Real)::Real = f(x,w)	
	
			h(v::Real, c::Real)::Real = f(c-v, w)
	
			h(x::AbstractArray)::Array = f.(x,w)
	
			h(v::AbstractArray, c::Real)::Array = f.(c.-v, w)	
	
			h(v::Real, c::AbstractArray)::Array = f.(c.-v, w)
	
			h(v::AbstractVector, c::AbstractVector)::Matrix = f.(reshape(c,1,:) .- reshape(v,:,1), w)	
	
	
			return h 
	
		end 

#		g(a::Tuple) = g(a...)
	
		return new(g)
#	return g
	
	end 
	
	
	function myDistrib(str::String; kwargs...)
	
		E = Meta.parse(str); 
		
		return myDistrib(@eval (x::Real, w::Real) -> $(E); kwargs...)
	
	end 


end 
# f = @eval ... ;return myDistrib(f)
#  Remembers only the last f and overwrites the previous ones
#
# myDistrib(@eval ...) might create strange "var is not defined" errors


#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#


normDistrib(W::Real)::Real = W 
normDistrib(W::AbstractVector)::Real = sum(W) + 1e-20
normDistrib(W::AbstractMatrix)::Matrix = sum(W, dims=2) .+ 1e-20

normDistrib(F::Function)::Function = normDistrib ∘ F
normDistrib(D::myDistrib)::Function = normDistrib(D.F)


normalizeDistrib(W::Real, N::Union{Float64,Int64})::Real = W/N 
normalizeDistrib(W::AbstractVector, N::Union{Int64,Float64})::Vector = W/N 
normalizeDistrib(W::AbstractMatrix, N::AbstractMatrix)::Matrix = W./N 


function normalizeDistrib(W::Tw, N::Tn, n::Bool) where Tw<:T where Tn<:T where T<:Union{Real, AbstractVecOrMat}

	n ? normalizeDistrib(W, N) : W 

end 


function normalizeDistrib(W::T) where T<:Union{Real, AbstractVecOrMat}

	normalizeDistrib(W, normDistrib(W))

end 


function normalizeDistrib(W::T, n::Bool) where T<:Union{<:Real, <:AbstractVecOrMat}

	n ? normalizeDistrib(W) : W 

end 


normalizeDistrib(F::Function)::Function = normalizeDistrib ∘ F 
normalizeDistrib(F::Function, n::Bool)::Function = n ? normalizeDistrib(F) : F

normalizeDistrib(D::myDistrib)::Function = normalizeDistrib(D.F) 
normalizeDistrib(D::myDistrib, n::Bool)::Function = normalizeDistrib(D.F, n)





function normalizedDistribAndNorm(W::Union{Real, AbstractVecOrMat}, n::Bool=true)::Tuple

	N = normDistrib(W) 

	return (normalizeDistrib(W, N, n), N)

end 


function normalizedDistribAndNorm(F::Function, n::Bool=true)::Function 

	FN(args...) = normalizedDistribAndNorm(F(args...), n)

end 

function normalizedDistribAndNorm(D::myDistrib, n::Bool=true)::Function 

	normalizedDistribAndNorm(D.F, n)

end 





function (D::myDistrib)(args...; normalize::Bool=false) 

	normalizeDistrib(D.F(args...), normalize)

end 

#===========================================================================#
#
# Lorentzian, Gaussian functions etc
#
#---------------------------------------------------------------------------#


Lorentzian = myDistrib("w / (x^2 + w^2) / pi")

Gaussian = myDistrib("exp(-x^2/w^2)/(w*sqrt(pi))")

Heaviside = myDistrib("Float64(x>=w)")

Rectangle = myDistrib("Float64(-w/2<=x<=w/2)")


#w=rand()+0.1;x=rand(100,100);
#@show LA.norm(Gaussian(x,w) - Lorentzian(x,w))
#
#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function getDistrib(args::Tuple; kwargs...)#::Function

	getDistrib(args...; kwargs...)

end 

function getDistrib(name::String, args...; kwargs...)#::Function
	
	getDistrib(Symbol(name), args...; kwargs...)

end 

function getDistrib(F::Function; normalize=false)::Function
	
	normalizeDistrib(F, normalize)

end

function getDistrib(name::Symbol, args...; kwargs...) 

	getDistrib(getproperty(@__MODULE__, name), args...; kwargs...)
	
end 


function getDistrib(F::Function, args...; normalize=false) 
	
	normalizeDistrib(F(args...), normalize)

end 

function getDistrib(D::myDistrib, args...; kwargs...)

	getDistrib(D.F, args...; kwargs...)

end 






function getDistribAndNorm(args...; kwargs...)

	normalizedDistribAndNorm(getDistrib(args...; normalize=false);
													 kwargs...)

end 






#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function getCombinedDistrib(args; normalize=false) 

	length(args)==1 && return getDistrib(args[1]...; normalize=normalize)

	D(a) = getDistrib(a..., normalize=false)

	D1 = D(args[1]) # may be a function or numeric


	if typeof(D1)<:Union{Real,AbstractArray}

		return normalizeDistrib(
							mapreduce(D, Utils.multiply_elementwise, args[2:end], init=D1),
							normalize)


	elseif D1 isa Function 

		Ds = [D1; D.(args[2:end])]

		totalD(Args...) = mapreduce(d->d(Args...), Utils.multiply_elementwise, Ds)

		return normalizeDistrib(totalD, normalize) 

	else 

		error("Wrong D1")

	end 

end 


function getCombinedDistribAndNorm(args; normalize=false)

	normalizedDistribAndNorm(getCombinedDistrib(args, normalize=false), 
													 normalize)

end 


#	If length('args')>1, the quantities must be given in the proper order. 
#	Either of:
#				names, deltas 
#				names, values, deltas 
#				names, values, centers, deltas 


mktuple(a) = isa(a,Tuple) ? a : tuple(a)

function getCombinedDistrib(args...; kwargs...)

	getCombinedDistrib(Utils.Zipmap(mktuple, args); kwargs...)

end


function getCombinedDistribAndNorm(args...; kwargs...)
	
	getCombinedDistribAndNorm(Utils.Zipmap(mktuple, args); kwargs...)

end 



#===========================================================================#
#
# Normalize an array (full, on columns, on rows)
#
#---------------------------------------------------------------------------#


function Normalize(A::AbstractArray, p::Real=2; dims=nothing, tol=1e-12)

	isnothing(dims) && return A/LA.norm(A, p)

	N = LA.norm.(eachslice(A, dims=dims), p)

	N[N.<tol] .= tol

  return A./reshape(N, [d==dims ? Colon() : 1 for d in 1:ndims(A)]...)
										
end

Normalize_Columns(A::AbstractMatrix, p::Real=2) = Normalize(A, p, dims=2)
Normalize_Rows(A::AbstractMatrix, p::Real=2) = Normalize(A, p, dims=1)


#===========================================================================#
#
# Pauli Matrices
#
#---------------------------------------------------------------------------#



function PauliMatrices()::Dict{Int64,Array{Complex{Float64},2}}

  s0 = Complex{Float64}[[1 0]; [0 1]]
  s1 = Complex{Float64}[[0 1]; [1 0]]
  s2 = Complex{Float64}[[0 -1im]; [1im 0]]
  s3 = Complex{Float64}[[1 0]; [0 -1]]
  
  
  return Dict(0=>s0, 1 => s1, 2=>s2, 3=>s3)
end

function SpinMatrices(s)
	
    """
    Construct spin-s matrices for any half-integer spin.

    Parameters
    ----------

    s : float or int
        Spin representation to use, must be integer or half-integer.
    include_0 : bool (default False)
        If `include_0` is True, S[0] is the identity, indices 1, 2, 3
        correspond to x, y, z. Otherwise indices 0, 1, 2 are x, y, z.

    Returns
    -------

    ndarray
        Sequence of spin-s operators in the standard spin-z basis.
        Array of shape `(3, 2*s + 1, 2*s + 1)`, or if `include_0` is True
        `(4, 2*s + 1, 2*s + 1)`.
    """

	d = Int(2*s+1)

	Sz = 1/2 * LA.Diagonal(d-1:-2:-d)
												
  # first diagonal for general s from en.wikipedia.org/wiki/Spin_(physics)
	
	diag = [1/2*sqrt(i*(2s+1-i)) for i in 1:d-1]

	Sx = LA.diagm(1=>diag,-1=>diag)
	
	Sy = im*LA.diagm(1=>-diag,-1=>diag)

	return Dict(0=>Utils.UnitMatrix(d), 1=>Sx, 2=>Sy, 3=>Sz)

end




#===========================================================================#
#
# solve Quadratic equation
#
#---------------------------------------------------------------------------#

function QuadraticEq(a::Number,b::Number,c::Number)::Vector

  ( -b .+ [1,-1]*sqrt(b^2-4*a*c) ) ./ (2*a)

end



#===========================================================================#
#
# Cross product with vectors length 2
#
#---------------------------------------------------------------------------#

function cross(A::AbstractArray,B::AbstractArray)::AbstractArray


  good(r) = length(r) < 3 ? vcat(r,zeros(3-length(r))) : r[1:3]

  return LA.cross(good(A),good(B))


end


#===========================================================================#
#
# Dot product with vectors of different lengths
#
#---------------------------------------------------------------------------#

function dot(A::AbstractArray,B::AbstractArray)::Number

  i = 1:min(length(A),length(B))

  return LA.dot(A[i],B[i])

end

#===========================================================================#
# 
# Efficient outer sum of two lists of vectors U = 
#	Returns a list [X,Y,...], such that U[i]+V[j] = X[i,j] + Y[i,j] + ...
#
#---------------------------------------------------------------------------#



function OuterBinary(U::Tu, V::Tv, op::Function; flat=false)::AbstractArray where Tu<:T where Tv<:T where T<:Union{Number,AbstractVector}
	
	right_shape(a::Number) = hcat(a)
	right_shape(a::AbstractVector) = reshape(a,:,1)

	return dropdims(OuterBinary(right_shape(U), right_shape(V), op; flat=flat),
									dims=3-flat)

end 


function OuterBinary(U::Tu, V::Tv, op::Function; flat=false, dim=1) where Tu<:T where Tv<:T where T<:AbstractMatrix

	dim2 = 3 - dim # columns if dims=rows and reversed


	sizes = map(enumerate(zip(size(U),size(V)))) do (d,(su,sv))

				if d==dim 

					return flat ? su*sv : (su,sv)

				else  
	
					su!=sv && error("Wrong sizes")

					return su 

				end 

			end 



	typ = Base.return_types(op, typeof.(first.((U,V)))) |> function (ts)

					!isempty(ts) ? ts[1] : typeof(op(U[1],V[1])) # Any

		end 


	out = similar(Array{typ}, Utils.flat(sizes)...)
	
	i0 = fill(Colon(), 2-flat)


#	for (i,(u,v)) in enumerate(zip(eachslice(U,dims=dim2), eachslice(V,dims=dim2)))


	for i in 1:sizes[dim2]

		setindex!(out,
							view(op.(selectdim(U, dim2, i),
											 reshape(selectdim(V, dim2, i), 1, :)
											 ), i0...),
							[i0,i][dim]..., [i0,i][dim2]...)
		
	end 
	
	return out 

end


OuterSum(args...; kwargs...) = OuterBinary(args..., +; kwargs...)

OuterDiff(args...; kwargs...) = OuterBinary(args..., -; kwargs...)

function OuterDist(args...; kwargs...) 
	
	D = get(kwargs, :dim, 1)==1 ? 3 : 1

	return selectdim(sum(abs2, OuterDiff(args...; kwargs...), dims=D),
									 D, 1) .|> sqrt  
									
end 

#===========================================================================#
# 
# Flattened outer sum, each line S[k] is a sum of U[i_k] and V[j_k]
#	S[k] = U[i] + V[j], with (i,j) = unravel_index(k,(len(U),len(V)))
#
#---------------------------------------------------------------------------#
  

FlatOuterSum(args...; kwargs...) = OuterBinary(args..., +; kwargs..., flat=true)

FlatOuterDiff(args...; kwargs...) = OuterBinary(args..., -; kwargs..., flat=true)

function FlatOuterDist(args...; kwargs...) 
	
	LA.norm.(eachslice(FlatOuterDiff(args...; kwargs...),
										 dims=get(kwargs,:dim,1)))

end 


#eachslice(D,dims=dim)
#
#
#def FlatOuter_IndexConvert(U,V):
#
#  get_k = lambda ij: np.ravel_multi_index(ij,(len(U),len(V)))
#
#  get_ij = lambda k: np.unravel_index(k,(len(U),len(V)))
#
#  return get_k,get_ij
#
#



#===========================================================================#
#
# Distance equals 
#
#---------------------------------------------------------------------------#

function EuclDistEquals(d0; tol=1e-8)

	isd0(dist) = isapprox(dist, d0, atol=tol)

	vtm(A::AbstractVector) = reshape(A,1,:)
	vtm(A::AbstractMatrix) = A

	return function(A::a, B::a) where {a<:T,b<:T} where T<:AbstractVecOrMat

		a<:AbstractVector && b<:AbstractVector && return isd0(LA.norm(A-B))

		return isd0.(OuterDist(vtm(A),vtm(B)))

	end


end 




#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#

function get_Bonds(atoms::AbstractMatrix, bondlength::Number; 
									 tol=1e-8, kwargs...)

	get_Bonds(atoms, EuclDistEquals(bondlength; tol=tol); kwargs...)

end 

function get_Bonds(atoms::AbstractMatrix, isBond::Function; 
									 inds=true, pos=false, as_matrices=false)


  N = 5000    # for memory reasons

	nr_at = size(atoms,1)

	nr_batches = Int(ceil(nr_at/N))

	batch_size = Int(ceil(nr_at/nr_batches))

	batches = [(k-1)*batch_size+1:min(k*batch_size,nr_at) for k=1:nr_batches]



	function get_pairs(b,c)
		
		d = isBond(atoms[batches[b],:], atoms[batches[c],:])

		x0 = batch_size*([b,c].-1)

		return sort([Tuple(sort(Tuple(x).+x0)) 
										for x in findall(b!=c ? d : LA.triu(d,1))])

	end
	
	bond_indices = sort(vcat([get_pairs(b,c) for b=1:nr_batches for c=1:b]...))

	
	if as_matrices

		return get_Bonds_toMatrix(atoms, bond_indices; inds=inds, pos=pos)

	end 

	!pos && return bond_indices

	bond_Rs = [[atoms[i,:] for i in ab] for ab in bond_indices]

	return !inds ? bond_Rs : (bond_indices, bond_Rs)



	







	#	Rs[pair_index][member_index][coordinate_index]

#	Rs = atoms[hcat(vcat.(b...)...),:]
#	Rs[pair_index, member_index, coordinate_index]


	
#	return np.transpose(Rs[pairs,:dim],axes=(0,2,1)).reshape(-1,dim)


end


function get_Bonds_toMatrix(X; inds=false, pos=false)

	xor(inds,pos) || error("Specify either 'inds' or 'pos'")


	if inds 

		bond_indices = X 

		M_bond_indices = zeros(Int, length(bond_indices), 2)


		for (i,(a,b)) in enumerate(bond_indices)

			M_bond_indices[i, :] = [a,b]

		end 

		return M_bond_indices



	elseif pos 

		bond_Rs = X


		M_bond_Rs = zeros(Float64, length(bond_Rs), sum(length.(bond_Rs[1])))


		for (i,Rs) in enumerate(bond_Rs)

			M_bond_Rs[i] = vcat(Rs...)

		end 

		return M_bond_Rs

#	return !inds ? M_bond_Rs : (M_bond_indices, M_bond_Rs)
	end 

	return 

end 




function get_Bonds_toMatrix(atoms, bond_indices; inds=false, pos=false)

	M_bond_indices = get_Bonds_toMatrix(bond_indices; inds=true)

	if !pos 

		inds && return M_bond_indices

	else 
	
		M_bond_Rs = hcat(atoms[M_bond_indices[:,1],:], atoms[M_bond_indices[:,2],:])

		inds && return (M_bond_indices, M_bond_Rs)

		return M_bond_Rs

	end 

	return 

end 





function get_Bonds_fromMatrix(M; inds=false, pos=false)

	xor(inds,pos) || error("Specify either 'inds' or 'pos'")
	
	inds && return Tuple.(eachrow(convert(Array{Int},M)))


	pos && return map(eachrow(M)) do Rij 

						n = div(length(Rij),2)

						return [Rij[1:n],Rij[n+1:end]]
				
				end 

end 	


#############################################################################

end
