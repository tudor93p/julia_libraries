module Algebra

import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
#import PyCall
import Utils

#export #PauliMatrices,OuterSum







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


	return [Algebra.GramSchmidt(eigen.vectors[:,s],2) 
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

function get_CombinedDistribution(values, centers, delta;
																	weights = "Lorentzian",
																	normalize=true, get_norm=false)

	choose(A, i=1) = typeof(A)<:Tuple ? A[min(i,end)] : A

	len(A) = typeof(A)<:Tuple ? length(A) : 1

	imax = maximum(len.((weights, values, centers, delta)))

	function get_w(i=1; kwargs...)

		get_Distribution(choose(weights, i); kwargs...)(
								choose(values, i), choose(centers, i), choose(delta, i))
	end

	imax==1 && return get_w(;normalize=normalize, get_norm=get_norm)

	w = mapreduce(i->get_w(i; normalize=false), (w1,w2)->w1.*w2, 1:imax)

	N = 1e-20 .+ sum(w, dims=2)

	return normalize ? w./N : w |> W -> get_norm ? (W,N) : W

end




function ConvoluteVectorPacket(values, centers, 
															 delta, vectors::AbstractMatrix;
															 weights="Lorentzian", get_weights=false) 

	W = get_CombinedDistribution(values, centers, delta; weights=weights) 

	return W*vectors[axes(W,2),:] |> v -> get_weights ? (v,W) : v

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




#===========================================================================#
#
# Lorentzian, Gaussian functions
#
#---------------------------------------------------------------------------#


function Gaussian(value,centers,width)

  return exp.(-abs2.(reshape(centers,1,:) .- value)/width^2)/(width*sqrt(pi))

end

function Lorentzian(value,centers,delta)

  return 1/pi*delta./(abs2.(reshape(centers,1,:) .- value) .+ delta^2)

end

function Heaviside(value,centers,delta=0.0)

  return 1.0*((reshape(centers,1,:) .- value) .<= delta)

end


function Rectangle(value,centers,delta)

  return Heaviside(-value,-centers,delta/2) .* Heaviside(value,centers,delta/2) 
  
#  return (-delta/2 .<= (reshape(centers,1,:) .- value) .<= delta/2)

end


 

function get_Distribution(name; normalize=true, get_norm=false)

  options = ["Lorentzian", "Gaussian","Rectangle"]

  functions = [Lorentzian,Gaussian,Rectangle]

  name in options || error("Please choose: "*join(options,", "))

  f = functions[findfirst(name .== options)]

	!normalize && !get_norm && return f


	fnormf(args...) = f(args...) |> fa -> (fa, 1e-20 .+ sum(fa,dims=2) )


	!normalize && get_norm && return fnormf

	!get_norm && return (args...) -> fnormf(args...) |> fn -> fn[1]./fn[2]

	return (args...) -> fnormf(args...) |> fn -> (fn[1]./fn[2],fn[2])


end

 



#===========================================================================#
#
# Normalize an array (full, on columns, on rows)
#
#---------------------------------------------------------------------------#


function Normalize(x)

  return x/LA.norm(x)

end



function Normalize_Columns(A)

  return A./reshape(LA.norm.(eachcol(A)),1,:)

end


function Normalize_Rows(A)

  return A./reshape(LA.norm.(eachrow(A)),:,1)

end


#function MultiplyRows(matrix,vector)
#
#  return matrix .* reshape(vector,:,1)
#
#end
#
#
#function MultiplyCols(matrix,vector)
#
#  return matrix .* reshape(vector,1,:)
#
#end

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

function OuterBinary(U,V,op)


  if !(ndims(U) in (1,2)) || !(ndims(V) in (1,2))
    error("The arrays must be 1D or 2D!")
  end

	ndims(U) == 1 && return OuterBinary(reshape(U,:,1),V,op)
	
	ndims(V) == 1 && return OuterBinary(U,reshape(V,:,1),op)


#  if  size(U,2) != size(V,2) # will produce an error
#  end

#  return [(op).(u,transpose(v)) for (u,v) in zip(eachcol(U),eachcol(V))]

	
  return [(op).(u,reshape(v,1,:)) for (u,v) in zip(eachcol(U),eachcol(V))]

end

  
function OuterSum(args...)

  return OuterBinary(args...,+)

end

function OuterDiff(args...)

  return OuterBinary(args...,-)

end


function OuterDist(args...)

  return sqrt.(sum([abs2.(dXi) for dXi in OuterDiff(args...)]))

end


#===========================================================================#
# 
# Flattened outer sum, each line S[k] is a sum of U[i_k] and V[j_k]
#	S[k] = U[i] + V[j], with (i,j) = unravel_index(k,(len(U),len(V)))
#						(true at least for python)
#---------------------------------------------------------------------------#

function FlatOuterSum(args...)
  
  return hcat([reshape(s,length(s)) for s in OuterSum(args...)]...)

end
 

function FlatOuterDiff(args...)

  return hcat([s[:] for s in OuterDiff(args...)]...)

end

function FlatOuterDist(args...)

  return OuterDist(args...)[:]


end


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







#############################################################################

end
