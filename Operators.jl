module Operators

import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
#import DelimitedFiles; const DlmF = DelimitedFiles
#using Einsum

import Utils,TBmodel


# println("Loaded 'Operator' module")


#===========================================================================#
#
# Special operators usually used
# 
#---------------------------------------------------------------------------#




	# ----- probability of localization on individual atoms ----- #

function LDOS(Op=[1];kwargs...)

  return Operator(Op,sum_atoms=false;kwargs...)


end


	# ----- inverse participation ratio ----- #


function IPR(;kwargs...)

  ldos = LDOS(;kwargs...)

  return (P;kwargs...) -> 1 ./ sum(ldos(P).^2,dims=2)
		# special function, cannot be constructed from Operator(...)
		# = 1 if psi localized on a single site, 
		# = nr_atoms for psi spread on all atoms
end



	# ----- expectation of the position operator ----- #


function Position_Expectation(axis::Int64,Rs::Matrix{Float64};
		fpos::Function=identity,nr_at=size(Rs,1),kwargs...)


  return Operator(fpos.(Rs[:,axis]),diag="orbitals",nr_at=nr_at;kwargs...)


end



	# ----- local charge operator ----- #


function Charge(Op_UC,atom;kwargs...)

  return Operator(Op_UC,nr_orb=size(Op_UC,1),purpose="matrixelement",acts_on_atoms=vcat(atom);kwargs...)

end



function Norm(;kwargs...)

  return Operator([1];kwargs...)

end


function Overlap(;kwargs...)

  return Operator([1],purpose="matrixelement";kwargs...)

end


#===========================================================================#
#
# Process input and construct the operator with a specific structure
# 
#---------------------------------------------------------------------------#


function Operator(Op::AbstractArray;purpose="expectation",kwargs...)



  Op,(nr_at,nr_orb,size_H),diag = Understand_OperatorInput(Op;kwargs...)
		# checks if the size of Op coincides with 
		# 	the number of orbitals or atoms 
		#	or the Hamiltonian size 
		# 	and the consistency with 'diag'
		#	diag can be ["orbitals","atoms","all","no"]




  any(purpose.==["expectation","matrixelement"]) || error("'purpose' should be either 'expectation' or 'matrixelement'")



  function choose(operators)

    out = operators[argmax(purpose.==["expectation","matrixelement"])]

    isnothing(out) && error("This operator cannot be used for '"*purpose*"'")
 
    return out

  end



#  println(diag,"  ",purpose,"  ",nr_orb,Op)

  if diag == "no" 
 
#    println("full") 
#    exit()
    return Operator_Full(Op) |> choose

  elseif diag =="all"
 
#    println("1x1")
#    exit()
    return Operator_aNumber(Op,nr_orb=nr_orb,nr_at=nr_at;kwargs...) |> choose
 
 
  elseif diag in ["orbitals","atoms"]

#    println("decomposable")  
#    exit()
    return Operator_Decomposable(Op,diag=diag,nr_orb=nr_orb,nr_at=nr_at;kwargs...) |> choose


  end


end


#===========================================================================#
#
# Operator for the cases 	O = O_{atoms}⊗ 1_{orbitals} 
#			or 	O = 1_{atoms}⊗ O_{orbitals}	
# 
#---------------------------------------------------------------------------#

function Operator_Full(Op::AbstractArray)


  return [(P;kwargs...) -> ExpectationValues_FullOp(P,Op),
	(Pl,Pr=Pl;kwargs...) -> MatrixElements_FullOp(Pl,Pr,Op)]


end


function Operator_Decomposable(Op::AbstractArray;diag=nothing,
		sum_atoms::Bool=true,sum_orbitals::Bool=true,
		acts_on=nothing,
		nr_at=nothing,nr_orb=nothing,Kwargs...)

  if isnothing(diag) || !(diag in ["orbitals","atoms"])
    error("Please specify if the operator is diagonal in atoms or orbitals")
  end

  


  function get_acts_on(acts_on_,diag_,nr_at_,nr_orb_)
  
    !isnothing(acts_on_) && return acts_on_
  
    diag_=="atoms" && !isnothing(nr_at_) && return 1:nr_at_
  
    diag_=="orbitals" && !isnothing(nr_orb_) && return 1:nr_orb_
  
    error("It's unclear on what the operator acts")
  
  end


  function argsinds(acts_on_,diag_,nr_at_,nr_orb_)

  
    diag_=="atoms" && return (1:nr_orb_,acts_on_)

    diag_=="orbitals" && return (acts_on_,1:nr_at_)

  end


  function repeatf(acts_on_,diag_,nr_at_,nr_orb_)

    if diag_=="atoms" 

      return o -> Repeat_Operator_ManyAtoms(o,length(acts_on_),nr_orb_)
     
    elseif diag_=="orbitals" 

      return o -> Repeat_Operator_ManyOrbitals(o,nr_at_,length(acts_on_))

    end
  end




  acts_on = get_acts_on(acts_on,diag,nr_at,nr_orb)
		# if acts_on not specified, all orbitals/atoms are assumed




  inds = TBmodel.Hamilt_indices_all(argsinds(acts_on,diag,nr_at,nr_orb)...,nr_orb;iter=diag)
	# for each item in iter: gives the indices of wf components
	#		on which the operator acts

  all_inds = sort(vcat(inds...))
		# indices for all the wf components on which the operator acts


  Repeat_Many = repeatf(acts_on,diag,nr_at,nr_orb)
		# repeat the operator for several atoms/orbitals 




  if ((diag=="atoms") & sum_atoms) | ((diag=="orbitals") & sum_orbitals)
#  if sum_action

    FullOp = ndims(Op)==1 ? LA.diag(Repeat_Many(LA.diagm(0=>Op))) : Repeat_Many(Op)



    return [(P;kwargs...) -> ExpectationValues(P,FullOp,all_inds),
	(Pl,Pr=Pl;kwargs...) -> MatrixElements(Pl,Pr,FullOp,all_inds)]

  end


  return [(P;kwargs...) -> ExpectationValues(P,Op,inds,iterate=true),
	[(Pl,Pr=Pl;kwargs...) -> MatrixElements_DiagonalOp(Pl,Pr,Op,i) for i in inds]
			]

end



#===========================================================================#
#
# Operator for the case O = 1_{atoms}⊗ 1_{orbitals}
# 
#---------------------------------------------------------------------------#

function Operator_aNumber(Op::AbstractArray;
		sum_atoms::Bool=true,acts_on_atoms=nothing,
		sum_orbitals::Bool=true,acts_on_orbs=nothing,
		nr_at=nothing,nr_orb=nothing,Kwargs...)



  if isnothing(acts_on_orbs) && !isnothing(nr_orb)
    acts_on_orbs = 1:nr_orb
  end

  if isnothing(acts_on_atoms) && !isnothing(nr_at)
    acts_on_atoms = 1:nr_at
  end


  
  acts_on,iter_name = [acts_on_atoms,acts_on_orbs],["atoms","orbitals"]

  good_iter = (!isnothing).(acts_on) 

  acts_on,iter_name = acts_on[good_iter],iter_name[good_iter]

  shorter_iterator = iter_name[argmin(length.(acts_on))]

  inds = TBmodel.Hamilt_indices_all(acts_on_orbs,acts_on_atoms,nr_orb;iter=shorter_iterator)

  all_inds = sort(vcat(inds...))

  sum_iterator = ["atoms","orbitals"][[sum_atoms,sum_orbitals]]




	# if no sum must be performed
  if length(sum_iterator)==0
     return [(P;kwargs...) -> Op[1]*transpose(abs2.(P[all_inds,:])),
		nothing]
		# matrix element doesn't make sense

	# if both sums must be performed   
  elseif length(sum_iterator)==2

    return [(P;kwargs...) -> Op[1]*transpose(sum(abs2.(P[all_inds,:]),dims=1)),
		(Pl,Pr=Pl;kwargs...)-> Op[1]*Pl[all_inds,:]'*Pr[all_inds,:]
		]

  end




  if shorter_iterator == sum_iterator[1]
    
    return [(P;kwargs...) -> 
		Op[1]*transpose(mapreduce(i->abs2.(P[i,:]),+,inds)),
	nothing] 
  
  else         
    return [(P,kwargs...) -> 
		Op[1]*hcat(map(i->sum(abs2.(P[i,:]),dims=1)[:],inds)...),
	nothing]

  end
end



#===========================================================================#
#
# Various ways to act with operators
# 
#---------------------------------------------------------------------------#

function ExpectationValues(P,M,inds=:;iterate=false)



  ndims(M) == 1 && return ExpectationValues_DiagOp(P,M,inds;iterate=iterate) 

  return ExpectationValues_FullOp(P,M,inds;iterate=iterate)

end


function MatrixElements(Pl,Pr,M,inds=:)

  ndims(M) == 1 && return MatrixElements_DiagOp(Pl,Pr,M,inds)
 
  return MatrixElements_FullOp(Pl,Pr,M,inds)

end



# ---- 


function ExpectationValues_DiagOp(P,M,inds=:;iterate=false) 

  !iterate && return transpose(abs2.(P[inds,:]))*reshape(M,:,1)

  return hcat(map(i->ExpectationValues_DiagOp(P,M,i)[:],inds)...)
end


function MatrixElements_DiagOp(Pl,Pr,M,inds=:)

  return Pl[inds,:]'*(reshape(M,:,1).*Pr[inds,:])

end

function MatrixElements_FullOp(Pl,Pr,M,inds=:)

  return Pl[inds,:]'*M*Pr[inds,:]

end





function ExpectationValues_FullOp(P,M,inds=:;iterate=false)

  !iterate && return reshape(LA.diag(MatrixElements_FullOp(P,P,M,inds)),:,1)

  return hcat(map(i->ExpectationValues_FullOp(P,M,i)[:],inds)...)


# diag(P'Op*P) is faster, although <psi_i|Op|psi_j> are computed as well, with i!=j
end





#===========================================================================#
#
# Repeat an orbital-operator for several atoms
#	or an atom-operator for several orbitals
#
#---------------------------------------------------------------------------#


function Repeat_Operator_ManyAtoms(Op,nr_at,nr_orb)
	# Op simulates a Hamiltonian with nr_orb orbitals and one atom
	# 	will be reproduced for nr_at atoms

  i,j,v = SpA.findnz(SpA.sparse(Op))

  indsvals(a) = hcat(TBmodel.Hamilt_indices.([i,j],a,nr_orb)...,v)

  return TBmodel.indsvals_to_Tm(vcat(indsvals.(1:nr_at)...),nr_at*nr_orb)

end


function Repeat_Operator_ManyOrbitals(Op,nr_at,nr_orb)
	# Op simulates a Hamiltonian with nr_at atoms and one orbital
	# 	will be reproduced for nr_orb orbitals




  ai,aj,v = SpA.findnz(SpA.sparse(Op))


  indsvals(o) = hcat(TBmodel.Hamilt_indices.(o,[ai,aj],nr_orb)...,v)


  return TBmodel.indsvals_to_Tm(vcat(indsvals.(1:nr_orb)...),nr_at*nr_orb)

end

#===========================================================================#
#
# Understand how the operator can be simplified and provide the information
# 
#---------------------------------------------------------------------------#



function Understand_OperatorInput(A::AbstractArray;diag=nothing,kwargs...) 
 
  
  ndims(A) in [1,2] || error("Wrong kind of array.")

  if ndims(A)==2 

    if  size(A,1) != size(A,2)

      A = reshape(A,:)

    elseif maximum(abs.(A.-LA.diagm(0=>LA.diag(A)))) < 1e-5

      A = LA.diag(A)

    end
  
  end


  dim = size(A,1)

  numbers = get_nratorbH(;kwargs...)




  conditions = [prod(size(A))==1;
      		(n->!isnothing(n) && n==dim).(numbers)]

  if any(conditions)

    meaning = ["all","orbitals","atoms","no"][argmax(conditions)]

    if meaning != Utils.Assign_Value(diag,meaning)

      error("Strange operator input.")

    end
   
    return A,numbers,meaning

  end

  return A,numbers,diag
 



#  conditions = (n->!isnothing(n) && n==dim).(numbers)
#  conditions = [conditions;[all(isnothing.(numbers)),prod(size(A))==1]]
#  meanings = ["orbitals","atoms","no","no","all"][conditions]


#  length(meanings) == 0 && error("Strange operator input.")


#  (isnothing(diag) || diag==meanings[1]) && return A,numbers,meanings[1]


end

function get_nratorbH(;nr_at=nothing,nr_orb=nothing,size_H=nothing,kwargs...)


  nr_vars = count((!isnothing).([nr_at,nr_orb,size_H]))

  if nr_vars == 1 & isnothing(nr_orb)

    nr_orb = 1
    nr_vars = 2

  end


  nr_vars < 2 && return [nr_at,nr_orb,size_H]

  nr_vars==3 && nr_at*nr_orb!=size_H && error("Wrong number of atoms/orbitals")

  all((!isnothing).([nr_at,nr_orb])) && return [nr_at,nr_orb,nr_at*nr_orb]
  all((!isnothing).([nr_at,size_H])) && return [nr_at,div(size_H,nr_at),size_H]
  all((!isnothing).([nr_orb,size_H])) && return [div(size_H,nr_orb),nr_orb,size_H]

end



#===========================================================================#
#
# Normalize wavefunctions on columns
# 
#---------------------------------------------------------------------------#


function N(P)

  return P ./ reshape(LA.norm.(eachcol(P)),1,:)

end



#===========================================================================#
#
# Trace over orbitals/atoms with some projection operator
# 
#---------------------------------------------------------------------------#


function Trace_Orbitals(M::AbstractVector;Op=[1],kwargs...)

  nr_at,nr_orb, = get_nratorbH(;size_H=size(M,1),kwargs...)

  i(a) = TBmodel.Hamilt_indices(1:nr_orb,a,nr_orb)

  which_len = (length(Op) .== [nr_at,1,nr_orb,nr_at*nr_orb])

  any(which_len) || error("Wrong operator length")

  oper = [a->Op[a],a->Op,a->Op,a->Op[i(a)]][findfirst(which_len)]

  return [sum(M[i(a)].*oper(a)) for a in 1:nr_at]

end

function Trace_Orbitals(M::AbstractMatrix;Op=[1],kwargs...)

	return Trace_Orbitals(LA.diag(M);Op=Op,kwargs...)

end

function Trace_Atoms(M::AbstractVector;Op=[1],kwargs...)

  nr_at,nr_orb, = get_nratorbH(;size_H=size(M,1),kwargs...)

  i(o) = TBmodel.Hamilt_indices(o,1:nr_at,nr_orb)

  which_len = (length(Op) .== [nr_at,1,nr_orb,nr_at*nr_orb])

  any(which_len) || error("Wrong operator length")

  oper = [o->Op,o->Op,o->Op[o],o->Op[i(o)]][findfirst(which_len)]

  return [sum(M[i(o)].*oper(o)) for o in 1:nr_orb]


end


function Trace_Atoms(M::AbstractMatrix;Op=[1],kwargs...)

	return Trace_Atoms(LA.diag(M);Op=Op,kwargs...)

end







#############################################################################









end