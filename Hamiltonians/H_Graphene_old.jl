###############################################################################
# This file refers to the model in the seminal paper of Kane and Mele
#   PRL 95, 146802 (2005)  
#   "Z2 Topological Order and the Quantum Spin Hall Effect"
###############################################################################


module H_Graphene

import LinearAlgebra; const LA = LinearAlgebra
#import PyCall

import Utils,TBmodel,Algebra,BandStructure


export Hopping_function#,Dirac_Points


println("Loaded H_Graphene")  



#===========================================================================#
#
# Finds the Dirac point(s)
#
#---------------------------------------------------------------------------#

function Dirac_Point(H,K,nr_bands=nothing,tol=1e-8,window=0.05)

  !isnothing(nr_bands)  && (nr_bands *= 3)

  EnK = BandStructure.get_eigen(H,tol=tol,nr_bands=nr_bands)(K).Es

  EnK = EnK[abs.(EnK) .< window]

  D = Algebra.OuterDist(EnK,EnK)

  
	# Catesian Index [CI] containing the two closest eigenvalues,


  indices = findall(LA.triu(trues(size(D)...),1))
#  indices = findall(LA.triu(trues(length(EnK),length(EnK)),1))


#  indices[argmin(D[indices])]

  D =  D[indices] 

  CIs = indices[findall(abs.(D .- minimum(D)) .< tol)]

   
 
#  println(findall(LA.triu(abs.(D.-minimum(D[LA.triu(trues(size(D)...),1)])) .< tol,1)))


#  CIs = findall(LA.triu(abs.(D.-minimum(D[findall(LA.triu(D.>0))])) .< tol,1))

  return [EnK[i] for i in Tuple(CIs[1])] |> E -> real(sum(E)/length(E))

#  if length(CIs) in [1,2]
#
#    return [EnK[i] for p in map(Tuple,CIs) for i in p] |> E -> sum(E)/length(E)
#
#  end

#  for w in range(0,maximum(abs.(EnK)),length=100)
#
#  EnK = EnK[abs.(EnK) .< window]
#
#  d1,d2 = 
#	# Catesian Indices CIs
#
#  for maxendiff in collect(100:-5:1)*minimum()
#
#    CIs = findall(LA.triu(D.<maxendiff,1))
#
#
#    !mute && println(length(CIs), "diffs ",[-(EnK[collect(Tuple(ci))]...) for ci in CIs])
#
#
#  end

  error("The Dirac points couldn't be found/too many found.")




end

#===========================================================================#
#
# Describes hoppings between two atoms in the lattice, according to the paper
#
#---------------------------------------------------------------------------#


function Hopping_function(param_H_,pyLatt;dist_tol=1e-5)

  default_param_H = (
	Intralayer_Hopping	= 1.0,
	Interlayer_Hopping	= 0.0,
	Interlayer_HoppDecay	= 1.0/50.0,
	Layer_Bias		= 0.0,
	Lattice_Imbalance	= 0.0,
#        KaneMele_SOC		= 0.0,
#	Rashba_SOC		= 0.0,
	AntiHaldane		= 0.0,
	Hopp_Cutoff 		= 1e-6,
					)






  !issubset(keys(param_H_),keys(default_param_H)) && error("The parameters must be among "*join(map(string,keys(default_param_H)),", ")*".")

  param_H = merge(default_param_H,param_H_)


  dist = pyLatt.Distances(nr_neighbors=2,sublattices="1")[2:3]#::Vector{Float64}

  hopp_cutoff = param_H.Hopp_Cutoff



        # ------------------ all atoms --------------------------- #      


  nr_uc = 1 # minimum, provided the kind of hoppings


  UCs = pyLatt.get_UCs(nr_uc)

  Atoms = NamedTuple{(:A,:B)}(
	map(s->Algebra.FlatOuterSum(pyLatt.PosAtoms(s),UCs),["A","B"]) )


  AllAtoms = vcat(values(Atoms)...)


  if pyLatt.VectDim == 3
    zmax,zmin = AllAtoms[:,3] |> x-> (maximum(x),minimum(x))

    zsep = zmax-zmin

  else
    zmax=zmin=zsep=0.0

  end
 


  hoppings = []


        # -------------- nearest neighbor ------------------------ #      

  push!(hoppings, TBmodel.Hopping_Term( 
		(ri,rj) -> ((length(ri) == 2) || isapprox(ri[3],rj[3],atol=dist_tol)) && isapprox(LA.norm(ri.-rj),dist[1],atol=dist_tol), 
		param_H[:Intralayer_Hopping],
		tol = hopp_cutoff))



   
        # -------------- interlayer hopping ---------------------- #      

  pyLatt.VectDim == 3 && push!(hoppings, TBmodel.Hopping_Term(
		param_H[:Interlayer_Hopping],
		(ri,rj) -> abs2(ri[3]-rj[3])/(LA.norm(ri.-rj)+hopp_cutoff/10.0)^2*exp(-(LA.norm(ri.-rj)-zsep)/param_H.Interlayer_HoppDecay),
		tol = hopp_cutoff))


        # -------------- layer bias (electric field) ------------- #      

  pyLatt.VectDim == 3 && push!(hoppings, TBmodel.Hopping_Term( 
		(ri,rj) -> isapprox(ri,rj,atol=dist_tol),
		param_H[:Layer_Bias],
		(ri,rj) -> ri[3],
		tol = hopp_cutoff))




        # -------------- atom bias (lattice imbalance) ----------- #      


  D1(ri,Rs) = Algebra.FlatOuterDist(hcat(ri...),Rs)

  atom_bias(ri) = [-1,1][[any(D1(ri,Rs).< dist_tol) for Rs in Atoms]][1]



  push!(hoppings, TBmodel.Hopping_Term( 
		(ri,rj) -> isapprox(ri,rj,atol=dist_tol),
		param_H[:Lattice_Imbalance],
		(ri,rj) -> atom_bias(ri),
		tol = hopp_cutoff))


        # -------------- anti-Haldane nnn imaginary hopping ------ #      



  nn(ri) = findall(abs.(D1(ri,AllAtoms) .- dist[1]) .< dist_tol)

  middle(ri,rj) = AllAtoms[intersect(nn(ri),nn(rj))[1],:]

  chirality(ri,rj) = sign(Algebra.cross(eachcol(hcat(ri,rj) .- middle(ri,rj))...)[3])



  push!(hoppings, TBmodel.Hopping_Term( 
		(ri,rj) -> ((length(ri) == 2) || isapprox(ri[3],rj[3],atol=dist_tol)) && isapprox(LA.norm(ri.-rj),dist[2],atol=dist_tol), 
		param_H[:AntiHaldane],
		(ri,rj) -> 1im*atom_bias(ri)*chirality(ri,rj),
		tol = hopp_cutoff))



        # -------------- Kane-Mele spin-orbit coupling ----------- #      

#    if same(la.norm(dif[:2]),dist_nnn)*same(abs(dif[2]),0.)*neighbors.size == 0
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:KaneMele_SOC],
#		,
#		tol = hopp_cutoff))

        # -------------- Rashba spin-orbit coupling -------------- #      
#
#  push!(hoppings, TBmodel.Hopping_Term( 
#		(ri,rj) -> isapprox(LA.norm(ri.-rj),dist_nnn,atol=dist_tol),
#		param_H[:Rashba_SOC],
#		,
#		tol = hopp_cutoff))
#


  
  hopping_fij = Utils.Sum_functions(hoppings...)
  


  return nr_uc,40,hopping_fij

        # -------------- distance cutoff ------------------------- #      

  function find_dist_cutoff()


    min_dist = max(maximum(dist),zsep)



    
    for dist in Iterators.countfrom(min_dist,dist[1]/10)

      ri = [0.0,0.0]

      rj = [dist,0.0]

      (pyLatt.VectDim == 3) && ( ri = vcat(ri,zmax) )
      (pyLatt.VectDim == 3) && ( rj = vcat(rj,zmin) )



      if isapprox(maximum(abs.(hopping_fij(ri,rj))),0.0,atol=hopp_cutoff)

        dist_xyz = LA.norm(ri.-rj)
        dist_xy  = LA.norm((ri.-rj)[1:2])


        nr=last(last(pyLatt.Number_of_neighbors(dist_xy,tol=-log10(dist_tol))))
#        HL = PyCall.pyimport("Lattices").Honeycomb_Lattice()
#        nr = last(last(HL.Number_of_neighbors(dist_xy,tol=-log10(dist_tol))))

        return (Float64(dist_xyz),max(Int64(nr),50))


      end
    end

    error("The distance cutoff couldn't be calculated\n")

  end


  dist_cutoff, max_neighbors = find_dist_cutoff()

#  println("dist cutoff: ",dist_cutoff,",     neighbors: ",max_neighbors,",   hopp cutoff: ",hopp_cutoff)




		# it works both ways, but somehow the second is faster

#  hopping_fij = TBmodel.Sum_Hoppings_Cutoffs(hoppings,dist_cutoff,hopp_cutoff)
  hopping_fij = TBmodel.Hopping_Cutoffs(hopping_fij,dist_cutoff,hopp_cutoff)

  return nr_uc,max_neighbors,hopping_fij




#  dist,rs_AB = param_latt
#  z,dist_nn,*rest = dist
#  (dist_nn) = param_latt
#  sigma = Algebra.PauliMatrices()



#  def atom_type(ri):
# 
#    for rA,rB in zip(*rs_AB):
#
#      if same(ri,rA): return 'A'
#      if same(ri,rB): return 'B'

#  def bond_index(v):
#    theta = np.arctan2(v[1],v[0])
#    return np.argmax([same(np.abs(theta), angle) for angle in [np.pi/6., np.pi/2., 5*np.pi/6.]])
#    
#  def sxdz(dd):
#    dd = np.append(dd,0)/la.norm(dd)
#
#    return np.cross(np.array([sx,sy,sz]),dd,axisa=0,axisb=0,axisc=0)[2]
            

#  def nu(ri,rj):
#    
#    for rq in rs_tot:
#      vj = np.append(rj-rq,0)
#      vi = np.append(rq-ri,0)
#      if same(la.norm(vi),dist_nn) and same(la.norm(vj),dist_nn):
#
#        return np.sign(np.cross(vi,vj)[2])
#
#    print('Error: intermediate atom not found')
#    print(ri,rj)
#    return 0


#  def f0(ri,rj):
#    				# staggered potential
#      if   atom_type(ri) == 'A': return +lv*s0
#      elif atom_type(ri) == 'B': return -lv*s0

#  def f1(ri,rj):
		  # nearest neighbor hopping and Rashba SOC
#    return t[bond_index(rj-ri)]*s0 + 1.j*Rashba*sxdz(rj-ri)

#  def f2(ri,rj):
		  # next nearest neighbor hopping (SOC)
#    return 1.j*lSO*nu(ri,rj)*sz


#  return [f0,f1,f2]

#end


# nr_TB : number of UCs around the central UC
# AtomPos: position of atomns in the original uc

#def extra_param_KM(distances,SCell_SLatt,nr_TB,Lattice0):
#    
#  aux, AtomPos = Lattice0
#
#  ucs_in_UC, [LattVect,aux] = SCell_SLatt
#
#  UCs = np.matmul(vectors_of_integers(len(LattVect),nr_TB),LattVect)
#
#  ucs = np.array([U+u for U in UCs for u in ucs_in_UC])
#
#  rs_AB = [ucs + r for r in AtomPos]
#  
#        # distances between atoms (0, nn, nnn)
#        # the A and the B atoms considered
#  return [distances,rs_AB]
end


###############################################################################

end
