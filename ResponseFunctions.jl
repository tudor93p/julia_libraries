module ResponseFunctions




#using Distributed

#using BenchmarkTools

import Dates
import LinearAlgebra; const LA = LinearAlgebra
import SparseArrays; const SpA = SparseArrays
#import DelimitedFiles; const DlmF = DelimitedFiles
import GreensFunctions
import Utils,Algebra

#?export *#Hamilt_indices#Bloch_Hamilt,Hopping_Term

#println("Loaded ResponseFunctions")


#===========================================================================#
#
# Josephson current from Furusaki, Physica B 203 (1994) 214-218 -- Eq. (5)
# 		see also Asano, PRB 63, 052512 (2001) -- Eq. (20)
#				PRB 74, 064507 (2006) -- Eq. (27) 
#
#---------------------------------------------------------------------------#

function JosephsonCurrent(GD,i;f=LA.tr)

  return -1im*f(GD(i,i-1) - GD(i-1,i))

  # GD = Green's function with Matsubara frequency wn
  # Furusaki:-- "Although the summation over wn in Eq. (5) i
  #		is originally from n=-oo to n=oo, we can alternatively sum up
  #		over positive wn only and take the real part of the sum. 
  #		Thus, we assume wn > 0 in the following discussion."
  # 	     -- "(...) calculate the current (...) in the normal region,
  #		[where] the electric charge is always conserved. 
  #		If the current is calculated in a superconducting region, 
  #		the source term proportional to the order parameter 
  #		must be included."

#  Im = sum(abs.(imag(out)) ./ (abs.(out) .+1e-12))/length(out)

#  Im > 1e-6 && println("J has imaginary part. mean(imag/abs)=",round(Im,digits=7))

# return vcat(real(out)...)
end

#===========================================================================#
#
# Tunneling conductance according to PRB 86, 174512 (2012),
#	formula based on Lee-Fisher: PRL 47, 882 (1981)
#
#---------------------------------------------------------------------------#

function TunnelingConductance_LeeFisher(GD,i,j=i;f=LA.tr)


	G(n,m) = (g -> (g' - g)/(2im))(GD(n,m))

  return f(	G(i,j)*G(j-1,i-1) 
									+ G(i-1,j-1)*G(j,i) 
									- G(i,j-1)*G(j,i-1) 
									- G(i-1,j)*G(j-1,i)
					)
	# prefactor aside, Lee-Fisher formula (3) with j0->i, j0'->j
	# becomes formula (15) in PRB with i=j=x+1


#	Im = sum(abs.(imag(La.diag(out))) ./ (abs.(LA.diag(out)) .+1e-12))/length(LA.diag(out))

#  Im > 1e-6 && println("LeeFisher T has imaginary part. mean(imag/abs)=",round(Im,digits=7))

#  return vcat(real(out)...)

end

#===========================================================================#
#
# Quantum conductance
#
#---------------------------------------------------------------------------#


function CaroliConductance(G1,source,drain,G2=G1;f=LA.tr)

#  sp(A) = SpA.dropzeros(SpA.sparse(A))
#
	GammaS = GreensFunctions.DecayWidth(source)
	
	GammaD = GreensFunctions.DecayWidth(drain)


#  out = f(LA.diag(sp(GammaS)*sp(G1)*sp(GammaD)*sp(G2)'))
##  out = LA.tr(GammaS*G1*GammaD*G2')
	
	return f(GammaS*G1*GammaD*G2')

#  Im = sum(abs.(imag(out)) ./ (abs.(out) .+1e-12))/length(out)

#  Im > 1e-6 && println("Caroli T has imaginary part. mean(imag/abs)=",round(Im,digits=7))

#  any(abs.(imag(out)) .> 1e-9) && println("T has imaginary part")


#  return real(out)

end




#	A = G*W*G'

#	GammaS = GreensFunctions.DecayWidth(source)
	
#	GammaD = GreensFunctions.DecayWidth(drain)

#	return f(GammaS*G1*GammaD*G2')


#h = H ( (layer2,atom2)-> (layer1,atom1))





function BondTransmission(G, SE_lead, Bonds, RBonds, Hoppings; f=LA.tr, kwargs...)

	W = GreensFunctions.DecayWidth(SE_lead)

	atoms = sort(unique(vcat(vcat.(Bonds...)...)))

#	Gd = Dict(a=>G(a) for a in atoms)


	bondT = zeros(length(Bonds))


	for sector in Utils.IdentifySectors(first.(sort(Bonds)))

		Gi = G( Bonds[sector[1]][1] )

		for bond_index in sector 

			Gj = G( Bonds[bond_index][2] )

			bondT[bond_index] = -2imag(f(	
															  Gi*W*Gj'*Hoppings[bond_index]
															 -Gj*W*Gi'*Hoppings[bond_index]'
															 ))

#			bondT[bond_index] = Tij_minus_Tji

		end

	end 

	return bondT

#	return [-2imag(f(Gd[i]*W*Gd[j]'*h)) for ((i,j),h) in zip(Bonds,Hoppings)]
	
#		T[i,j] = -2imag(f(G(i)*W*G(j)'*h))

#		T[i,j] = f((G(i)*W*G(j)').*(H(j,i))) # ? Probably not element-wise



#	nr_at = 

																					
end


function SiteTransmission(BondT, Bonds, RBonds)

	nr_at = length(unique(vcat(vcat.(Bonds...)...)))

	siteT = zeros(nr_at * length(RBonds[1][1]) )

	for ((i,j),(Ri,Rj),Tij) in zip(Bonds, RBonds, BondT)

		siteT[ i .+ [0,nr_at] ] +=  Tij * (Rj-Ri)

		siteT[ j .+ [0,nr_at] ] += -Tij * (Ri-Rj)
		
	end

	return siteT 

end



#===========================================================================#
#
# Charge response function (density-density)
#
#---------------------------------------------------------------------------#


function DensityDensity(charge_oper::Function,
			ene_chi_i::AbstractVector{Float64},
			wf_chi_i::Matrix{Complex{Float64}},
			ene_chi_j::AbstractVector{Float64},
			wf_chi_j::Matrix{Complex{Float64}},
			ene_chi::AbstractVector{Float64}, 	# Energy for response
			ichi::Int64=1, jchi=1; 		# matrix element of chi
			temperature_chi::Float64=1e-7,	# enters the occupation factors
			delta::Float64=1e-4,		# analytic cntinuation
			tol::Float64=1e-7
					)::Matrix{Complex{Float64}}

  jchi = vcat(jchi...)


  charge_oper_i  = charge_oper(ichi)
  charge_oper_j  = charge_oper.(jchi)




#  first_matrix_element = conj(wf_chi_i[ichi,:]) .* wf_chi_j[ichi:ichi,:]
  println("computing matrix elements")

@time  first_matrix_element = charge_oper_i(wf_chi_i,wf_chi_j)


  

  occupation_factor = reshape(Occupation.(ene_chi_j,temperature_chi),1,:) .- Occupation.(ene_chi_i,temperature_chi) 




@time  ME_charge_oper_j = map(o->o(wf_chi_j,wf_chi_i),charge_oper_j)
		# too large


  second_matrix_element(a,b) = [ME[b,a] for ME in ME_charge_oper_j]
 
#  second_matrix_element(a,b) = conj(wf_chi_j[jchi,b]) .* wf_chi_i[jchi,a] 




  denominator(a,b) = (ene_chi_i[a] - ene_chi_j[b] + 1im*delta) .- reshape(ene_chi,1,:)



  chi_contribution(a,b) = occupation_factor[a,b]*first_matrix_element[a,b]*second_matrix_element(a,b)./denominator(a,b)


  chi_total = zeros(Complex{Float64},length(jchi),length(ene_chi))

#  test = (abs.(occupation_factor) .> tol)
  test = (abs.(occupation_factor.*first_matrix_element) .> tol)

  nr = sum(test) 

  start = Dates.now()

  println()
#  @show count(test)

  nprints = 10
  
  for (q,(a,b)) in enumerate(Tuple.(findall(test)))

    if q%div(nr,nprints)==0
  
      tot= Dates.value(Dates.now()-start)/1000*nprints

      tot = tot/60 > 1 ? [round(tot/60,digits=1),"minutes"] : [round(tot,digits=1),"seconds"]

      println(q,"/",nr,"/",prod(size(test)),"  total: ",join(tot," "))

      start = Dates.now()
    end


      chi_total .+= chi_contribution(a,b)

  end

  return chi_total

end


#===========================================================================#
#
# Occupation number
#
#---------------------------------------------------------------------------#








function Occupation(energy::Float64,temperature::Float64)::Float64
			# gives the occupation at a given temperature

  energy < -temperature && return 1.0
			# if below threshold

  energy > temperature && return 0.0
			# if above threshold

  return 0.5 - 0.5*energy/temperature
			# if in interval
end


















#===========================================================================#
#
# Bare translation of the charge response function from fortran90
#
#---------------------------------------------------------------------------#

#
## norbitals is number of spinless orbitals
## i and j correspond to different spin flavours!!!
#function elementchi_old(	wf_chi_i,ene_chi_i,
#			wf_chi_j,ene_chi_j,
#			ene_chi,
#			ichi, jchi,
#			temperature_chi, delta)
#			
##subroutine elementchi(wf_chi_i,ene_chi_i,&
##                           wf_chi_j,ene_chi_j, &
##                           ene_chi, &
##                           ichi, &
##                           jchi, &
##                           temperature_chi, &
##                           delta, &
##                           chi_total, &
##                           norbitals, &
##                           num_wf_i, &
##                           num_wf_j, &
##                           num_ene_chi)
##  implicit none
##  !!!!!!!!!!!!!!!!!!!!!!!!!!!
##  ! input variables
##  !!!!!!!!!!!!!!!!!!!!!!!!!!!
##
##  integer, intent(in) :: num_wf_i,num_wf_j ! dimensions of wavefunctions
##  integer, intent(in) :: norbitals ! number of operators
##  integer, intent(in) :: ichi,jchi ! matrix elements
##  integer, intent(in) :: num_ene_chi ! number of energies
##  complex (kind=8),intent(in) :: wf_chi_i(num_wf_i,norbitals)  ! WF i
##  complex (kind=8),intent(in) :: wf_chi_j(num_wf_j,norbitals)  ! WF j
##  real (kind=8),intent(in) :: ene_chi_i(num_wf_i)   ! Energy i
##  real (kind=8),intent(in) :: ene_chi_j(num_wf_j)   ! Energy j
##  real (kind=8),intent(in) :: ene_chi(num_ene_chi)   ! Energy for response
##  real (kind=8),intent(in) :: temperature_chi   ! Energy for response
##  real (kind=8),intent(in) :: delta   ! analytic cntinuation
##
##  complex (kind=8),intent(out) :: chi_total(num_ene_chi)  ! WF i
##
##
##  ! index for energy wave i and wave j
##  integer :: ie,iwf1,iwf2
##  real (kind=8) :: enetmp,etmp1,etmp2 ! temporal energies
##  ! temporal wavefunctions
##  complex (kind=8) :: wftmp1(norbitals),wftmp2(norbitals)
##  complex (kind=8) :: chitmp  ! temporal chi
##  complex (kind=8) :: ieps,im  ! smearing of Chi
##  complex (kind=8) :: vawwbv,holdab  ! matrix element of the response
##  real (kind=8) :: occ_fac,occ1,occ2 ! occupation factor 
##  real (kind=8) :: den_res ! denominator of the response function
##
##  im = (0.d00,1.d00)
##
##
##  chi_total = 0.d00
#  chi_total = zeros(Complex{Float64},length(ene_chi))
##  ieps = im*delta  ! complex infinitesimal
#  ieps = 1im*delta
##
#
#
#  for (iwf1,(etmp1,wftmp1)) in enumerate(zip(ene_chi_i,eachrow(wf_chi_i)))
##    do iwf1 = 1, num_wf_i ! first loop over states
##      etmp1 = ene_chi_i(iwf1)
##      wftmp1(:) = wf_chi_i(iwf1,:)
#
#
##    for (etmp2,wftmp2) in zip(ene_chi_j,eachrow(wf_chi_j))
#    for (iwf2,(etmp2,wftmp2)) in enumerate(zip(ene_chi_j,eachrow(wf_chi_j)))
#
#
##    for (iwf2,etmp2) in enumerate(ene_chi_j)
##      do iwf2 = 1, num_wf_j  ! second loop over states
##        etmp2 = ene_chi_j(iwf2)
#
##        ! fermi energy has been put in 0
##        call occupation(etmp1,temperature_chi,occ1)
##        call occupation(etmp2,temperature_chi,occ2)
##        occ_fac = occ2-occ1  ! occupation factor
#
#      occ_fac = occupation(etmp2,temperature_chi) - occupation(etmp1,temperature_chi)
#
##!        if (dabs(occ_fac).lt.1.d-04)  cycle  ! next iteration if too far
##        ! if contribution is different from zero continue
##        wftmp2(:) = wf_chi_j(iwf2,:)
##! compute element of linear response
##    ! calculate only at zero temperature
##        ! calculate the matrix elementes <A><B>
#
#
#
#      vawwbv = conj(wftmp1[ichi])*wftmp2[ichi] 		# first matrix element
#      vawwbv = vawwbv * conj(wftmp2[jchi])*wftmp1[jchi]	# second
#
##          vawwbv = conjg(wftmp1(ichi))*wftmp2(ichi) ! first matrix element
##          vawwbv = vawwbv * conjg(wftmp2(jchi))*wftmp1(jchi) ! second
#
##  ! save in more accesible variables
#
#      
#      for (ie,enetmp) in enumerate(ene_chi)
#
#       
##          do ie = 1, num_ene_chi
##            enetmp = ene_chi(ie)
##            den_res = ((etmp1-etmp2) - enetmp)
#        den_res = (etmp1-etmp2) - enetmp
#
##            chitmp = occ_fac*vawwbv/(den_res+ieps)  ! add contribution
#        chitmp = occ_fac*vawwbv/(den_res+ieps)  #! add contribution
#
##            ! add to the chi matrix array
#        chi_total[ie] += chitmp
#
#
#
##            chi_total(ie) = chi_total(ie) + chitmp 
##          enddo ! close loop over energies
##    enddo  ! close loop over energies
##  enddo ! close loop over the different matrices of linear response
#
#      end
# 
#    end
#  end
#
#  return chi_total
##  return
##end subroutine elementchi
#
#end
#



#############################################################################

end
