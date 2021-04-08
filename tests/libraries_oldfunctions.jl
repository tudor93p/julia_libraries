
#function GF_Decimation(Energy,H_intracell,H_intercell,NrLayers,gRightLeft,include_leads=2)
#				"""
#	Energy -> value of the energy at which the GF is evaluated
#
#	Inhomogeneous hopping:
#
#		H_intracell(n) -> intra-layer Hamiltonian of layer n
#
#		H_intercell(n,m) -> inter-layer hopping between layer n and layer m
#
#	NrLayers -> upper bound of n,m  (lower bound is 1)
#
#	gRightLeft = gR,gL -> functions which give the GF of the semi-infinite leads
#
#				"""
#
#	
#	function coupling(source::Int64,dest::Int64)
#
#		dest > source && return (H_intercell(dest,source),GR(dest,dest))
#
#		dest < source && return (H_intercell(dest,source),GL(dest,dest))
#
#
#	end
#
#	function G_lead_layer(lead,m,getGF)
#
#		lead > NrLayers && return GR(lead,m+1)*H_intercell(m+1,m)*getGF(m,m)
#
#		lead < 1 && return GL(lead,m-1)*H_intercell(m-1,m)*getGF(m,m)
#
#
#	end
#
#
#
#	gRs,gLs = [g(Energy) for g in gRightLeft]
#	
#
#	# methods to calculate #
#	
#	
#	function GL_(n,m)
#
#		n==m && return	GF(Energy,H_intracell(n),coupling(n,n-1))
#
#		n==0<m && return	G_lead_layer(n,m,GL)
#
#		return nothing
#
#	end
#
#	function GR_(n,m)  
#
#		n==m && return GF(Energy,H_intracell(m),coupling(m,m+1))
#
#		n==NrLayers+1>m && return G_lead_layer(n,m,GR)
#
#		return nothing
#	end
#
#
#
#	# combine sweeps into actual GFs 
#
#	function G_(n,m=n)
#
#		n==m>NrLayers && return GF(GR(n,n),coupling(n,n-1))
#
#		n==m<1 && return GF(GL(n,n),coupling(n,n+1))
#
#		n==m && return GF(Energy,H_intracell(n),coupling(n,n-1),coupling(n,n+1))
#
#		(( n<1 ) | (n>NrLayers )) && return G_lead_layer(n,m,G)
#
#		abs(n-m)==1 && return G(n)*H_intercell(n,m)*(m>n ? GR : GL)(m,m)
#
#		n in [0,NrLayers+1] && return G_lead_layer(n,m,G)
#
#	end
#
#
#	# calculate *and* store #
#
#	GLd,GRd,Gd = Dict(),Dict(),Dict()
#
#	GL(n,m) = n==m<1 ? gLs[min(end,1-n)] : get!(GLd,(n,m)) do  
#																						GL_(n,m) end
#
#	GR(n,m) = n==m>NrLayers ? gRs[min(end,n-NrLayers)] : get!(GRd,(n,m)) do 
#																						GR_(n,m) end
#
#
#	G(n,m=n) = get!(Gd,(n,m)) do 
#																G_(n,m) end
#
#
#
#	return GL,GR,G
#
#end
#

#	for (n,m) in zip(1:NrLayers+1+include_leads,NrLayers:-1:0-include_leads)

#	  GLd[(n,n)] = GF(Energy,h(n)+U(n,n-1)*GL(n-1,n-1)*U(n-1,n))

#		GLd[(n,n)] = GF(Energy,h(n),(U(n-1,n),GL(n-1,n-1)))

#	  GLd[(0,n)] = GL(0,n-1)*U(n-1,n)*GL(n,n)

#	  GRd[(m,m)] = GF(Energy, h(m) + U(m,m+1)*GR(m+1,m+1)*U(m+1,m))

#		GRd[(m,m)] = GF(Energy,h(m),(U(m+1,m),GR(m+1,m+1)))

#	  GRd[(NrLayers+1,m)] = GR(NrLayers+1,m+1)*U(m+1,m)*GR(m,m)

#	end

#


#===========================================================================#
#
# Reconstruct the full GF matrix 
#
#---------------------------------------------------------------------------#


#function fullGF_fromGD(GD,NrLayers;iLayer,nr_orb,nr_at=nothing,include_leads=false)
#
#  iterLayers = include_leads ? collect(0:NrLayers+1) : collect(1:NrLayers)
#
#  if isnothing(nr_at)
#    nr_at = sum(length.(iLayer.(iterLayers)))
#  end
#
#  out = SpA.spzeros(Complex{Float64},nr_at*nr_orb,nr_at*nr_orb)
#
#  i(n) = TBmodel.Hamilt_indices_all(1:nr_orb,iLayer(n),nr_orb,flat=true) 
#
#  for n in iterLayers
#
#    for (a,b) in [(n,n),(n-1,n),(n,n+1)]
#      
#      if (a in iterLayers) && (b in iterLayers)
#
#        out[i(a),i(b)] = GD(a,b)
#
#      end
#
#    end
#
#  end
#
#  return out 
#
#end
#

#===========================================================================#
#
#
#
#---------------------------------------------------------------------------#





