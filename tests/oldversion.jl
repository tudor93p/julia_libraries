
LeftLead = get(VirtualLeads,:LeftLead,nothing)
RightLead = get(VirtualLeads,:RightLead,nothing)

function h(n)

	n<1 && return LeftLead[:intracell][min(end,1-n)]

	n>NrLayers && return RightLead[:intracell][min(end,n-NrLayers)]

	return Hscatt(n)

end



function U(n,m)

	n>m && return U(m,n)'

	#:coupling defined from the lead to the layer
	#:intercell is defined in the lead direction (opposite)

	n==0==m-1 && return LeftLead[:coupling]
	n==N==m-1 && return RightLead[:coupling]'

	n<m<1 && return LeftLead[:intercell][min(1-m,end)]'
	N<n<m && return RightLead[:intercell][min(n-N,end)]

	return Hscatt(n,m)

end


function oldGFDecim(Energy)





function GF_DecimationOld(Energy,H_intracell,H_intercell,NrLayers;gRight=nothing,gLeft=nothing)


	function exists(n)

		1<=n<=NrLayers && return true

		n<1 && !isnothing(gLeft) && return true

		n>NrLayers && !isnothing(gRight) && return true

		return false


	end

	
	function coupling(source::Int64,dest::Int64)

		!exists(dest) && return ()

		return ((H_intercell(dest,source),(dest > source ? GR : GL)(dest,dest)),)

	end



	function G_lead_layer(lead,m,getGF)

		lead > NrLayers && return GR(lead,m+1)*H_intercell(m+1,m)*getGF(m,m)

		lead < 1 && return GL(lead,m-1)*H_intercell(m-1,m)*getGF(m,m)


	end


	gRs = isnothing(gRight) ? [nothing] : [g(Energy) for g in gRight]

	gLs = isnothing(gLeft) ? [nothing] : [g(Energy) for g in gLeft]


	
	function GL_(n,m)

		n==m && return	RGF(Energy,H_intracell(n),coupling(n,n-1)...)

		n==0<m && return	G_lead_layer(n,m,GL)

		return nothing

	end


	function GR_(n,m)  

		n==m && return RGF(Energy,H_intracell(m),coupling(m,m+1)...)

		n==NrLayers+1>m && return G_lead_layer(n,m,GR)

		return nothing
	end



	# combine sweeps into actual GFs 

	function G_(n,m=n)

		n==m>NrLayers && return RGF(GR(n,n),coupling(n,n-1)...)

		n==m<1 && return RGF(GL(n,n),coupling(n,n+1)...)

		n==m && return RGF(Energy,H_intracell(n),coupling(n,n-1)...,coupling(n,n+1)...)

		(( n<1 ) | (n>NrLayers )) && return G_lead_layer(n,m,G)

		abs(n-m)==1 && return G(n)*H_intercell(n,m)*(m>n ? GR : GL)(m,m)

		n in [0,NrLayers+1] && return G_lead_layer(n,m,G)

	end


	# calculate *and* store #

	GLd,GRd,Gd = Dict(),Dict(),Dict()


	GL(n,m) = n==m<1 ? gLs[min(end,1-n)] : get!(GLd,(n,m)) do  
																						GL_(n,m) end

	GR(n,m) = n==m>NrLayers ? gRs[min(end,n-NrLayers)] : get!(GRd,(n,m)) do 
																						GR_(n,m) end


	G(n,m=n) = get!(Gd,(n,m)) do 
																G_(n,m) end



	return GL,GR,G

end





	
return GF_DecimationOld(Energy,h,U,NrLayers,
												gRight = isnothing(RightLead) ? nothing : RightLead[:GF],
												gLeft = isnothing(LeftLead) ? nothing : LeftLead[:GF],
												
												)



end

#########################################################################

