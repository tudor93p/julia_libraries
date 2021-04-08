#nr_leads = 3


lead_sizes = [2,3,2][1:nr_leads]

leadlabels = ["A","B","C"][1:nr_leads]

LeadHeads = [[[0 div(Ny,2)] ; [0 div(Ny,2)+1]], [ [Nx+1 2] ; [Nx+1 3] ; [Nx+1 4]],  [[div(Nx,2) 0];[div(Nx,2)+1 0]] ][1:nr_leads]
										# on the surface already

gfs = [r(size(LH,1)*d0) for LH in LeadHeads]

# each dict is a natural user minimal input
userleads = [
				Dict(
					:label => leadlabels[i],
					:coupling => HoppMatr, 

					:head => [LeadHeads[i]],


					:intercell => [HoppMatr(LeadHeads[i]) + r(lead_sizes[i]*d0)*0.1],
					
					:intracell =>	[HoppMatr(LeadHeads[i])],

					:GF => E-> [gfs[i] +  E*LA.I*0.1],

							)
				for i in 1:nr_leads]

