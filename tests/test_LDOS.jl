#!/home/tudor/apps/julia/julia-1.1.1/bin/julia
import LightGraphs,MetaGraphs
import Graph,Algebra,Utils,Random


const LG = LightGraphs
const MG = MetaGraphs

using TikzGraphs
using TikzPictures
import TBmodel
include("./../GreensFunctions.jl")
import LinearAlgebra;const LA=LinearAlgebra

Random.seed!(1234);

r(n,m=n) = rand(Complex{Float64},n,m)


#===========================================================================#
#
#	Atom list & bond function
#
#---------------------------------------------------------------------------#

Nx = 10

Ny = 5

dist_tol=1e-5

Atoms = Algebra.FlatOuterSum(hcat(zeros(Ny),1:Ny),hcat(1:Nx,zeros(Nx)))


function bond(A::AbstractArray,B::AbstractArray)

	ndims(A)==ndims(B) == 1 && return isapprox(LA.norm(A-B),1,atol=dist_tol)

	ndims(A)==1 && return bond(hcat(A...),B)
	ndims(B)==1 && return bond(A,hcat(B...))

	return isapprox.(Algebra.OuterDist(A,B),1,atol=dist_tol)

end


#===========================================================================#
#
#	Lead heads and Semi-infinite GF
#
#---------------------------------------------------------------------------#


include("fakehopping.jl")


nr_leads=2

include("fakeleads.jl")









#===========================================================================#
#
#	Rule atom-layer as input or inferred
#
#---------------------------------------------------------------------------#

LeadContacts = [findall(any.(eachcol(bond(L[:head],Atoms)))) for L in userleads]



NrLayers,LayerAtomRelations = 
		Nx, Utils.FindPartners(	map(i->(i,Int64(Atoms[i,1])),axes(Atoms,1)),
													 	sortfirst=true)

q=1												

if true

NrLayers,LayerAtomRelations = 
	GreensFunctions.DistribAtoms_toLayers(Atoms, bond, LeadContacts)

q=2
end

LayerOfAtom, IndsAtomsOfLayer = LayerAtomRelations 

Graph.PlotAtoms_asGraph(Atoms,bond,fname="lattice$q",colorrule=LayerOfAtom)




#===========================================================================#
#
#	prepare hopping
#
#---------------------------------------------------------------------------#




N = NrLayers

AtomsOfLayer(l) = Atoms[IndsAtomsOfLayer(l),:]

println([size(AtomsOfLayer(l)) for l=1:NrLayers])


Hscatt(n,m=n) = HoppMatr(AtomsOfLayer(n),AtomsOfLayer(m))




#===========================================================================#
#
#	Make single left/right lead
#
#---------------------------------------------------------------------------#



VirtualLeads,RealVirtLeadRelations,LeadSizes = 
				GreensFunctions.DistribLeads_toVirtLeads(userleads,
																					 LeadContacts,
																					 Atoms,
																					 N,
																					 LayerAtomRelations)



#===========================================================================#
#
#	path graph => GF 
#
#---------------------------------------------------------------------------#



g = GreensFunctions.LayeredSystem_toGraph(Hscatt,NrLayers;VirtualLeads...)



nodelabel(i) = join(g[i,:name]," ")* (1<=i<=N ? (": ("*join(Utils.IdentifyRanges(IndsAtomsOfLayer(i))," ")*")") : "")

Graph.Plot_Graph(Graph.SimpleDiGraph(g),fname="layers$q",nodelabel=nodelabel,colorrule= i-> 1<=i<=N)






Energy = rand()


translator = GreensFunctions.AtomLead_toLayerVirtLead(
								LayerAtomRelations,d0,
								RealVirtLeadRelations,
								LeadSizes,
								)	


G = GreensFunctions.GF_Decimation_fromGraph(Energy,g,translator)


println()

#a = ("Atom",div(size(Atoms,1),4),"Atom",3*div(size(Atoms,1),4))
#for a in [ ("C",1,"Atom",22), ("A",1,"A",1)]
for a in [("Atom",1,"Atom",1)]
println(a)

G(a...) |> eachrow .|> println
println()
end



a = GreensFunctions.LDOS_Decimation(G,NrLayers,IndsAtomsOfLayer;[k=>v[:head] for (k,v) in VirtualLeads]...) 


show(stdout, "text/plain", round.(a,digits=3))






















