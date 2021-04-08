#!/home/tudor/apps/julia/julia-1.1.1/bin/julia
import LightGraphs,MetaGraphs
import Graph,Algebra,Utils,Random


const LG = LightGraphs
const MG = MetaGraphs

using TikzGraphs
using TikzPictures
import GreensFunctions,TBmodel
import LinearAlgebra;const LA=LinearAlgebra

Random.seed!(1234);

r(n,m=n) = rand(Complex{Float64},n,m)

RGF = GreensFunctions.RGF

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
																					 N,
																					 LayerAtomRelations)



#===========================================================================#
#
#	path graph => GF 
#
#---------------------------------------------------------------------------#



g = GreensFunctions.LayeredSystem_toGraph(Hscatt,NrLayers;VirtualLeads...)



nodelabel(i) = join(g[i,:name]," ")* (1<=i<=N ? (": ("*join(IndsAtomsOfLayer(i)," ")*")") : "")

Graph.Plot_Graph(Graph.SimpleDiGraph(g),fname="layers$q",nodelabel=nodelabel,colorrule= i-> 1<=i<=N )






include("oldversion.jl")


Energy =rand()

##### ------- a certain energy has been given --------- ######3


## evaluate lead GF at the given Energy; replaces props of g




H	= GreensFunctions.GraphLayeredSystem_Utils(g)[5]


node = GreensFunctions.GraphLayeredSystem_Utils(g)[end]

oldGss = oldGFDecim(Energy)

translator = GreensFunctions.AtomLead_toLayerVirtLead(
								LayerAtomRelations,d0,
								RealVirtLeadRelations,
								LeadSizes,
								)	
println()
@show translator("Atom",5)
@show translator("A",1)
@show translator("B",2,"Atom",12)

G = GreensFunctions.GF_Decimation_fromGraph(Energy,g,translator)


									

println()

#a = ("Atom",div(size(Atoms,1),4),"Atom",3*div(size(Atoms,1),4))
for a in [ ("C",1,"Atom",22), ("A",1,"A",1)]
println(a)

G(a...) |> eachrow .|> println
println()
end

leads = [(String(k),1) for (k,v) in VirtualLeads]
				 
				 

layers = [("Layer",n) for n=1:N]


attached_atoms = leads

atoms_bigH = layers


include("testbigmatrix.jl")



if true##false

#
# 	testing



function oldindex((name,index)) 

	name=="Layer" && return index

	ii = findfirst(name.==["LeftLead","RightLead"])

	!isnothing(ii) && return [1-index,N+index][ii]


end








for (dir,oldf) in collect(zip(["left","right","both"],oldGss))[3:3]
for newsrc in [leads;layers]

for newdest in layers




oldsrc = oldindex(newsrc)
olddest = oldindex(newdest)



if newsrc in leads || LA.norm(H(newsrc...,newdest...))>1e-8

println("G $dir($oldsrc,$olddest)\t $newsrc->$newdest")


new = G(newsrc...,newdest...,dir=dir)


old = oldf(oldsrc,olddest)	


#end
if isnothing(new) | isnothing(old) 
	
#	println("G $dir($oldsrc,$olddest)\t $newsrc->$newdest\t" ,isnothing(new) ?  "new is nothing " : " ",isnothing(old) ? "old is nothing " : " " ))

elseif LA.norm(old-new)>1e-8 
	
#	println("G $dir($oldsrc,$olddest)\t $newsrc->$newdest\t\t")
	println("################################# G $dir : results are different !!!! ", LA.norm(old-new) )

else
#	println("G $dir($oldsrc,$olddest)\t $newsrc->$newdest\t\t alright")
	println("\t\t\t\t\t\talright")

end
																 
end
end

end

println()

end























end

















