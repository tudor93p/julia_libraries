using OrderedCollections:OrderedDict
include("../Lattices.jl")

import Algebra, Utils


println()

#Lattices.Lattice(rand(3,3), Dict("A"=>rand(3,2),"B"=>rand(3)))

#Lattices.Lattice(rand(3,3),OrderedDict("A"=>rand(3,2),"B"=>rand(3,1)))

#L = Lattices.Lattice([[sqrt(3) sqrt(3)];[1 -1]]/2, ["B"=>[0,0],"A"=>[1/3,1/3]], mode=:fractional)

L = Lattices.Lattice(sqrt(3)*1*[[cos(2π/12),sin(2π/12)] [cos(2π/12),-sin(2π/12)]], ["A"=> [0,0], "B" => [1/3,1/3]], mode=:fractional)

#a = Algebra.OuterBinary(unique(["1","1","2"]),string.(keys(L.Sublattices)),occursin)



#Lattices.PosAtoms(L)


#@show Lattices.Distances(L)

0


v = rand(2,2)


b1 = Lattices.BodyVertices_fromVectors(v, dim=1)
b2 = Lattices.BodyVertices_fromVectors(transpose(v), dim=2)

#b1 |> eachrow .|> println

println()

#b2 |> eachcol .|> println




#@show extrema(b1,dims=1)
#@show extrema(b2,dims=2)




#@show b1 .- transpose(b2)


#@show Lattices.ucs_in_UC([4,4])



@show Lattices.Superlattice(L,[1,2]) |> Lattices.AddAtoms! |> Lattices.AddAtoms!

println()

@show Lattices.Superlattice([L,L],[1,1])


@show L 
println()
println()


Lattices.AddAtoms!(L)

@show L 

println()

L = Lattices.Lattice(rand(2,2))

println()


@show L


println()


Lattices.AddAtoms!(L)


@show L

#Lattices.AddAtoms!(L,rand(2))


println()

Lattices.ShiftAtoms!(L, n=1)


@show L 

println()


println()


println()


println()


println()


































nothing



