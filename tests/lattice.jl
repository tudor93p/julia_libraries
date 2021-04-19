using OrderedCollections:OrderedDict
include("../Lattices.jl")

import Algebra, Utils, Plots


println()

#L = Lattices.Lattice(sqrt(3)*1*[[cos(2π/12),sin(2π/12)] [cos(2π/12),-sin(2π/12)]], ["A"=> [0,0], "B" => [1/3,1/3]], mode=:fractional)


println()



L = Lattices.Lattice(Utils.UnitMatrix(2))

Lattices.AddAtoms!(L)

Lattices.Superlattice!(L, [5,3])





#@show L


#Lattices.PosAtoms(Lattices.Superlattice(Lattices.Lattice(Utils.UnitMatrix(2))|>Lattices.AddAtoms, [10,10])).-Lattices.PosAtoms(L) .|> abs|> maximum |> println


#Lattices.ReduceDim!(L)


#@show L

L1 = L



L2 = Utils.UnitMatrix(2) |> Lattices.Lattice |> Lattices.AddAtoms 

Lattices.ShiftAtoms!(L2, n=[0,1.6])

Lattices.Superlattice!(L2, (2,5))


Lattices.ReduceDim!(L2, 2)




Lattices.Align_toAtoms(L2, Lattices.PosAtoms(L1))

Lattices.plot(L2,L1)























