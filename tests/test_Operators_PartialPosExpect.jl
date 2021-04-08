#!/home/tudor/apps/julia/julia-1.1.1/bin/julia

import Operators

nr_at = 70*50

E = rand(500,nr_at)

R = rand(1:Int(ceil(sqrt(nr_at))),nr_at,2)

A = Operators.Partial_PositExpect_fromLDOS(E,R,2,convolute=true)
@time A = Operators.Partial_PositExpect_fromLDOS(E,R,2,convolute=true)

B = Operators.Partial_PositExpect_fromLDOS(E,R,2,convolute=false)
@time B = Operators.Partial_PositExpect_fromLDOS(E,R,2,convolute=false)




x1=A[1] |>collect 
y1 = A[2][1,:]
x2 = B[1] |> collect
y2 = B[2][1,:]


@show x1 x2 y1 y2

