include("../Algebra.jl")

using BenchmarkTools
using Profile
using Plots
import Utils

X = range(-2,2,length=3000);

#
#f1,f2,f3=Algebra.Heaviside, Algebra.Heaviside2,Algebra.Heaviside3;
#
#c =rand()
#
#n=100 
#
#for a in [rand(2n),rand(1)],b in [rand(n),rand(1)]
#
#	println()
#	println()
#	println()
#
#
#F1(x) = f1(a,b,c)
#
#
#F2(x) = f2(a,b,c)
#F3(x) = f3(a,b,c)
##F3(x) = f2(c)(a,b)
#
#F1(1)
#F2(1)
#F3(1)
#
#F1(2)
#F2(2)
#F3(2)
#
##@btime F1(1); @btime F2(1);@btime F3(1);
#
#println()
#
#
#
#@time  foreach(F1, 1:1000)
#@time  foreach(F2, 1:1000)
#@time  foreach(F3, 1:1000)
#
#
#@show Algebra.LA.norm(F1(1)-F2(2))
#@show Algebra.LA.norm(F1(1)-F3(2))
##Profile.print(format=:flat)
#
#
#
#
##g1(x) = (x^3+x+1)/(x^2+abs(x)+2)
##
##E  = Meta.parse("(x^3+x+1)/(x^2+abs(x)+2)")
##
##@eval g2(x) = $(E)
##
##
##
###@btime g1(rand());@btime g2(rand());
##
##
##
#
#
#
#
#
#
#
#
#
#end
#





values = rand(2)

centers = rand(3)

delta = rand()


@show sum(Algebra.Lorentzian(values, centers, delta, normalize=true),dims=2)


values = rand(2) 

centers = rand(3) 
println()

@show Algebra.normalizeDistrib(Algebra.Lorentzian(values, centers, delta).*Algebra.Gaussian(values, centers, delta))


println() 


@show Algebra.getCombinedDistrib((:Lorentzian,:Gaussian), values, centers, delta; normalize=true)


















































































































































































































































































































































































nothing
