import Utils 


fn(x) = string(x)

Write!,d = Utils.Write_NamesVals(fn)

@show d

Write!("testobs", rand(10,3), d)

@show d
Write!("testobs2", rand(10,3), d)

@show d
