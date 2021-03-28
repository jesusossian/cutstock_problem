using JuMP
using LinearAlgebra

W = 10
println(W)

tamM = 5
println(tamM)

M = [1:tamM]
println(M)

#A = Array{Int64}(I,tamM,tamM)
#A = zeros(tamM,tamM)
A = Diagonal(ones(tamM))
println(A)

p = zeros(5)
println(p)

b = [45; 38; 25; 11; 12]
println(b)

w = [22; 42; 52; 53; 78]
println(w)

master = Model() # modelo for the master problem

Jp = [1:size(A,2)] # initial number of variables
println(Jp)

@variable(master, 0 <= x[Jp] <= Inf) # setting variables
@objective(master, Min, sum(x[j] for j in Jp)) # setting objective function
@constraint(master, rest[i=1:tamM], sum( A[i,j] for j in Jp) == b[i])
#@constraint(master, rest[i=1:tamM], sum(A[i,j]*x[j] for j in Jp) == b[i])

println(master)
