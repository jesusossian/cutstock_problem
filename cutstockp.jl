using JuMP
using LinearAlgebra
using Gurobi

W = 100
println(W)

tamM = 5
println(tamM)

M = 1:tamM
println(M)

A = Array{Int64}(I,tamM,tamM)
#A = zeros(tamM,tamM)
#A = Diagonal(ones(tamM))
println(A)

p = zeros(5)
println(p)

b = [45; 38; 25; 11; 12]
println(b)

w = [22; 42; 52; 53; 78]
println(w)

master = Model(solver = GurobiSolver()) # modelo for the master problem

Jp = 1:size(A,2) # initial number of variables
println(Jp)

@variable(master, 0 <= x[Jp] <= Inf) # setting variables
@objective(master, Min, sum(1*x[j] for j in Jp)) # setting objective function
@constraint(master, rest[i=1:tamM], sum(A[i,j] * x[j] for  j in Jp) == b[i])

println(master)

# Solving the master problem with the initial data
solve(master)
println("Current solution of the master problem is ", getvalue(x))
println("Current objective value of the master problem is ", getobjectivevalue(master))

#Collect the dual variables for the equality constraints and store them in an array p
for i in M
  p[i] = getdual(rest[i]) # These p[i] are the input data for the subproblem
end 
println("The array storing the dual variables is ", p)

# Describe the sub problem
# ------------------------
subprob = Model(solver = GurobiSolver()) # Model for the subproblem
@variable(subprob, 0 <= y[M] <= Inf, Int )
@objective(subprob, Min, 1 - sum(p[i] * y[i] for i in M) )
@constraint(subprob, sum(w[i] * y[i] for i in M) <= W)
print(subprob)

# Solve the sub problem
# ---------------------
solve(subprob)
minreducedCost = getobjectivevalue(subprob)
println("The minimum component of the reduced cost vector is ", minreducedCost)

if minreducedCost >= 0
    println("We are done, current solution of the master problem is optimal")
else
    println("We have a cost reducing column ", getvalue(y))
end

println(typeof(y))

#Anew=Float64[] # This Anew correspond to the newly added column to the A matrix
#for (i in 1:tamM)
#    push!(Anew, getvalue(y)[i])
#end
