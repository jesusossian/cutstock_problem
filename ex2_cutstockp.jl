using JuMP
using GLPK
using Gurobi
using SparseArrays

function solve_pricing(dual_demand_satisfaction, maxwidth, widths, rollcost, demand, prices)

    reduced_costs = dual_demand_satisfaction + prices
    n = length(reduced_costs)
    
    # The actual pricing model.
    submodel = Model(Gurobi.Optimizer)
    set_silent(submodel)
    @variable(submodel, xs[1:n] >= 0, Int)
    @constraint(submodel, sum(xs .* widths) <= maxwidth)
    @objective(submodel, Max, sum(xs .* reduced_costs))
    
    println("\n---------- The actual pricing problem ----------\n")
    println(submodel)
    
    optimize!(submodel)
    
    new_pattern = round.(Int, value.(xs))
    
    net_cost = rollcost - sum(new_pattern .* (dual_demand_satisfaction .+ prices))
    println("net_cost = $(net_cost)")
    
    # If the net cost of this new pattern is nonnegative, no more patterns to add.
    return net_cost >= 0 ? nothing : new_pattern

end # end function solve_pricing

function example_cutting_stock(; max_gen_cols::Int = 5_000)
    maxwidth = 100.0
    rollcost = 500.0
    
    println("maxwith = $(maxwidth), rollcost = $(rollcost)")
    
    prices = [
        167.0, 197.0, 281.0, 212.0, 225.0, 111.0, 93.0, 129.0, 108.0, 106.0,
        55.0, 85.0, 66.0, 44.0, 47.0, 15.0, 24.0, 13.0, 16.0, 14.0,
    ]
    
    widths = [
        75.0, 75.0, 75.0, 75.0, 75.0, 53.8, 53.0, 51.0, 50.2, 32.2,
        30.8, 29.8, 20.1, 16.2, 14.5, 11.0, 8.6, 8.2, 6.6, 5.1,
    ]
    
    demand = [
        38, 44, 30, 41, 36, 33, 36, 41, 35, 37,
        44, 49, 37, 36, 42, 33, 47, 35, 49, 42,
    ]
    
    nwidths = length(prices)
    n = length(widths)
    ncols = length(widths)
    
    println("nwidths = $(nwidths), n = $(n), ncols = $(ncols)")
    
    readline()
    
    # Initial set of patterns (stored in a sparse matrix: 
    # a pattern won't include many different cuts).
    patterns = SparseArrays.spzeros(UInt16, n, ncols)
    for i = 1:n
        patterns[i, i] = min(floor(Int, maxwidth / widths[i]), round(Int, demand[i]))
        tmp = patterns[i, i]
        println("pattern[$(i), $(i)] = $(tmp)")       
    end
        
    # Write the master problem with this "reduced" set of patterns.
    # Not yet integer variables: otherwise, the dual values may make no sense
    # (actually, Guobi(or GLPK) will yell at you if you're trying to get duals for
    # integer problems).
    m = Model(Gurobi.Optimizer)
    set_silent(m)
    
    @variable(m, x[1:ncols] >= 0)
    @objective(
        m, 
        Min,
        sum(
            x[p] * (rollcost - sum(patterns[j, p] * prices[j] for j = 1:n))
            for p = 1:ncols
        )
    )
    @constraint(
        m, 
        demand_satisfaction[j=1:n],
        sum(patterns[j, p] * x[p] for p = 1:ncols) >= demand[j]
    )

    println("\n---------- First the master problem(reduced) ---------- \n")    
    println(m)
    
    # First solve of the master problem.
    optimize!(m)
    
    if termination_status(m) != MOI.OPTIMAL
        warn("Master not optimal ($ncols patterns so far)")
    end
    
    # Then, generate new patterns, based on the dual information.
    while ncols - n <= max_gen_cols ## Generate at most max_gen_cols columns.
        
        println("&&&&&&&&&& new pattern: $(ncols+1) &&&&&&&&&&")
        
        readline()
                
        if ! has_duals(m)
            break
        end
        
        new_pattern = solve_pricing(
            dual.(demand_satisfaction),
            maxwidth,
            widths,
            rollcost,
            demand,
            prices,
        )
        
        # No new pattern to add to the formulation: done!
        if new_pattern === nothing
            break
        end
        
        # Otherwise, add the new pattern to the master problem, recompute the
        # duals, and go on waltzing one more time with the pricing problem.
        ncols += 1
        patterns = hcat(patterns, new_pattern)
        
        println("ncols = $(ncols)")
        
        # One new variable.
        new_var = @variable(m, [ncols], base_name = "x", lower_bound = 0)
        push!(x, new_var[ncols])
        
        # Update the objective function.
        set_objective_coefficient(
            m,
            x[ncols],
            rollcost - sum(patterns[j, ncols] * prices[j] for j = 1:n)
        )
        
        # Update the constraint number j if the new pattern impacts this production.
        for j = 1:n
            if new_pattern[j] > 0
                set_normalized_coefficient(
                    demand_satisfaction[j], new_var[ncols], new_pattern[j]
                )
            end
        end
        
        println("\n---------- The actual master problem(reduced) ---------- \n")
        println(m)
                
        # Solve the new master problem to update the dual variables.
        optimize!(m)
        
        if termination_status(m) != MOI.OPTIMAL
            @warn("Master not optimal ($ncols patterns so far)")
        end
    end
    
    # Finally, impose the master variables to be integer and resolve.
    # To be exact, at each node in the branch-and-bound tree, we would need to
    # restart the column generation process (just in case a new column would be
    # interesting to add). This way, we only get an upper bound (a feasible
    # solution).
    set_integer.(x)
    
    optimize!(m)
    
    if termination_status(m) != MOI.OPTIMAL
        @warn("Final master not optimal ($ncols patterns)")
        return
    end
    
    println("\n********** Final Solution **********\n")
    for i = 1:length(x)
        if value(x[i]) > 0.5
            println("$(round(Int, value(x[i]))) units of pattern $(i)")
        end
    end

    return

end

example_cutting_stock()
