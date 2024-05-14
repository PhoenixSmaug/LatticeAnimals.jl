using Plots
using DataStructures
using Dates
using ProgressMeter

const d = 2  # number of dimensions

mutable struct Poly
    tiles::Set{NTuple{d, Int64}}
end


"""
Return d-dimensional bounding box
"""
function boundingBox(tiles::Set{NTuple{d, Int64}})
    # Initialize vectors for min and max with extreme values
    minVec = fill(typemax(Int), d)
    maxVec = fill(typemin(Int), d)
    
    for tile in tiles
        for i in 1:d
            minVec[i] = min(minVec[i], tile[i])
            maxVec[i] = max(maxVec[i], tile[i])
        end
    end
    
    return [(minVec[i] => maxVec[i]) for i in 1:d]
end


"""
Plot 2D polyforms
"""
function prettyPrint(tiles::Set{NTuple{d, Int64}}, basis::Matrix{Float64}, p::Float64)
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    coords = [basis * [tile...] for tile in tiles]

    # Extract transformed x and y coordinates
    xCoords = [c[1] for c in coords]
    yCoords = [c[2] for c in coords]

    # Create the plot
    pl = scatter(xCoords, yCoords, title="Polyform", legend=false, xlabel="", ylabel="", aspect_ratio=:equal)

    savefig(pl, "plots/$p-$(Dates.format(now(), "HH-MM-SS-MS"))-plot.png")
end


function prettyPrint(tiles::Set{NTuple{d, Int64}}, perimeter::Set{NTuple{d, Int64}}, basis::Matrix{Float64}, p::Float64)
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    tileCoords = [basis * [tile...] for tile in tiles]
    perimeterCoords = [basis * [adj...] for adj in perimeter]

    # Extract transformed x and y coordinates for tiles
    xCoordsTiles = [c[1] for c in tileCoords]
    yCoordsTiles = [c[2] for c in tileCoords]

    # Extract transformed x and y coordinates for perimeter
    xCoordsAdjacent = [c[1] for c in perimeterCoords]
    yCoordsAdjacent = [c[2] for c in perimeterCoords]

    # Create the plot with tiles in blue and perimeter in red
    pl = scatter(xCoordsTiles, yCoordsTiles, color=:blue, label="Tiles", marker=:circle, aspect_ratio=:equal)
    scatter!(pl, xCoordsAdjacent, yCoordsAdjacent, color=:red, label="Adjacent", marker=:square)

    title!(pl, "Polyform")
    xlabel!(pl, "")
    ylabel!(pl, "")

    savefig(pl, "plots/$p-$(Dates.format(now(), "HH-MM-SS-MS"))-plot-perimeter.png")
end


"""
Breath-First Search to check connectivity
"""
@inline function bfs(tiles::Set{NTuple{d, Int64}}, start::NTuple{d, Int64}, finish::NTuple{d, Int64}, neighbours::Vector{NTuple{d, Int64}})
    q = Queue{NTuple{d, Int64}}()
    done = Set{NTuple{d, Int64}}()
    enqueue!(q, start)
    push!(done, start)

    while !isempty(q)
        tile = dequeue!(q)
        if tile == finish
            return true
        end

        for neighbour in neighbours
            neighbourTile = tile .+ neighbour

            if (neighbourTile in tiles) && !(neighbourTile in done)
                enqueue!(q, neighbourTile)
                push!(done, neighbourTile)
            end
        end
    end

    return false
end


@inline function shuffle(tiles::Set{NTuple{d, Int64}}, perimeter::Set{NTuple{d, Int64}}, p::Float64, neighbours::Vector{NTuple{d, Int64}})
    tA = length(perimeter)

    # Step 1: Uniformly at random, select a cell, x, in A

    x = rand(tiles)
    delete!(tiles, x)
    push!(perimeter, x)  # removed cell will always become part of side parameter

    delPerimeter = Set{NTuple{d, Int64}}() # tiles to be removed from side perimeter

    for neighbour in neighbours
        neighbourTile = x .+ neighbour
        # neighbours of removed tile without any other neighbours in polyform
        if !any(neighbourTile .+ n in tiles for n in neighbours) && (neighbourTile in perimeter)
            push!(delPerimeter, neighbourTile)
        end
    end

    for tile in delPerimeter
        delete!(perimeter, tile)
    end

    # Step 2: Uniformly at random, select a site, y, on the site perimeter of A \ {x}

    y = rand(perimeter)
    delete!(perimeter, y)
    push!(tiles, y)

    # tiles to be added to the side perimeter
    addPerimeter = Set{NTuple{d, Int64}}()

    for neighbour in neighbours
        neighbourTile = y .+ neighbour
        # neighbours of added tile which aren't in polyform and not already in perimeter
        if !(neighbourTile in tiles) && !(neighbourTile in perimeter)
            push!(addPerimeter, neighbourTile)
        end
    end

    for tile in addPerimeter
        push!(perimeter, tile)
    end

    # Step 3: If B = (A \ {x}) U {y} is a polyform, then accept it as the next polyform with probability min(1, (1-p)^(t_B-t_A)), where
    # t_A and t_B are the site perimeters of A and B, respectively. Otherwise, keep the current polynomino A.

    # Find x's neighboring tiles in tiles
    xNeighbours = [x .+ n for n in neighbours if (x .+ n in tiles)]

    # Check for circular connectivity among the selected neighbours
    allConnected = true
    if length(xNeighbours) > 1
        for i in 1:length(xNeighbours)
            if !bfs(tiles, xNeighbours[i], xNeighbours[i % length(xNeighbours) + 1], neighbours)
                allConnected = false
                break
            end
        end
    end

    if allConnected
        # difference in side perimeter
        diff = length(perimeter) - tA

        accepted = true  # accepted shuffle with probability (1 - p)^diff
        if !iszero(p)
            accepted = rand() < (1 - p)^diff;
        end

        if accepted
            return true
        end
    end

    # Undo all moves if not connected or not accepted
    
    # Revert to A \ {x}
    delete!(tiles, y)

    for tile in addPerimeter
        delete!(perimeter, tile)
    end
    push!(perimeter, y)

    # Revert to A
    push!(tiles, x)

    for tile in delPerimeter
        push!(perimeter, tile)
    end
    delete!(perimeter, x)

    return false
end


"""
Calculate number of holes as the number of connected components in the perimeter set minus 1. For some tesselations like the square this BFS
needs different neighbours (all 8 vertex perimeter) then for the polyform connectivity.
"""
function holes(perimeter::Set{NTuple{d, Int64}}, neighbours::Vector{NTuple{d, Int64}})
    # BFS to count connected components of empty points
    function countComponents(points)
        visited = Set{NTuple{d, Int64}}()
        components = 0
        
        while !isempty(points)
            startPoint = pop!(points)
            queue = Queue{NTuple{d, Int64}}()
            enqueue!(queue, startPoint)
            push!(visited, startPoint)

            while !isempty(queue)
                current = dequeue!(queue)
                for neighbour in neighbours
                    neighbourPoint = current .+ neighbour
                    if neighbourPoint in points && !(neighbourPoint in visited)
                        enqueue!(queue, neighbourPoint)
                        push!(visited, neighbourPoint)
                    end
                end
            end

            components += 1
            points = setdiff(points, visited)
        end

        return components
    end

    # Count the connected components in the empty points
    numberOfHoles = countComponents(copy(perimeter)) - 1
    return numberOfHoles
end


"""
Generate polyform. Neighbours should be sorted such that they are circular to improve performance of connectivity check
"""
function Poly(n::Int64, p::Float64, basis::Matrix{Float64}, neighbours::Vector{NTuple{d, Int64}})
    Poly(n, p, basis, neighbours, neighbours)
end

"""
Generate polyform and plot holes
"""
function Poly(n::Int64, p::Float64, basis::Matrix{Float64}, neighbours::Vector{NTuple{d, Int64}}, neighboursOutside::Vector{NTuple{d, Int64}})
    # Initialize tiles as line
    tiles = Set{NTuple{d, Int64}}()
    for i in 0 : n - 1
        push!(tiles, Tuple(fill(0, d)) .+ i .* neighbours[1])
    end

    # Determine initial perimeter
    perimeter = Set{NTuple{d, Int64}}()
    for tile in tiles
        for neighbour in neighbours
            neighbour_tile = tile .+ neighbour
            if !(neighbour_tile in tiles)
                push!(perimeter, neighbour_tile)
            end
        end
    end

    holeData = Vector{Int64}()

    @showprogress 1 "Shuffling..." for i in 1 : floor(Int, n^2 * log(n))
        while !shuffle(tiles, perimeter, p, neighbours)
            # Keep trying shuffle until it succeeds
        end

        if i % n == 0
            currentHoles = holes(perimeter, neighboursOutside)
            push!(holeData, currentHoles)
        end
    end

    scatter((Vector(1:length(holeData)) .* n^2), holeData, title="Development of Holes Over Time", xlabel="Iterations", legend=false)
    savefig("plots/$p-$(Dates.format(now(), "HH-MM-SS-MS"))-holes.png")

    prettyPrint(tiles, basis, p)
    #prettyPrint(tiles, perimeter, basis, p)

    println(holes(perimeter, neighboursOutside))

    open("plots/$p-$(Dates.format(now(), "HH-MM-SS-MS"))-data.txt", "w") do file
        write(file, repr(Poly(tiles)))
    end

    Poly(tiles)
end


function Polyomino(n::Int64, p::Float64)
    Poly(n, p, [1. 0.; 0. 1.], [(1, 0), (0, -1), (-1, 0), (0, 1)], [(1, 0), (0, -1), (-1, 0), (0, 1), (1, -1), (-1, -1), (-1, 1), (1, 1)])
end

function Polyhex(n::Int64, p::Float64)
    Poly(n, p, [1. -1/2; 0. sqrt(3)/2], [(1, 0), (0, -1), (-1, -1), (-1, 0), (0, 1), (1, 1)])
end


if length(ARGS) == 1
    while true
        #Polyhex(1000, parse(Float64, ARGS[1]))
        Polyomino(900, parse(Float64, ARGS[1]))
    end
end