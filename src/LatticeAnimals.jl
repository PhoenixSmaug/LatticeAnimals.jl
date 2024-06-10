module LatticeAnimals

using Plots
using DataStructures
using ProgressMeter
using Preferences

export setDimension, Poly, Polyomino, Polyhex, polyPlot


"""
    setDimension(d)

Set the number of dimensions for the polyforms as a disk-persistent setting. A rebuilt of the package is necessary afterwards

# Arguments
    * `d`: Number of dimensions
"""
function setDimension(d::Int64)
    if !(d >= 2)
        throw(ArgumentError("Invalid dimension: $d"))
    end

    # Set it in our runtime values, as well as saving it to disk
    @set_preferences!("dimension" => d)
    @info("New dimension set; restart your Julia session for this change to take effect!")
end

const d = @load_preference("dimension", 2)


"""
    Poly(tiles, perimeter, holes)

Struct to represent a d-dimensional polyform (Animal lattice).

# Arguments
    * `tiles`: Lattice points which are part of the polyform
    * `perimeter`: Lattice points which are in the side perimeter of the polyform, so are neighboring a tile in the polyform
    * `holes`: Number of d-dimensional holes
"""
mutable struct Poly
    tiles::Set{NTuple{d, Int64}}
    perimeter::Set{NTuple{d, Int64}}
    holes::Int64
end


"""
    Poly(n, p, basis, neighbours; doPlot)

Generate a random d-dimensional polyform of size n from the percolation distribution at p using the Metropolis algorithm. neighbours
is a list of coordinate differences to each respective neighbouring lattice point defined by the lattice and should be sorted
clock-wise to improve performance of the connectivity checks during generation. The default mixing time is 2*n^2

# Arguments
    * `n`: Size of the polyform
    * `p`: Percolation factor in [0, 1)
    * `basis`: Lattice basis of the polyform
    * `neighbours`: List of difference to the respective neighbouring tiles for the current polyform 
    * `doPlot`: If a scatter plot visualizing the polyform should be created (only available for 2d-polyforms)
"""
function Poly(n::Int64, p::Float64, basis::Matrix{Float64}, neighbours::Vector{NTuple{d, Int64}}; doPlot = false)
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

    @showprogress 1 "Shuffling..." for _ in 1 : floor(Int, 2 * n^2)
        while !shuffle(tiles, perimeter, p, neighbours)
            # Keep trying shuffle until it succeeds
        end
    end

    if doPlot
        polyPlot(tiles, basis, "$p-$(Dates.format(now(), "HH-MM-SS-MS"))-plot.png")
        #polyPlot(tiles, perimeter, basis, "$p-$(Dates.format(now(), "HH-MM-SS-MS"))-plot.png")
    end

    Poly(tiles, perimeter, holes(tiles, neighbours))
end


"""
    Polyomino(n, p)

Generate a random polyomino of size n from the percolation distribution at p using the Metropolis algorithm. The basis consits
of the two unit-vectors as columns in a matrix. The four neighbours are then specified each as a linear combination of the basis
vectors.

# Arguments
    * `n`: Size of the polyomino
    * `p`: Percolation factor in [0, 1)
"""
function Polyomino(n::Int64, p::Float64)
    Poly(n, p, [1. 0.; 0. 1.], [(1, 0), (0, -1), (-1, 0), (0, 1)])
end


"""
    Polyhex(n, p)

Generate a random polyhex of size n from the percolation distribution at p using the Metropolis algorithm. The basis consits
of the first unit vector and the vector of length 1 at 120 degrees to the left. The six neighbours are then specified each as a
linear combination of the basis vectors.

# Arguments
    * `n`: Size of the polyhex
    * `p`: Percolation factor in [0, 1)
"""
function Polyhex(n::Int64, p::Float64)
    Poly(n, p, [1. -1/2; 0. sqrt(3)/2], [(1, 0), (0, -1), (-1, -1), (-1, 0), (0, 1), (1, 1)])
end


"""
    boundingBox(tiles)

Calculate the smallest d-dimensional bounding box enclosing the polyform and return it as (min, max) pair for each dimension.

# Arguments
    * `tiles`: Lattice points in polyform
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
    polyPlot(tiles, basis, p)

Represent 2-dimensional polyform with a scatter plot.

# Arguments
    * `tiles`: Lattice points in polyform
    * `basis`: Lattice basis of the polyform
    * `path`: Output path
"""
function polyPlot(tiles::Set{NTuple{d, Int64}}, basis::Matrix{Float64}, path::String)
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    coords = [basis * [tile...] for tile in tiles]

    # Extract transformed x and y coordinates
    xCoords = [c[1] for c in coords]
    yCoords = [c[2] for c in coords]

    # Create the plot
    pl = scatter(xCoords, yCoords, title="Polyform", legend=false, xlabel="", ylabel="", aspect_ratio=:equal)

    savefig(pl, path)
end


"""
    polyPlot(tiles, basis, p)

Represent 2-dimensional polyform (blue) and its side perimeter (red) with a scatter plot.

# Arguments
    * `tiles`: Lattice points in polyform
    * `perimeter`: Lattice points in side perimeter
    * `basis`: Lattice basis of the polyform
    * `path`: Output path
"""
function polyPlot(tiles::Set{NTuple{d, Int64}}, perimeter::Set{NTuple{d, Int64}}, basis::Matrix{Float64}, path::String)
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    tileCoords = [basis * [tile...] for tile in tiles]
    perimeterCoords = [basis * [adj...] for adj in perimeter]

    # Extract transformed x and y coordinates for tiles
    xCoordsTiles = [c[1] for c in tileCoords]
    yCoordsTiles = [c[2] for c in tileCoords]

    # Extract transformed x and y coordinates for perimeter
    xCoordsPeri = [c[1] for c in perimeterCoords]
    yCoordsPeri = [c[2] for c in perimeterCoords]

    # Create the plot with tiles in blue and perimeter in red
    pl = scatter(xCoordsTiles, yCoordsTiles, color=:blue, label="Tiles", marker=:circle, aspect_ratio=:equal)
    scatter!(pl, xCoordsPeri, yCoordsPeri, color=:red, label="Perimeter", marker=:square)

    title!(pl, "Polyform")
    xlabel!(pl, "")
    ylabel!(pl, "")

    savefig(pl, path)
end


"""
    bfs(tiles, start, finish, neighbours)

Determine if a path in the polyform exists connecting the lattice point start and the lattice point end using Breath-First Search.

# Arguments
    * `tiles`: Lattice points in polyform
    * `start`: Starting point of search
    * `finish`: Ending point of search
    * `neighbours`: List of difference to the respective neighbouring tiles for the current polyform 
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

        for n in neighbours
            nTile = tile .+ n

            if (nTile in tiles) && !(nTile in done)
                enqueue!(q, nTile)
                push!(done, nTile)
            end
        end
    end

    return false
end


"""
    shuffle(tiles, perimeter, p, neighbours)

Execute one Markov step of the Metropolis algorithm:
1) Uniformly at random, select a tile, x, in the polyform A
2) Uniformly at random, select a tile, y, on the site perimeter of A \\setminus {x}
3) If B = (A \\setminus {x}) U {y} is a still connected, then accept it as the next polyform with probability min(1, (1-p)^(t_B-t_A)), where
   t_A and t_B are the site perimeters of A and B, respectively. Otherwise, keep the current polyform A.

# Arguments
    * `tiles`: Lattice points in polyform
    * `perimeter`: Lattice points in side perimeter
    * `p`: Percolation factor in [0, 1)
    * `neighbours`: List of difference to the respective neighbouring tiles for the current polyform 
"""
@inline function shuffle(tiles::Set{NTuple{d, Int64}}, perimeter::Set{NTuple{d, Int64}}, p::Float64, neighbours::Vector{NTuple{d, Int64}})
    tA = length(perimeter)

    # Step 1: Uniformly at random, select a tile, x, in the polyform A

    x = rand(tiles)
    delete!(tiles, x)
    push!(perimeter, x)  # removed cell will always become part of side parameter

    delPerimeter = Vector{NTuple{d, Int64}}() # tiles to be removed from side perimeter

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

    # Step 2: Uniformly at random, select a tile, y, on the site perimeter of A \ {x}

    y = rand(perimeter)
    delete!(perimeter, y)
    push!(tiles, y)

    # tiles to be added to the side perimeter
    addPerimeter = Vector{NTuple{d, Int64}}()

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
    holes(tiles, neighbours)

Count the number of d-dimensional holes of the polyform. This is done by determining the number of connected components of the complement, so
all lattice points not in the polyform.

# Arguments
    * `tiles`: Lattice points in polyform
    * `neighbours`: List of difference to the respective neighbouring tiles for the current polyform 
"""
function holes(tiles::Set{NTuple{d, Int64}}, neighbours::Vector{NTuple{d, Int64}})
    bounds = boundingBox(tiles)
    m = length(bounds)

    # Adjust bounds to include the boundary of the bounding box
    expandedBounds = [(b.first - 1 => b.second + 1) for b in bounds]

    # Generate all points within the expanded bounding box
    ranges = [b.first : b.second for b in expandedBounds]
    gridPoints = Set{NTuple{m, Int64}}(Iterators.product(ranges...))

    # Difference set to find the points not occupied by tiles in the bounding box
    emptyPoints = setdiff(gridPoints, tiles)

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
    numberOfHoles = countComponents(emptyPoints) - 1
    return numberOfHoles
end

# (c) Mia Muessig


end # module LatticeAnimals
