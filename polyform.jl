using Plots
using DataStructures
using Dates
using ProgressMeter

const d = 2  # number of dimensions

struct Poly
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
function prettyPrint(tiles::Set{NTuple{d, Int64}}, basis::Matrix{Float64})
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    coords = [basis * [tile...] for tile in tiles]

    # Extract transformed x and y coordinates
    xCoords = [c[1] for c in coords]
    yCoords = [c[2] for c in coords]

    # Create the plot
    pl = scatter(xCoords, yCoords, title="Polyform", legend=false, xlabel="", ylabel="", aspect_ratio=:equal)

    timestamp = Dates.format(now(), "HH-MM-SS-MS")
    savefig(pl, "plots/plot-$(timestamp).png")
end

function prettyPrint(tiles::Set{NTuple{d, Int64}}, adjacent::Set{NTuple{d, Int64}}, basis::Matrix{Float64})
    @assert d == 2 "Only 2D-Plotting is supported"

    # Transform tile coordinates using lattice basis
    tileCoords = [basis * [tile...] for tile in tiles]
    adjacentCoords = [basis * [adj...] for adj in adjacent]

    # Extract transformed x and y coordinates for tiles
    xCoordsTiles = [c[1] for c in tileCoords]
    yCoordsTiles = [c[2] for c in tileCoords]

    # Extract transformed x and y coordinates for adjacent
    xCoordsAdjacent = [c[1] for c in adjacentCoords]
    yCoordsAdjacent = [c[2] for c in adjacentCoords]

    # Create the plot with tiles in blue and adjacent in red
    pl = scatter(xCoordsTiles, yCoordsTiles, color=:blue, label="Tiles", marker=:circle, aspect_ratio=:equal)
    scatter!(pl, xCoordsAdjacent, yCoordsAdjacent, color=:red, label="Adjacent", marker=:square)

    title!(pl, "Polyform")
    xlabel!(pl, "")
    ylabel!(pl, "")

    timestamp = Dates.format(now(), "HH-MM-SS-MS")
    savefig(pl, "plots/plot-$(timestamp).png")
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


@inline function shuffle(tiles::Set{NTuple{d, Int64}}, adjacent::Set{NTuple{d, Int64}}, p::Float64, neighbours::Vector{NTuple{d, Int64}})
    diff = 0  # circumfrence difference

    ranTile = rand(tiles)
    ranAdj = rand(adjacent)

    # Move random tile to new position
    delete!(tiles, ranTile)
    push!(tiles, ranAdj)

    # Find ranTile's neighboring tiles in tiles
    ranTileNeighbours = [ranTile .+ n for n in neighbours if (ranTile .+ n in tiles)]

    # add edges from neighbours of removed tile to circumfrence
    diff += length(ranTileNeighbours)

    # Check for circular connectivity among the selected neighbours
    allConnected = true
    if length(ranTileNeighbours) > 1
        for i in 1:length(ranTileNeighbours)
            if !bfs(tiles, ranTileNeighbours[i], ranTileNeighbours[i % length(ranTileNeighbours) + 1], neighbours)
                allConnected = false
                break
            end
        end
    end

    if allConnected
        # remove edges from neighbours of removed tile from circumfrence
        diff -= length([ranAdj .+ n for n in neighbours if (ranAdj .+ n in tiles)])

        accepted = true  # accepted shuffle with probability (1 - p)^diff
        if !iszero(p)
            accepted = rand() < (1 - p)^diff;
        end

        if accepted
            # Update adjacent set
            push!(adjacent, ranTile)

            for neighbour in neighbours
                neighbourTile = ranTile .+ neighbour
                if !any(neighbourTile .+ n in tiles for n in neighbours)
                    delete!(adjacent, neighbourTile)
                end
            end

            delete!(adjacent, ranAdj)

            for neighbour in neighbours
                neighbourTile = ranAdj .+ neighbour
                if !(neighbourTile in tiles)
                    push!(adjacent, neighbourTile)
                end
            end

            return true

        end
    end

    # Undo move if not connected or not accepted
    delete!(tiles, ranAdj)
    push!(tiles, ranTile)
    return false
end


function Poly(size::Int64, p::Float64, basis::Matrix{Float64}, neighbours::Vector{NTuple{d, Int64}})
    # Initialize tiles linearly along the first dimension from 1 to size
    tiles = Set([(i, fill(0, d-1)...) for i in 1:size])

    # Determine initial adjacent
    adjacent = Set{NTuple{d, Int64}}()
    for tile in tiles
        for neighbour in neighbours
            neighbour_tile = tile .+ neighbour
            if !(neighbour_tile in tiles)
                push!(adjacent, neighbour_tile)
            end
        end
    end

    @showprogress 1 "Shuffling..." for i in 1 : size^3
        while !shuffle(tiles, adjacent, p, neighbours)
            # Keep trying shuffle until it succeeds
        end
    end

    prettyPrint(tiles, basis)

    Poly(tiles)
end

# hexagon example
# p = Poly(30, 0.6, [1 -1/2; 0 sqrt(3)/2], [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (-1, -1)])