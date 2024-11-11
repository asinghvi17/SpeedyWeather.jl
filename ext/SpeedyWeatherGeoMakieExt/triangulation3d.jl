using GeoMakie.Makie.LinearAlgebra
using GeoMakie.GeometryBasics
using GeoMakie.GO.ExactPredicates
using GeoMakie.GO.ExactPredicates.Codegen
using GeoMakie.GO.ExactPredicates.StaticArrays: SVector

# This function is an "exact predicate" that is robust to floating point error
# May not be necessary, but it's here if we need it.
@genpredicate function _collinear(p :: 3, q :: 3, r :: 3)
    pq = q - p
    pr = r - p
    Codegen.group!(pq...)
    Codegen.group!(pr...)
    # Cross product of vectors will be zero if points are collinear
    cross = SVector{3}(
        pq[2]*pr[3] - pq[3]*pr[2],
        pq[3]*pr[1] - pq[1]*pr[3], 
        pq[1]*pr[2] - pq[2]*pr[1]
    )
    return ExactPredicates.inp(cross, cross) # Will be 0 if collinear, positive otherwise
end

function GeometryBasics.earcut_triangulate(polygon::Vector{<: Vector{<: Makie.VecTypes{3,T}}}) where T
    # Here, we are assuming that the polygon is actually planar,
    # but the plane is in 3D, and not necessarily the XY plane.
    # So, we can find the plane of best fit using the first three points of the polygon!
    p1, p2, p3 = polygon[1][1], polygon[1][2], polygon[1][3]
    # Account for equal points
    if p1 == p2 || p1 == p3 || p2 == p3
        if length(polygon[1]) <= 3
            error("Polygon has only three points and they are all the same, we can't triangulate this!")
        elseif p1 == p2 == p3
            new_point_idx = findfirst(p -> p != p1, polygon[1])
            if new_point === nothing
                error("All points in the polygon are the same, we can't triangulate this!")
            end
            p2 = polygon[1][new_point_idx]
            new_point_idx = findfirst(p -> p != p1 && p != p2, polygon[1])
            if new_point_idx === nothing
                error("Only found two unique points in the polygon, we can't triangulate this!")
            end
            p3 = polygon[1][new_point_idx]
        elseif p1 == p2
            p2 = polygon[1][4]
        elseif p1 == p3
            p3 = polygon[1][4]
        elseif p2 == p3
            p2 = polygon[1][4]
        end
    end

    # Account for collinear points
    if _collinear(p1, p2, p3) == 0 # collinear, all the points lie on the same line
        if length(polygon[1]) <= 3
            error("Polygon has only three points and they are all collinear, we can't triangulate this!")
        end
        new_point_idx = findfirst(p -> _collinear(p1, p2, p) != 0, polygon[1])
        if new_point_idx === nothing
            error("All points in the polygon are collinear, we can't triangulate this!")
        end
        p3 = polygon[1][new_point_idx]
    end

    # Define a plane that can be used to project the polygon into 2D
    v1 = p2 - p1
    v2 = p3 - p1
    normal = cross(v1, v2)
    x = v1
    y = cross(normal, x)

    # Project the polygon into 2D
    projected_polygon = map(ring -> map(p -> Point2{Float64}(dot(p, x), dot(p, y)), ring), polygon)

    # Now, call earcut_triangulate on the projected polygon, which is 2D
    return GeometryBasics.earcut_triangulate(projected_polygon)
end