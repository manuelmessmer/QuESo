/*
  ____        ______  _____
 / __ \      |  ____|/ ____|
| |  | |_   _| |__  | (___   ___
| |  | | | | |  __|  \___ \ / _ \
| |__| | |_| | |____ ____) | (_) |
 \___\_\\__,_|______|_____/ \___/
        Quadrature for Embedded Solids

 License:    BSD 4-Clause License
             See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE

 Authors:    Manuel Messmer
*/

//// STL includes
#include <algorithm>
#include <map>

//// Project includes
#include "queso/utilities/mesh_utilities.h"
#include "queso/utilities/math_utilities.hpp"
#include "queso/utilities/triangle_utilities.hpp"

namespace queso {
namespace MeshUtilities {

void Refine(TriangleMesh &rTriangleMesh, IndexType MinNumberOfTriangles)
{
    while (rTriangleMesh.NumOfTriangles() < MinNumberOfTriangles) {
        double max_area = 0.0;
        TriangleMeshView(rTriangleMesh).VisitEachTriangle<WithoutNormals>([&](const auto &triangle) {
            max_area = std::max(max_area, TriangleUtilities::Area(triangle));
        });

        if (max_area <= ZEROTOL) {
            break;
        }

        const double split_threshold = 0.9 * max_area;
        const auto old_vertices = rTriangleMesh.Vertices();
        const auto old_triangles = rTriangleMesh.TriangleIndices();
        const auto old_normals = rTriangleMesh.Normals();

        TriangleMesh refined;
        refined.Reserve(old_vertices.size() + 3 * old_triangles.size());
        for (const auto &v : old_vertices) {
            refined.AddVertex(v);
        }

        std::map<std::pair<IndexType, IndexType>, IndexType> midpoint_cache;
        auto midpoint_id = [&](IndexType i, IndexType j) {
            if (j < i) {
                std::swap(i, j);
            }
            const auto key = std::make_pair(i, j);
            const auto it = midpoint_cache.find(key);
            if (it != midpoint_cache.end()) {
                return it->second;
            }

            const auto midpoint = (0.5 * (old_vertices[i] + old_vertices[j]));
            const IndexType id = refined.AddVertex(midpoint);
            midpoint_cache.emplace(key, id);
            return id;
        };

        bool did_split = false;
        for (IndexType tid = 0; tid < old_triangles.size(); ++tid) {
            const auto tri = old_triangles[tid];
            const auto proxy = rTriangleMesh.Triangle<WithoutNormals>(tid);
            const double area = TriangleUtilities::Area(proxy);
            const auto &parent_normal = old_normals[tid];

            if (area > split_threshold) {
                did_split = true;
                const IndexType e1 = midpoint_id(tri[0], tri[1]);
                const IndexType e2 = midpoint_id(tri[1], tri[2]);
                const IndexType e3 = midpoint_id(tri[2], tri[0]);

                refined.AddTriangle({ tri[0], e1, e3 }, parent_normal);
                refined.AddTriangle({ e1, tri[1], e2 }, parent_normal);
                refined.AddTriangle({ e2, tri[2], e3 }, parent_normal);
                refined.AddTriangle({ e1, e2, e3 }, parent_normal);
            } else {
                refined.AddTriangle(tri, parent_normal);
            }
        }

        if (!did_split) {
            break;
        }

        rTriangleMesh = std::move(refined);
    }
}

void Append(TriangleMesh &rTriangleMesh, const TriangleMesh &rNewMesh)
{
    const IndexType initial_vertex_count = rTriangleMesh.NumOfVertices();
    const IndexType target_capacity = rTriangleMesh.NumOfTriangles() + rNewMesh.NumOfTriangles();
    rTriangleMesh.Reserve(target_capacity);

    for (const auto &vertex : rNewMesh.Vertices()) {
        rTriangleMesh.AddVertex(vertex);
    }

    const auto source_triangles = rNewMesh.TriangleIndices();
    const auto source_normals = rNewMesh.Normals();
    for (IndexType tid = 0; tid < source_triangles.size(); ++tid) {
        const auto tri = source_triangles[tid];
        const Vector3i offset_triangle{
            tri[0] + initial_vertex_count,
            tri[1] + initial_vertex_count,
            tri[2] + initial_vertex_count
        };
        rTriangleMesh.AddTriangle(offset_triangle, source_normals[tid]);
    }
}

// void Append(TriangleMesh &rTriangleMesh, const TriangleMesh &rNewMesh, const std::vector<IndexType> &rIndices)
// {
//     const auto source_vertices = rNewMesh.Vertices();
//     const auto source_triangles = rNewMesh.TriangleIndices();
//     const auto source_normals = rNewMesh.Normals();
//
//     const IndexType triangle_capacity = rTriangleMesh.NumOfTriangles() + rIndices.size();
//     const IndexType vertex_capacity = rTriangleMesh.NumOfVertices() + 3 * rIndices.size();
//     rTriangleMesh.Reserve(std::max(vertex_capacity, triangle_capacity));
//
//     constexpr IndexType invalid_id = std::numeric_limits<IndexType>::max();
//     std::vector<IndexType> remapped_vertex_ids(source_vertices.size(), invalid_id);
//
//     for (const IndexType tid : rIndices) {
//         const auto tri = source_triangles[tid];
//         Vector3i remapped{};
//
//         for (IndexType j = 0; j < 3; ++j) {
//             const IndexType source_vid = tri[j];
//             IndexType mapped_id = remapped_vertex_ids[source_vid];
//             if (mapped_id == invalid_id) {
//                 mapped_id = rTriangleMesh.AddVertex(source_vertices[source_vid]);
//                 remapped_vertex_ids[source_vid] = mapped_id;
//             }
//             remapped[j] = mapped_id;
//         }
//
//         rTriangleMesh.AddTriangle(remapped, source_normals[tid]);
//     }
// }

TriangleMesh MakeMeshBox(const PointType &rLowerPoint, const PointType &rUpperPoint)
{
	TriangleMesh new_triangle_mesh{};

	new_triangle_mesh.AddVertex({ rLowerPoint[0], rLowerPoint[1], rLowerPoint[2] }); // 0
	new_triangle_mesh.AddVertex({ rUpperPoint[0], rLowerPoint[1], rLowerPoint[2] }); // 1
	new_triangle_mesh.AddVertex({ rLowerPoint[0], rUpperPoint[1], rLowerPoint[2] }); // 2
	new_triangle_mesh.AddVertex({ rUpperPoint[0], rUpperPoint[1], rLowerPoint[2] }); // 3
	new_triangle_mesh.AddVertex({ rLowerPoint[0], rLowerPoint[1], rUpperPoint[2] }); // 4
	new_triangle_mesh.AddVertex({ rUpperPoint[0], rLowerPoint[1], rUpperPoint[2] }); // 5
	new_triangle_mesh.AddVertex({ rLowerPoint[0], rUpperPoint[1], rUpperPoint[2] }); // 6
	new_triangle_mesh.AddVertex({ rUpperPoint[0], rUpperPoint[1], rUpperPoint[2] }); // 7

	new_triangle_mesh.AddTriangle({ 0, 6, 2 }, { -1.0, 0.0, 0.0 });
	new_triangle_mesh.AddTriangle({ 0, 4, 6 }, { -1.0, 0.0, 0.0 });

	new_triangle_mesh.AddTriangle({ 1, 7, 5 }, { 1.0, 0.0, 0.0 });
	new_triangle_mesh.AddTriangle({ 1, 3, 7 }, { 1.0, 0.0, 0.0 });

	new_triangle_mesh.AddTriangle({ 4, 1, 5 }, { 0.0, -1.0, 0.0 });
	new_triangle_mesh.AddTriangle({ 4, 0, 1 }, { 0.0, -1.0, 0.0 });

	new_triangle_mesh.AddTriangle({ 6, 7, 3 }, { 0.0, 1.0, 0.0 });
	new_triangle_mesh.AddTriangle({ 6, 3, 2 }, { 0.0, 1.0, 0.0 });

	new_triangle_mesh.AddTriangle({ 1, 0, 3 }, { 0.0, 0.0, -1.0 });
	new_triangle_mesh.AddTriangle({ 0, 2, 3 }, { 0.0, 0.0, -1.0 });

	new_triangle_mesh.AddTriangle({ 4, 5, 7 }, { 0.0, 0.0, 1.0 });
	new_triangle_mesh.AddTriangle({ 4, 7, 6 }, { 0.0, 0.0, 1.0 });

	return new_triangle_mesh;
}

double Area(const TriangleMeshView &rTriangleMeshView)
{
    double area = 0.0;
    rTriangleMeshView.VisitEachTriangle<WithoutNormals>([&](const auto &triangle) {
        area += TriangleUtilities::Area(triangle);
    });
    return area;
}

double AreaOMP(const TriangleMeshView &rTriangleMeshView)
{
    double area = 0.0;
    const IndexType num_triangles = rTriangleMeshView.NumOfTriangles();
    #pragma omp parallel for reduction(+ : area)
    for (int i = 0; i < static_cast<int>(num_triangles); ++i) {
        area += TriangleUtilities::Area(rTriangleMeshView.Triangle<WithoutNormals>(static_cast<IndexType>(i)));
    }
    return area;
}

double Volume(const TriangleMeshView &rTriangleMeshView)
{
    double volume = 0.0;
    rTriangleMeshView.VisitEachTriangle<WithNormals>([&](const auto &triangle) {
        const auto center = TriangleUtilities::Center(triangle);
        volume += Math::Dot(triangle.Normal, center) * TriangleUtilities::Area(triangle);
    });
    return std::abs(volume / 3.0);
}

double VolumeOMP(const TriangleMeshView &rTriangleMeshView)
{
    double volume = 0.0;
    const IndexType num_triangles = rTriangleMeshView.NumOfTriangles();
    #pragma omp parallel for reduction(+ : volume)
    for (int i = 0; i < static_cast<int>(num_triangles); ++i) {
        const auto triangle = rTriangleMeshView.Triangle<WithNormals>(static_cast<IndexType>(i));
        const auto center = TriangleUtilities::Center(triangle);
        volume += Math::Dot(triangle.Normal, center) * TriangleUtilities::Area(triangle);
    }
    return std::abs(volume / 3.0);
}

double Volume(const TriangleMeshView &rTriangleMeshView, IndexType Dir)
{
    QuESo_ERROR_IF(Dir > 2UL) << "Directional Index is out-of-range.\n";

    double volume = 0.0;
    rTriangleMeshView.VisitEachTriangle<WithNormals>([&](const auto &triangle) {
        const auto center = TriangleUtilities::Center(triangle);
        volume += triangle.Normal[Dir] * center[Dir] * TriangleUtilities::Area(triangle);
    });
    return std::abs(volume);
}

double MaxAspectRatio(const TriangleMeshView &rTriangleMeshView)
{
    double max_aspect_ratio = MIND;
    rTriangleMeshView.VisitEachTriangle<WithoutNormals>([&](const auto &triangle) {
        max_aspect_ratio = std::max(max_aspect_ratio, TriangleUtilities::AspectRatio(triangle));
    });
    return max_aspect_ratio;
}

double AverageAspectRatio(const TriangleMeshView &rTriangleMeshView)
{
    if (rTriangleMeshView.NumOfTriangles() == 0) {
        return 0.0;
    }

    double sum = 0.0;
    rTriangleMeshView.VisitEachTriangle<WithoutNormals>([&](const auto &triangle) {
        sum += TriangleUtilities::AspectRatio(triangle);
    });
    return sum / static_cast<double>(rTriangleMeshView.NumOfTriangles());
}

double EstimateQuality(const TriangleMeshView &rTriangleMeshView)
{
    double total_volume_1 = 0.0;
    double total_volume_2 = 0.0;
    double total_volume_3 = 0.0;
    double total_area = 0.0;
    PointType directional_areas = { 0.0, 0.0, 0.0 };

    rTriangleMeshView.VisitEachTriangle<WithNormals>([&](const auto &triangle) {
        const double area = TriangleUtilities::Area(triangle);
        const auto center = TriangleUtilities::Center(triangle);

        total_area += area;
        directional_areas[0] += area * triangle.Normal[0];
        directional_areas[1] += area * triangle.Normal[1];
        directional_areas[2] += area * triangle.Normal[2];

        total_volume_1 += triangle.Normal[0] * center[0] * area;
        total_volume_2 += triangle.Normal[1] * center[1] * area;
        total_volume_3 += triangle.Normal[2] * center[2] * area;
    });

    const double total_volume = std::abs((total_volume_1 + total_volume_2 + total_volume_3) / 3.0);
    const double denom = std::max(total_volume, EPS4);
    const double error_v1 = std::abs(total_volume_1 - total_volume) / denom;
    const double error_v2 = std::abs(total_volume_2 - total_volume) / denom;
    const double error_v3 = std::abs(total_volume_3 - total_volume) / denom;

    const double error_area = Math::Norm(directional_areas) / std::max(std::abs(total_area), EPS4);
    return std::max(std::max(std::max(error_v1, error_v2), error_v3), error_area);
}

std::pair<PointType, PointType> BoundingBox(const TriangleMeshView &rTriangleMeshView)
{
    PointType lower_bound = { MAXD, MAXD, MAXD };
    PointType upper_bound = { LOWESTD, LOWESTD, LOWESTD };

    rTriangleMeshView.VisitEachTriangle<WithoutNormals>([&](const auto &triangle) {
        const PointType x_values = { triangle.P1[0], triangle.P2[0], triangle.P3[0] };
        const PointType y_values = { triangle.P1[1], triangle.P2[1], triangle.P3[1] };
        const PointType z_values = { triangle.P1[2], triangle.P2[2], triangle.P3[2] };

        const auto x_min_max = std::minmax_element(x_values.begin(), x_values.end());
        const auto y_min_max = std::minmax_element(y_values.begin(), y_values.end());
        const auto z_min_max = std::minmax_element(z_values.begin(), z_values.end());

        lower_bound[0] = std::min<double>(*x_min_max.first, lower_bound[0]);
        upper_bound[0] = std::max<double>(*x_min_max.second, upper_bound[0]);

        lower_bound[1] = std::min<double>(*y_min_max.first, lower_bound[1]);
        upper_bound[1] = std::max<double>(*y_min_max.second, upper_bound[1]);

        lower_bound[2] = std::min<double>(*z_min_max.first, lower_bound[2]);
        upper_bound[2] = std::max<double>(*z_min_max.second, upper_bound[2]);
    });

    return std::make_pair(lower_bound, upper_bound);
}

} // namespace MeshUtilities
} // namespace queso



