//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

//// STL includes
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stack>

//// Project includes
#include "queso/embedding/domain_topology.hpp"
#include "queso/embedding/flood_fill.h"

namespace queso::embedding {
namespace {

    using Direction = GridIndexer::Direction;
    using PartitionBox = PartitionBoxType;
    using PartitionBoxes = std::vector<PartitionBox>;

    /// @brief Tracks whether a fillable cell has been assigned to a partition-local component.
    enum class VisitState : std::uint8_t {
        unvisited,  ///< Cell remains available as a component seed.
        visited  ///< Cell is trimmed or belongs to one discovered component.
    };

    /// @brief Connected non-trimmed cells discovered globally or within one partition.
    /// @details Partition-local fragments retain signed classification evidence and physical exterior contact so these
    ///          properties can be combined when interfaces are reconciled.
    struct Component
    {
        std::vector<IndexType> cell_indices;  ///< Zero-based cells belonging to this component or fragment.
        std::int64_t classification_vote{};  ///< Signed sum of blocked-face classification evidence.
        bool touches_grid_exterior{};  ///< Whether any member cell reaches the physical grid exterior.
    };

    using Components = std::vector<Component>;

    /// @brief Union-find structure used to merge partition-local component fragments.
    class DisjointSet
    {
    public:
        /// @brief Creates Size singleton sets.
        explicit DisjointSet(IndexType Size) : mParents(Size), mRanks(Size, 0)
        {
            for (IndexType index = 0; index < Size; ++index) { mParents[index] = index; }
        }

        /// @brief Returns the representative of one set and compresses the traversed path.
        [[nodiscard]] IndexType Find(IndexType Index)
        {
            IndexType root = Index;
            while (mParents[root] != root) { root = mParents[root]; }
            while (mParents[Index] != Index) {
                const IndexType parent = mParents[Index];
                mParents[Index] = root;
                Index = parent;
            }
            return root;
        }

        /// @brief Merges two sets using rank to keep the parent trees shallow.
        void Union(IndexType First, IndexType Second)
        {
            IndexType first_root = Find(First);
            IndexType second_root = Find(Second);
            if (first_root == second_root) { return; }
            if (mRanks[first_root] < mRanks[second_root]) { std::swap(first_root, second_root); }
            mParents[second_root] = first_root;
            if (mRanks[first_root] == mRanks[second_root]) { ++mRanks[first_root]; }
        }

    private:
        std::vector<IndexType> mParents;  ///< Parent index for each component.
        std::vector<std::uint8_t> mRanks;  ///< Upper bound on each root tree's height.
    };

    /// @brief Selects the longest grid axis for slab partitioning.
    /// @details Equal axis sizes retain the first axis in x, y, z order.
    [[nodiscard]] IndexType PartitionAxis(const Vector3i& rNumberOfElements) noexcept
    {
        IndexType axis = 0;
        for (IndexType i = 1; i < 3; ++i) {
            if (rNumberOfElements[i] > rNumberOfElements[axis]) { axis = i; }
        }
        return axis;
    }

    /// @brief Resolves the number of flood-fill slabs.
    /// @details Explicit options are validated against the selected axis. Automatic selection uses the available
    ///          OpenMP thread count when enabled and one partition otherwise.
    [[nodiscard]] IndexType PartitionCount(IndexType LargestAxisSize, const FloodFillOptions& rOptions)
    {
        if (rOptions.partition_count.has_value()) {
            QuESo_ERROR_IF(*rOptions.partition_count == 0) << "Flood-fill partition count must be positive.\n";
            QuESo_ERROR_IF(*rOptions.partition_count > LargestAxisSize)
                << "Flood-fill partition count cannot exceed the largest grid direction.\n";
            return *rOptions.partition_count;
        }
#ifdef _OPENMP
        return std::min(LargestAxisSize, static_cast<IndexType>(omp_get_max_threads()));
#else
        static_cast<void>(LargestAxisSize);
        return 1;
#endif
    }

    /// @brief Divides the grid into contiguous, exhaustive, non-overlapping slabs.
    /// @details Slab sizes differ by at most one cell, with larger slabs ordered first.
    [[nodiscard]] PartitionBoxes
        CreatePartitions(const Vector3i& rNumberOfElements, IndexType PartitionAxisIndex, IndexType NumberOfPartitions)
    {
        PartitionBoxes partitions;
        partitions.reserve(NumberOfPartitions);
        const IndexType axis_size = rNumberOfElements[PartitionAxisIndex];
        const IndexType base_size = axis_size / NumberOfPartitions;
        const IndexType remainder = axis_size % NumberOfPartitions;
        IndexType begin = 0;
        for (IndexType partition_index = 0; partition_index < NumberOfPartitions; ++partition_index) {
            const IndexType size = base_size + (partition_index < remainder ? 1 : 0);
            Vector3i lower{ 0, 0, 0 };
            Vector3i upper{ rNumberOfElements[0] - 1, rNumberOfElements[1] - 1, rNumberOfElements[2] - 1 };
            lower[PartitionAxisIndex] = begin;
            upper[PartitionAxisIndex] = begin + size - 1;
            partitions.emplace_back(lower, upper);
            begin += size;
        }
        return partitions;
    }

    /// @brief Returns whether one cell is excluded from fill as trimmed.
    [[nodiscard]] bool IsTrimmed(const detail::DomainTopology& rTopology, IndexType CellIndex) noexcept
    { return rTopology.cell_topology[CellIndex].is_trimmed; }

    /// @brief Validates topology invariants required by serial and partitioned flood fill.
    /// @details Validation is performed serially before OpenMP traversal. The topology must reference the supplied
    ///          indexer, contain one cell record per grid cell, block every exterior face, and block every face
    ///          adjacent to a trimmed cell.
    void ValidateTopology(const GridIndexer& rGridIndexer, const detail::DomainTopology& rTopology)
    {
        QuESo_ERROR_IF(rTopology.cell_topology.size() != rGridIndexer.NumberOfElements())
            << "Domain topology and GridIndexer must describe the same number of cells.\n";

        for (IndexType cell_index = 0; cell_index < rGridIndexer.NumberOfElements(); ++cell_index) {
            for (const Direction direction : EnumRange<Direction>()) {
                const auto face_state = rTopology.face_topology.Get(cell_index, direction);
                QuESo_ERROR_IF(
                    rGridIndexer.IsEnd(cell_index, direction) && face_state != detail::GridFaceState::blocked
                ) << "Exterior grid faces must be blocked before flood fill.\n";
                QuESo_ERROR_IF(IsTrimmed(rTopology, cell_index) && face_state != detail::GridFaceState::blocked)
                    << "Every face adjacent to a trimmed cell must be blocked before flood fill.\n";
            }
        }
    }

    /// @brief Processes one directed cell face and optionally discovers a fillable neighbor.
    /// @details Physical exterior faces record conclusive exterior contact. Internal partition boundaries are deferred
    ///          to reconciliation. Blocked faces contribute their directed vote, while traversable faces may add one
    ///          previously unvisited neighbor to the component.
    [[nodiscard]] std::optional<IndexType> Move(
        const GridIndexer& rGridIndexer,
        const detail::DomainTopology& rTopology,
        IndexType CellIndex,
        Direction FaceDirection,
        const PartitionBox& rPartition,
        Component& rComponent,
        std::vector<VisitState>& rVisitStates
    )
    {
        // A physical grid boundary is also a slab boundary, but it must first mark the component as exterior-connected.
        if (rGridIndexer.IsEnd(CellIndex, FaceDirection)) {
            QuESo_ASSERT(
                rTopology.face_topology.Get(CellIndex, FaceDirection) == detail::GridFaceState::blocked,
                "Exterior grid faces must be blocked."
            );
            rComponent.touches_grid_exterior = true;
            return std::nullopt;
        }
        if (rGridIndexer.IsEnd(CellIndex, FaceDirection, rPartition)) { return std::nullopt; }

        const auto [next_index, index_info] = rGridIndexer.GetNextIndex(CellIndex, FaceDirection);
        QuESo_ASSERT(index_info == GridIndexer::IndexInfo::middle, "Expected a physical neighboring grid cell.");
        if (rTopology.face_topology.Get(CellIndex, FaceDirection) == detail::GridFaceState::blocked) {
            rComponent.classification_vote += rTopology.face_votes.Get(CellIndex, FaceDirection);
            return std::nullopt;
        }

        QuESo_ASSERT(!IsTrimmed(rTopology, next_index), "A traversable face cannot lead into a trimmed cell.");
        if (rVisitStates[next_index] == VisitState::visited) { return std::nullopt; }
        rVisitStates[next_index] = VisitState::visited;
        rComponent.cell_indices.push_back(next_index);
        return next_index;
    }

    /// @brief Discovers one connected component within a single partition.
    /// @details Uses iterative depth-first traversal and never crosses the supplied slab boundary.
    void Fill(
        const GridIndexer& rGridIndexer,
        const detail::DomainTopology& rTopology,
        IndexType CellIndex,
        const PartitionBox& rPartition,
        Component& rComponent,
        std::vector<VisitState>& rVisitStates
    )
    {
        rVisitStates[CellIndex] = VisitState::visited;
        rComponent.cell_indices.push_back(CellIndex);
        std::stack<IndexType> pending_indices;
        pending_indices.push(CellIndex);
        while (!pending_indices.empty()) {
            const IndexType current_index = pending_indices.top();
            pending_indices.pop();
            for (const Direction direction : EnumRange<Direction>()) {
                if (const auto next_index =
                        Move(rGridIndexer, rTopology, current_index, direction, rPartition, rComponent, rVisitStates)) {
                    pending_indices.push(*next_index);
                }
            }
        }
    }

    /// @brief Discovers independent component fragments inside all flood-fill partitions.
    /// @details Partitions are processed independently and may run concurrently because traversal cannot cross slab
    ///          interfaces. Fragment vectors are concatenated deterministically in partition order.
    [[nodiscard]] Components PartitionedFill(
        const GridIndexer& rGridIndexer,
        const detail::DomainTopology& rTopology,
        const PartitionBoxes& rPartitions,
        std::vector<VisitState>& rVisitStates
    )
    {
        std::vector<Components> components_by_partition(rPartitions.size());
        // Slabs are disjoint and Move never crosses a partition boundary, so workers mutate disjoint visit-state cells.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1)
#endif
        for (std::int64_t raw_partition_index = 0; raw_partition_index < static_cast<std::int64_t>(rPartitions.size());
             ++raw_partition_index) {
            const IndexType partition_index = static_cast<IndexType>(raw_partition_index);
            const auto& r_partition = rPartitions[partition_index];
            auto& r_components = components_by_partition[partition_index];
            for (IndexType z = r_partition.first[2]; z <= r_partition.second[2]; ++z) {
                for (IndexType y = r_partition.first[1]; y <= r_partition.second[1]; ++y) {
                    for (IndexType x = r_partition.first[0]; x <= r_partition.second[0]; ++x) {
                        const IndexType cell_index = rGridIndexer.GetVectorIndexFromMatrixIndices(x, y, z);
                        if (rVisitStates[cell_index] == VisitState::visited) { continue; }
                        Component component;
                        Fill(rGridIndexer, rTopology, cell_index, r_partition, component, rVisitStates);
                        r_components.push_back(std::move(component));
                    }
                }
            }
        }

        Components components;
        for (auto& r_partition_components : components_by_partition) {
            components.insert(
                components.end(),
                std::make_move_iterator(r_partition_components.begin()),
                std::make_move_iterator(r_partition_components.end())
            );
        }
        return components;
    }

    /// @brief Reconciles partition-local fragments into global connected components.
    /// @details Traversable interfaces union fragments. Blocked interfaces keep fragments separate and add each
    ///          directed side's vote to its adjacent component before merged vote and exterior state are accumulated.
    [[nodiscard]] Components MergeComponents(
        const GridIndexer& rGridIndexer,
        const detail::DomainTopology& rTopology,
        Components& rComponents,
        IndexType PartitionAxisIndex,
        const PartitionBoxes& rPartitions
    )
    {
        const IndexType invalid_component = std::numeric_limits<IndexType>::max();
        std::vector<IndexType> component_indices(rGridIndexer.NumberOfElements(), invalid_component);
        for (IndexType component_index = 0; component_index < rComponents.size(); ++component_index) {
            for (const IndexType cell_index : rComponents[component_index].cell_indices) {
                component_indices[cell_index] = component_index;
            }
        }

        DisjointSet disjoint_set(rComponents.size());
        const Direction forward_direction = static_cast<Direction>(2 * PartitionAxisIndex);
        const Direction backward_direction = GridIndexer::ReverseDirection(forward_direction);
        // Traversable interfaces join local fragments. Blocked interfaces retain both directed votes separately.
        for (IndexType partition_index = 0; partition_index + 1 < rPartitions.size(); ++partition_index) {
            const auto& r_partition = rPartitions[partition_index];
            for (IndexType z = r_partition.first[2]; z <= r_partition.second[2]; ++z) {
                for (IndexType y = r_partition.first[1]; y <= r_partition.second[1]; ++y) {
                    for (IndexType x = r_partition.first[0]; x <= r_partition.second[0]; ++x) {
                        Vector3i indices{ x, y, z };
                        if (indices[PartitionAxisIndex] != r_partition.second[PartitionAxisIndex]) { continue; }

                        const IndexType left_index = rGridIndexer.GetVectorIndexFromMatrixIndices(indices);
                        const auto [right_index, index_info] = rGridIndexer.GetNextIndex(left_index, forward_direction);
                        QuESo_ASSERT(
                            index_info == GridIndexer::IndexInfo::middle,
                            "Adjacent flood-fill partitions must share a physical grid face."
                        );
                        const IndexType left_component = component_indices[left_index];
                        const IndexType right_component = component_indices[right_index];
                        if (rTopology.face_topology.Get(left_index, forward_direction)
                            == detail::GridFaceState::traversable) {
                            QuESo_ASSERT(
                                left_component != invalid_component && right_component != invalid_component,
                                "A traversable partition face must connect two fillable cells."
                            );
                            disjoint_set.Union(left_component, right_component);
                            continue;
                        }

                        if (left_component != invalid_component) {
                            rComponents[left_component].classification_vote +=
                                rTopology.face_votes.Get(left_index, forward_direction);
                        }
                        if (right_component != invalid_component) {
                            rComponents[right_component].classification_vote +=
                                rTopology.face_votes.Get(right_index, backward_direction);
                        }
                    }
                }
            }
        }

        Components merged_components;
        std::vector<IndexType> merged_indices(rComponents.size(), invalid_component);
        for (IndexType component_index = 0; component_index < rComponents.size(); ++component_index) {
            const IndexType root = disjoint_set.Find(component_index);
            if (merged_indices[root] == invalid_component) {
                merged_indices[root] = merged_components.size();
                merged_components.emplace_back();
            }
            auto& r_merged = merged_components[merged_indices[root]];
            auto& r_component = rComponents[component_index];
            r_merged.cell_indices.insert(
                r_merged.cell_indices.end(), r_component.cell_indices.begin(), r_component.cell_indices.end()
            );
            r_merged.classification_vote += r_component.classification_vote;
            r_merged.touches_grid_exterior = r_merged.touches_grid_exterior || r_component.touches_grid_exterior;
        }
        return merged_components;
    }

}  // namespace

namespace detail {

    ElementStates FloodFill(const DomainTopology& rTopology, FloodFillOptions Options)
    {
        const GridIndexer& rGridIndexer = rTopology.GetGridIndexer();
        const IndexType number_of_cells = rGridIndexer.NumberOfElements();
        QuESo_ERROR_IF(number_of_cells == 0) << "Flood fill requires a non-empty grid.\n";
        ValidateTopology(rGridIndexer, rTopology);

        ElementStates states(number_of_cells, IntersectionState::outside);
        std::vector<VisitState> visit_states(number_of_cells, VisitState::unvisited);
        // Trimmed cells are final states, not fillable component members. Topology validation guarantees blocked
        // neighbors.
        for (IndexType cell_index = 0; cell_index < number_of_cells; ++cell_index) {
            if (IsTrimmed(rTopology, cell_index)) {
                states[cell_index] = IntersectionState::trimmed;
                visit_states[cell_index] = VisitState::visited;
            }
        }

        const Vector3i number_of_elements = rGridIndexer.ElementCounts();
        const IndexType partition_axis = PartitionAxis(number_of_elements);
        const auto partitions = CreatePartitions(
            number_of_elements, partition_axis, PartitionCount(number_of_elements[partition_axis], Options)
        );
        auto components = PartitionedFill(rGridIndexer, rTopology, partitions, visit_states);
        const auto merged_components = MergeComponents(rGridIndexer, rTopology, components, partition_axis, partitions);

        // Physical exterior contact is conclusive outside evidence. Only enclosed components use their signed vote sum.
        for (const auto& r_component : merged_components) {
            const IntersectionState state = !r_component.touches_grid_exterior && r_component.classification_vote > 0
                                                ? IntersectionState::inside
                                                : IntersectionState::outside;
            for (const IndexType cell_index : r_component.cell_indices) { states[cell_index] = state; }
        }
        return states;
    }

}  // namespace detail

}  // namespace queso::embedding
