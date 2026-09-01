// Copyright 2022 PyMatching Contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef SPECMATCHING_FILL_MATCH_ARENA_H
#define SPECMATCHING_FILL_MATCH_ARENA_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

namespace pm {

/// Slot value meaning "this object is not currently in its arena's live vector".
constexpr uint32_t ARENA_NOT_LIVE = UINT32_MAX;

/// Detects types that opt into arena liveness tracking (specmatching §M2.9.1) by declaring
///
///     uint32_t arena_live_index = pm::ARENA_NOT_LIVE;
///
/// Arenas of types that do not declare it are byte-for-byte and instruction-for-instruction what
/// they were before, because every use below is behind `if constexpr`.
template <typename T, typename = void>
struct arena_tracks_live : std::false_type {};
template <typename T>
struct arena_tracks_live<T, std::void_t<decltype(std::declval<T &>().arena_live_index)>> : std::true_type {};

/// World's simplest bulk memory owner.
///
/// Memory allocated by the arena is free'd when the arena is destructed.
///
/// For types that opt in (see `arena_tracks_live`), the arena additionally keeps `live`: a compact
/// vector of the objects currently checked out, so that "everything alive right now" is enumerable
/// in one linear scan, with no traversal of whatever structure those objects happen to form. Every
/// allocation and every deallocation of such an object funnels through this type, which is what
/// makes the vector impossible to leave stale: there is no mutation site to miss.
///
/// Membership is maintained by swap-remove, which is O(1) and needs the object to remember its own
/// slot — one `uint32_t`, which on both tracked types fits in padding the object already had. An
/// intrusive doubly-linked list would have been the more obvious shape and is what §M2.9.1
/// suggests, but it costs a pointer pair per object, and on this machine those 16 extra bytes cost
/// a measurable ~3% of the stock decode path purely in cache footprint. A vector is also the better
/// answer for what the list is *for*: enumeration is a linear scan of contiguous pointers rather
/// than a chase.
///
/// Note that `alloc_unconstructed` deliberately does *not* enrol the object, because the caller has
/// not constructed it yet and the placement-new would overwrite the slot; callers that want a
/// tracked object use `alloc_constructed` or `alloc_default_constructed`. `del` handles both cases.
template <typename T>
struct Arena {
    std::vector<T *> allocated;
    std::vector<T *> available;
    /// The objects currently checked out, in no particular order. Always empty for types that do
    /// not opt into liveness tracking. Keeps its capacity across shots, so a steady-state shot does
    /// not allocate here.
    std::vector<T *> live;

    Arena() : allocated(), available() {
    }
    Arena(const Arena &) = delete;
    Arena(Arena &&other)
        : allocated(std::move(other.allocated)),
          available(std::move(other.available)),
          live(std::move(other.live)) {
    }

    T *alloc_unconstructed() {
        if (available.empty()) {
            T *p = (T *)malloc(sizeof(T));
            allocated.push_back(p);
            available.push_back(p);
        }
        T *result = available.back();
        available.pop_back();
        return result;
    }

    T *alloc_default_constructed() {
        T *result = alloc_unconstructed();
        new (result) T();
        link_live(result);
        return result;
    }

    /// Allocates and constructs in one step, so that the liveness link can be installed after the
    /// object's lifetime has begun. This is the tracked equivalent of `alloc_unconstructed` followed
    /// by a placement-new.
    template <typename... Args>
    T *alloc_constructed(Args &&...args) {
        T *result = alloc_unconstructed();
        new (result) T(std::forward<Args>(args)...);
        link_live(result);
        return result;
    }

    void del(T *p) {
        unlink_live(p);
        available.push_back(p);
        p->~T();
    }

    inline void link_live(T *p) {
        if constexpr (arena_tracks_live<T>::value) {
            p->arena_live_index = (uint32_t)live.size();
            live.push_back(p);
        } else {
            (void)p;
        }
    }

    inline void unlink_live(T *p) {
        if constexpr (arena_tracks_live<T>::value) {
            // Objects handed out by `alloc_unconstructed` were never enrolled, so deleting one is a
            // no-op here rather than a corruption of somebody else's slot.
            if (p->arena_live_index == ARENA_NOT_LIVE)
                return;
            uint32_t slot = p->arena_live_index;
            T *moved = live.back();
            live[slot] = moved;
            moved->arena_live_index = slot;
            live.pop_back();
            // Last, so that the `moved == p` case ends up marked not-live rather than pointing at
            // the slot it was just removed from.
            p->arena_live_index = ARENA_NOT_LIVE;
        } else {
            (void)p;
        }
    }

    ~Arena() {
        // `Mwpm::reset()` invokes this destructor explicitly and then lets the object be destroyed
        // again on scope exit, so every member has to survive being destructed twice. `allocated`
        // and `available` do, because they are moved out below and a moved-from vector owns no
        // buffer. `live` is emptied the same way, for the same reason — clearing it would keep the
        // capacity and the second destruction would then double-free it.
        live = std::vector<T *>();
        std::vector<T *> to_free = std::move(allocated);
        std::vector<T *> not_in_use = std::move(available);

        // Destruct the objects that were still in use.
        std::sort(to_free.begin(), to_free.end());
        std::sort(not_in_use.begin(), not_in_use.end());
        size_t kf = 0;
        size_t kn = 0;
        while (kf < to_free.size()) {
            if (kn < not_in_use.size() && not_in_use[kn] == to_free[kf]) {
                kf++;
                kn++;
            } else {
                to_free[kf]->~T();
                kf++;
            }
        }

        // Free all allocated memory.
        for (T *v : to_free) {
            free(v);
        }
    }
};
}  // namespace pm

#endif