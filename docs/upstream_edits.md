# Edits to vendored upstream files

`src/pyrematching/` below `two_phase/` is vendored PyMatching, renamed (M0). Everything under
`src/pyrematching/two_phase/`, `tests/two_phase/` and `benchmarks/two_phase/` is ours and is not
listed here.

This file exists so that "which lines are ours" stays answerable by path alone, for the licence,
for the paper, and for anyone who ever wants to diff this fork against upstream again. **Every edit
to a file outside `two_phase/` gets a row.** The design budgets three such files; anything beyond
that needs a justification in the PR as well as a row here.

The budget is exceeded, deliberately and twice. M1 needed four files and said so; M2.9 needs five
more (`arena.h`, `alternating_tree.h`, `graph_fill_region.h`, `mwpm.cc`, `graph_flooder.{h,cc}`),
against a design that anticipated "a fourth vendored-file edit (`matcher/alternating_tree.h`,
possibly `matcher/mwpm.cc`)". The trade taken in each case is *fewer mutation sites*, not fewer
files: liveness is maintained inside `Arena<T>`, which is the single funnel both allocations and
deallocations already pass through, rather than at the seven `add_match` / `Match::clear` sites in
`mwpm.cc` that §M2.9.2's literal reading would need. The design's own guard rail is the call-site
count — "if the intrusive list turns out to need edits at more than two call sites, prefer an arena
sweep with a liveness bit" — and it is the count, not the file count, that the silent-failure risk
scales with. See `docs/two_phase_m29_exit.md`.

| File | Milestone | Edit | Why |
|---|---|---|---|
| `sparse_blossom/tracker/queued_event_tracker.h` | M1.1 | `set_desired_event` gains `cumulative_event_time` and `horizon` arguments, both defaulted so existing call sites are unchanged; adds the `pm::NO_HORIZON` sentinel and the debug-only `pm::horizon_gate_stats` counters | This is the funnel every scheduled event passes through, so it is the one place a horizon gate can guarantee that nothing past `T` is *inserted* rather than inserted and later skipped. The cumulative time has to come from the caller: `ev.time` is a `cyclic_time_int` and cannot be ordered against a cumulative bound |
| `sparse_blossom/flooder/graph_flooder.h` | M1.1 | Adds the `cumulative_time_int horizon` field, initialised to `pm::NO_HORIZON` | Planned by the design. The flooder is what the scheduling sites have in hand; the sentinel means the stock path never observes a finite horizon |
| `sparse_blossom/flooder/graph_flooder.cc` | M1.1 | Passes the horizon at the three node-event scheduling sites (`reschedule_events_at_detector_node`, both sites in `do_look_at_node_event`); leaves `schedule_tentative_shrink_event` un-gated; copies `horizon` in the move constructor | Planned by the design. Each site carries a one-line comment saying gate-or-exempt and why. Shrink events are exempt so that a shrinking region's shell-area bookkeeping can always complete |
| `sparse_blossom/flooder/graph_flooder.cc` | M1.1 | **Beyond the design's three rows.** `run_until_next_mwpm_notification` stops and reports `NO_EVENT` once `queue.cur_time > horizon` | Consequence of the shrink exemption. Exempt shrink events *can* be scheduled past `T`, and processing one would carry the timeline — and therefore every surviving dual — past the horizon, breaking `Y(u) == T`. With the sentinel horizon the comparison is never true, and the stock-path perf gate is green |
| `sparse_blossom/driver/mwpm_decoding.cc` | M1.2 | The detection-event/negative-weight preamble of `process_timeline_until_completion` is factored out into `pm::begin_timeline`, which it now calls | Required by M1.2: the truncated timeline reuses the preamble verbatim rather than copying it. `process_timeline_until_completion`'s behaviour is unchanged |
| `sparse_blossom/driver/mwpm_decoding.h` | M1.2 | Declares `pm::begin_timeline` | Same edit, other half |
| `CMakeLists.txt` | M1 | Adds the `two_phase/` sources to `SOURCE_FILES_NO_MAIN`, the `tests/two_phase/` suites to `TEST_FILES`, the repo root to the test target's include path, `-UNDEBUG` on the test target, and the `two_phase_m1_artifact` target | The file enumerates sources explicitly, so new files have to be listed. `-UNDEBUG` keeps the debug-build invariants (asserts, horizon gate counters) live in the test binary even when the tree is configured `Release` |
| `sparse_blossom/arena.h` | M2.9.1 | Adds the `pm::arena_tracks_live<T>` trait, a `std::vector<T *> live` of the currently checked-out objects maintained by swap-remove, `alloc_constructed(Args&&...)`, and the link/unlink calls in `alloc_default_constructed` and `del`. Also empties `live` by move in `~Arena` rather than by `clear()`. Every use is behind `if constexpr`, so arenas of types that do not opt in are unchanged | §M2.9.1 needs to enumerate the live `AltTreeNode`s, and §M2.9.2 the live `GraphFillRegion`s, without traversing the structures they form. The design's fallback for that is "an arena sweep with a liveness bit", and this is that: the arena is the **single funnel** every allocation and deallocation of both types already passes through, so the list cannot be left stale by a missed mutation site — which is the one failure mode the design flags as silent. Linking has to happen after construction (a placement-new would clobber the slot), which is why `alloc_constructed` exists rather than a hook inside `alloc_unconstructed`. The destructor detail is not cosmetic: `Mwpm::reset()` destructs the arena explicitly and then lets it be destructed again, which the pre-existing members survive only because they are moved-from — a `clear()` keeps the capacity and double-frees on the second pass (caught by ASan in `MwpmDecoding.NoValidSolutionForLineGraph`) |
| `sparse_blossom/matcher/alternating_tree.h` | M2.9.1 | Adds one `uint32_t arena_live_index = pm::ARENA_NOT_LIVE` to `AltTreeNode`, with a default member initializer | Declaring it is what opts the type into the arena's liveness tracking, and it is the slot swap-remove needs to make removal O(1). **Not** the pointer pair §M2.9.1 budgets: `prev`/`next` grew `AltTreeNode` from 104 bytes to 128 and `GraphFillRegion` from 128 to 152, and that cost a consistent 2–4% of the stock decode path in cache footprint alone, failing the perf gate. The single slot lands in padding both structs already had, so both are byte-for-byte their previous size. A vector is also the better shape for what the list is *for* — enumeration is a linear scan of contiguous pointers, not a chase. Nothing in the matcher reads or writes the field, and no constructor mentions it |
| `sparse_blossom/flooder/graph_fill_region.h` | M2.9.2 | The same field on `GraphFillRegion`, declared at offset 44 between `shrink_event_tracker` and `match` | §M2.9.2 asks for a maintained list of top-level matched regions, appended at the match sites. Doing it that way needs edits at seven sites in `mwpm.cc` (four `add_match` calls, the direct boundary match, and the two `Match::clear`s in `handle_tree_hitting_match`), which is three and a half times the design's own "more than two call sites" threshold, and a missed one is silent. Tracking *every* region instead reduces "top-level" to one local `blossom_parent == nullptr` test, costs no mutation sites at all, and additionally turns `dual_sum_at_truncation` from a descent through blossom nesting into a flat reduction. The declaration position is load-bearing — that hole is the existing padding the field has to land in, so do not tidy it to the bottom of the struct |
| `sparse_blossom/matcher/mwpm.cc` | M2.9.1 | `make_child` and `create_detection_event` allocate through `node_arena.alloc_constructed(...)` instead of `alloc_unconstructed()` followed by a placement-new | The two sites where an `AltTreeNode` comes into existence. Behaviourally identical; it exists so the arena can install the liveness link once the object's lifetime has begun |
| `sparse_blossom/flooder/graph_flooder.h` | M2.9.6 | Adds `pm::BlossomFormationStats` and the process-wide `pm::blossom_formation_stats` counter | §M2.9.6 measurement 3. §M2.9.4's eager-cached-base option is worth its field only if blossom formations are much rarer than base readouts, and that ratio is not recoverable after the fact — a blossom that forms and then shatters leaves nothing to count. Follows the `pm::horizon_gate_stats` precedent, except that this one is live in release builds, because the measurement it feeds is a release measurement |
| `sparse_blossom/flooder/graph_flooder.cc` | M2.9.6 | `create_blossom` increments `blossom_formation_stats.formations` | The one place a blossom is created. One increment, against a formation that already sweeps the blossom's whole area to reschedule events |
| `CMakeLists.txt` | M2.9 | Adds `two_phase/truncation/split_experiment.cc`, the two new `tests/two_phase/` suites, and the `two_phase_m29_artifact` target | Same reason as the M1 row: the file enumerates sources explicitly |

## Not edited, deliberately

- The `pm` / `pm_pybind` namespaces are untouched, as the design requires. New code lives in
  `pm::two_phase`.
- `sparse_blossom/matcher/mwpm.{h,cc}` is untouched. `shatter_exposed_blossom_and_extract_matches`
  is a *sibling* of the upstream shatter routines, in `two_phase/truncation/exposed_blossom.cc`,
  rather than a modification of them.
- `sparse_blossom/ints.h` is untouched. `pm::NO_HORIZON` went into `queued_event_tracker.h`, which
  was already being edited, rather than opening a fourth vendored file.
