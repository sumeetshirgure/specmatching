# Edits to vendored upstream files

`src/pyrematching/` below `two_phase/` is vendored PyMatching, renamed (M0). Everything under
`src/pyrematching/two_phase/`, `tests/two_phase/` and `benchmarks/two_phase/` is ours and is not
listed here.

This file exists so that "which lines are ours" stays answerable by path alone, for the licence,
for the paper, and for anyone who ever wants to diff this fork against upstream again. **Every edit
to a file outside `two_phase/` gets a row.** The design budgets three such files; anything beyond
that needs a justification in the PR as well as a row here.

| File | Milestone | Edit | Why |
|---|---|---|---|
| `sparse_blossom/tracker/queued_event_tracker.h` | M1.1 | `set_desired_event` gains `cumulative_event_time` and `horizon` arguments, both defaulted so existing call sites are unchanged; adds the `pm::NO_HORIZON` sentinel and the debug-only `pm::horizon_gate_stats` counters | This is the funnel every scheduled event passes through, so it is the one place a horizon gate can guarantee that nothing past `T` is *inserted* rather than inserted and later skipped. The cumulative time has to come from the caller: `ev.time` is a `cyclic_time_int` and cannot be ordered against a cumulative bound |
| `sparse_blossom/flooder/graph_flooder.h` | M1.1 | Adds the `cumulative_time_int horizon` field, initialised to `pm::NO_HORIZON` | Planned by the design. The flooder is what the scheduling sites have in hand; the sentinel means the stock path never observes a finite horizon |
| `sparse_blossom/flooder/graph_flooder.cc` | M1.1 | Passes the horizon at the three node-event scheduling sites (`reschedule_events_at_detector_node`, both sites in `do_look_at_node_event`); leaves `schedule_tentative_shrink_event` un-gated; copies `horizon` in the move constructor | Planned by the design. Each site carries a one-line comment saying gate-or-exempt and why. Shrink events are exempt so that a shrinking region's shell-area bookkeeping can always complete |
| `sparse_blossom/flooder/graph_flooder.cc` | M1.1 | **Beyond the design's three rows.** `run_until_next_mwpm_notification` stops and reports `NO_EVENT` once `queue.cur_time > horizon` | Consequence of the shrink exemption. Exempt shrink events *can* be scheduled past `T`, and processing one would carry the timeline — and therefore every surviving dual — past the horizon, breaking `Y(u) == T`. With the sentinel horizon the comparison is never true, and the stock-path perf gate is green |
| `sparse_blossom/driver/mwpm_decoding.cc` | M1.2 | The detection-event/negative-weight preamble of `process_timeline_until_completion` is factored out into `pm::begin_timeline`, which it now calls | Required by M1.2: the truncated timeline reuses the preamble verbatim rather than copying it. `process_timeline_until_completion`'s behaviour is unchanged |
| `sparse_blossom/driver/mwpm_decoding.h` | M1.2 | Declares `pm::begin_timeline` | Same edit, other half |
| `CMakeLists.txt` | M1 | Adds the `two_phase/` sources to `SOURCE_FILES_NO_MAIN`, the `tests/two_phase/` suites to `TEST_FILES`, the repo root to the test target's include path, `-UNDEBUG` on the test target, and the `two_phase_m1_artifact` target | The file enumerates sources explicitly, so new files have to be listed. `-UNDEBUG` keeps the debug-build invariants (asserts, horizon gate counters) live in the test binary even when the tree is configured `Release` |

## Not edited, deliberately

- The `pm` / `pm_pybind` namespaces are untouched, as the design requires. New code lives in
  `pm::two_phase`.
- `sparse_blossom/matcher/mwpm.{h,cc}` is untouched. `shatter_exposed_blossom_and_extract_matches`
  is a *sibling* of the upstream shatter routines, in `two_phase/truncation/exposed_blossom.cc`,
  rather than a modification of them.
- `sparse_blossom/ints.h` is untouched. `pm::NO_HORIZON` went into `queued_event_tracker.h`, which
  was already being edited, rather than opening a fourth vendored file.
