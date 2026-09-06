# Shot-level edge list construction in hardware

## Who this is for

You are going to build a hardware block (Verilog, VHDL, Chisel, or a high-level
synthesis description — the choice is yours) that does one job that is currently
done in software. This document tells you what that job is, what data goes in,
and what data must come out. It does not tell you how to build it.

Everything described here already exists as C++ code. The single function that
this hardware block must replace is:

    build_ball_graph(...)

- Declared in `src/specmatching/spec_matching/manifold/ball_graph.h:332`
- Implemented in `src/specmatching/spec_matching/manifold/ball_graph.cc:149`

The data types it reads and writes are in the same two headers:

- `src/specmatching/spec_matching/manifold/ball_tables.h:71` — the big
  precomputed input tables (`struct BallTables`).
- `src/specmatching/spec_matching/manifold/ball_graph.h:36` — the output edge
  record (`struct BallGraphEdge`).
- `src/specmatching/spec_matching/manifold/ball_graph.h:46` — the output
  boundary-edge record (`struct BallBoundaryEdge`).
- `src/specmatching/spec_matching/manifold/ball_graph.h:57` — the container
  holding the whole output (`struct BallGraph`).

Read those four structs first. They are short, and they are the contract.

## The job in plain words

There is a large fixed graph called `G`. Its vertices are called **detectors**
and they are numbered `0, 1, 2, ... num_nodes - 1`. This graph never changes
while the machine runs. It is known before the machine starts.

Every time the quantum computer is measured, we get one **shot**. A shot tells
us that a small subset of the detectors "fired". A shot typically has $d^3$ detectors.
A detector that fired is called
a **detection event**, or a **defect**. Typically only a small fraction of all
detectors fire in any one shot.

For each shot we must build a small temporary graph called `H`:

- The vertices of `H` are exactly the detectors that fired in this shot.
- There is an edge in `H` between two fired detectors if and only if the
  distance between them in the big graph `G` is at most a threshold `2T`.
- Each such edge carries the exact distance between the two detectors as its
  weight.
- Separately, a fired detector also gets a **boundary edge** if its own
  precomputed distance to the boundary is at most `T`.

That is the whole job. Nothing else in the decoder is in scope.

The reason this can be done fast in hardware is that all the distances are
precomputed. For each detector `v`, we already have a stored list of every
detector within a radius `R` of `v`, together with the exact distance. So the
per-shot work is not a distance computation. It is an **intersection**: take
detector `v`'s precomputed neighbour list, intersect it with the set of
detectors that fired in this shot, and keep the ones that survive and are close
enough. Then repeat for every detector that fired.

## Inputs

There are two kinds of input. The first kind is loaded once and then never
changes. The second kind arrives fresh with every shot.

### Input group 1 — the precomputed tables (loaded once, read-only)

These come from `struct BallTables` in
`src/specmatching/spec_matching/manifold/ball_tables.h:71`. They are built
offline by software, before the hardware starts running, and they are the same
for every shot.

They are stored in "CSR" form. CSR just means: one long flat array holding all
the data back to back, plus an offsets array that says where each detector's
slice starts and ends. Detector `v`'s slice of the flat array `X` is
`X[offsets[v]]` up to but not including `X[offsets[v + 1]]`.

Scalar values:

| Name | C++ type | Meaning |
| --- | --- | --- |
| `num_nodes` | `size_t` | How many detectors exist in the big graph `G`. |
| `num_observables` | `size_t` | How many logical observables exist. Used only to size the mask ids below. |
| `r_int` | `int64_t` | The radius the tables were built with, in integer distance units. |
| `t_max_int` | `int64_t` | The largest horizon `T` the tables support. A shot asking for more than this is an error. |

The neighbour lists ("balls"). For every detector `v`, its ball is the list of
all other detectors within distance `R` of `v`. Detector `v` itself is not in
its own ball.

| Name | C++ type | Meaning |
| --- | --- | --- |
| `ball_offsets` | `vector<uint64_t>` | Length `num_nodes + 1`. Where detector `v`'s ball starts and ends in the two arrays below. |
| `ball_target` | `vector<uint32_t>` | The detector id of each neighbour. |
| `ball_w_int` | `vector<weight_int>` (`uint32_t`) | The exact integer distance from `v` to that neighbour. |

Important ordering fact: within one detector's ball, the entries are sorted by
`(distance, target id)`. Distance ascending first. This matters, because it
means a walk down the list can stop at the first entry whose distance exceeds
the threshold — everything after it is also too far.

The same balls, stored a second time as bitmaps. This second copy exists because
it makes the intersection a plain bitwise AND.

| Name | C++ type | Meaning |
| --- | --- | --- |
| `ball_word_offsets` | `vector<uint64_t>` | Length `num_nodes + 1`. Where detector `v`'s bitmap starts and ends in `ball_words`. |
| `ball_word_base` | `vector<uint32_t>` | Which 64-bit word of the global detector space detector `v`'s bitmap starts at. So `v`'s bitmap covers detector ids `64 * ball_word_base[v]` upward. |
| `ball_words` | `vector<uint64_t>` | The bitmap itself. Bit `k` of a word is set if that detector is in `v`'s ball. |
| `ball_word_rank` | `vector<uint32_t>` | For each word of `v`'s bitmap, how many bits are set in all the *earlier* words of `v`'s bitmap. This turns a set bit into a position number. |
| `ball_entry_by_rank` | `vector<uint64_t>` | Maps that position number back to the index into `ball_target` / `ball_w_int`. Needed because the bitmap yields hits in detector-id order, while the arrays above are in distance order. |

The boundary tables. These say how far each detector is from the boundary of the
code.

| Name | C++ type | Meaning |
| --- | --- | --- |
| `has_bcost` | `vector<uint8_t>` | 1 if detector `v` has a boundary path within radius `R`, 0 if not. |
| `bcost_w_int` | `vector<weight_int>` (`uint32_t`) | The exact integer distance from `v` to the boundary. Only meaningful when `has_bcost[v]` is 1. |

The observable masks. Every ball entry also records which logical observables
the shortest path between that pair crosses. **The hardware block described here
does not need to read these.** It only needs to output the index that points at
them, so that a later stage can fetch them. They are listed here so you know
what the `entry` output field below is for.

| Name | C++ type | Meaning |
| --- | --- | --- |
| `ball_mask_offsets` | `vector<uint64_t>` | Where entry `e`'s observable id list starts and ends. |
| `ball_mask_ids` | `vector<uint32_t>` | The observable ids themselves. |
| `bcost_mask_offsets`, `bcost_mask_ids` | same | The same thing for boundary paths, indexed by detector id. |

### Input group 2 — the per-shot inputs (new every shot)

| Name | C++ type | Meaning |
| --- | --- | --- |
| `seeded_dets` | `vector<uint64_t>` | The detectors that fired in this shot. **Sorted ascending, no duplicates.** Every value is less than `num_nodes`. |
| `horizon` | `int64_t` (`horizon_int`) | The threshold `T`, in the same integer distance units as `ball_w_int`. Must be greater than zero and not greater than `t_max_int`. |

The software derives `two_t = 2 * horizon` once at
`src/specmatching/spec_matching/manifold/ball_graph.cc:174` and compares
distances against that.

A note on where `seeded_dets` comes from: it is not the raw measurement record.
Software preprocesses the raw shot before this point (see
`src/specmatching/spec_matching/manifold/ball_decoding.cc:349`). Whether that
preprocessing also moves into hardware is a separate question. For this block,
assume you are handed the already-cleaned, already-sorted list.

## Outputs

The output is the small per-shot graph `H`. It has three parts.

### Output part 1 — the node list

| Name | C++ type | Meaning |
| --- | --- | --- |
| `h_to_det` | `vector<uint64_t>` | The detector id of each node of `H`, ascending. |

This is just a copy of `seeded_dets`. Its real purpose is to define the
**local numbering**: node `0` of `H` is the first fired detector, node `1` is the
second, and so on. **_All the edges below are expressed in this local numbering,
not in detector ids._**

### Output part 2 — the edge list

A list of records, each one being `struct BallGraphEdge`
(`src/specmatching/spec_matching/manifold/ball_graph.h:36`):

| Field | C++ type | Meaning |
| --- | --- | --- |
| `i` | `uint32_t` | Local index of the first endpoint (a position in `h_to_det`). |
| `j` | `uint32_t` | Local index of the second endpoint. |
| `w_int` | `weight_int` (`uint32_t`) | The exact integer distance between them. Copied straight out of `ball_w_int`. |
| `entry` | `uint64_t` | The index into the ball arrays that this edge came from. A later stage uses it to look up the observable mask. |

### Output part 3 — the boundary edge list

A list of records, each one being `struct BallBoundaryEdge`
(`src/specmatching/spec_matching/manifold/ball_graph.h:46`):

| Field | C++ type | Meaning |
| --- | --- | --- |
| `i` | `uint32_t` | Local index of the node this boundary edge belongs to. |
| `w_int` | `weight_int` (`uint32_t`) | The distance from that node to the boundary. Copied from `bcost_w_int`. |
| `det` | `uint64_t` | The detector id of that node. A later stage uses it to look up the boundary observable mask. |

### Optional output — counters

The software also fills a `struct BallGraphCounts`
(`src/specmatching/spec_matching/manifold/ball_graph.h:302`) with counts of how
many table bytes were touched and how many edges were written. These are
measurement aids, not part of the result. The hardware does not have to produce
them.

## The rules the output must obey

These rules are what make the output well defined. Two different implementations
that both follow them will produce byte-identical results, and that is the
property the rest of the system depends on.

1. **A pair appears at most once.** If detectors `a` and `b` both fired and are
   close enough, there is exactly one edge record for them, not two. The
   software achieves this by only emitting a pair when walking from the *lower*
   detector id — see the `if (target <= det) continue;` test at
   `src/specmatching/spec_matching/manifold/ball_graph.cc:221` and again at
   `:248`.

2. **A pair's data is always looked up from the lower detector id.** The
   distance `w_int` and the index `entry` must be taken from the lower-id
   endpoint's ball, not the higher-id one. This is not a matter of taste: the
   stored shortest paths are not symmetric under the tie-breaking rule used to
   pick them, so reading the pair from the other side can name a different path
   and therefore a different observable mask. Fixing the direction is what makes
   the edge a function of the pair alone. The comment at
   `src/specmatching/spec_matching/manifold/ball_graph.h:40` states this.

3. **The distance test is `w_int <= 2 * horizon`.** Strictly greater than that
   means no edge.

4. **The boundary test is `bcost_w_int[det] <= horizon`,** and only when
   `has_bcost[det]` is 1. Note the threshold here is `T`, not `2T`. See
   `src/specmatching/spec_matching/manifold/ball_graph.cc:271`.

5. **A node gets at most one boundary edge.** There is only one boundary
   distance stored per detector, so this follows automatically.

6. **A detector is never joined to itself.**

7. **The final edge list is sorted by `(i, j)`,** with `i` ascending as the
   primary key and `j` ascending as the tie-break. The software does this with an
   explicit sort at `src/specmatching/spec_matching/manifold/ball_graph.cc:315`.
   The boundary edge list is sorted by `i`. If your hardware emits edges in some
   other order, something downstream has to sort them, and you should say so
   clearly in your interface description, because the sortedness is relied on.

## The two ways the software already does it

The existing code contains two different implementations that must produce
identical output. They are selected by `enum class BallGraphBuildMode` at
`src/specmatching/spec_matching/manifold/ball_graph.h:28`. You should read both,
because they read *different* input tables, and which one you take as your
starting point decides which tables your hardware needs access to.

**`SCAN` mode** (`ball_graph.cc:212`): for each fired detector `v`, walk down
`v`'s ball entry list. For each entry, first check the distance; because the list
is sorted by distance, the first entry past `2T` ends the walk for that detector.
Then check whether the neighbour's id is greater than `v` (rule 1), and whether
the neighbour also fired. If both hold, emit an edge. This mode reads
`ball_offsets`, `ball_target`, `ball_w_int`.

**`BITSET` mode** (`ball_graph.cc:231`): for each fired detector `v`, take `v`'s
bitmap and AND it, word by word, against a bitmap of the whole shot's fired
detectors. Every set bit in the result is a neighbour that also fired. For each
such bit, count the set bits below it to get its position number, use
`ball_entry_by_rank` to turn that into an entry index, read the distance, and
apply the same two tests. This mode reads `ball_word_offsets`, `ball_word_base`,
`ball_words`, `ball_word_rank`, `ball_entry_by_rank`, `ball_w_int`.

The AND in `BITSET` mode needs a bitmap of the shot. The software builds one at
`src/specmatching/spec_matching/manifold/ball_graph.cc:191` by setting one bit
per fired detector in an array of `(num_nodes + 63) / 64` words, and clears it
again afterwards by remembering which words it touched
(`ball_graph.cc:320`). If your hardware takes `seeded_dets` as a list, it will
need to do the same conversion; if it can take the shot as a bitmap directly,
that step disappears.

One detail that trips people up in `BITSET` mode: a detector's bitmap window can
run off the end of the shot bitmap, and the loop stops there
(`ball_graph.cc:238`).

## A tiny worked example

Suppose `num_nodes = 100`, `horizon = 5`, so `2T = 10`.

Suppose detectors `7`, `12` and `40` fired. Then `seeded_dets = [7, 12, 40]` and
the node list is:

    h_to_det = [7, 12, 40]

so local index `0` means detector `7`, local index `1` means detector `12`, and
local index `2` means detector `40`.

Suppose detector `7`'s ball contains detector `12` at distance `6` and detector
`40` at distance `30`, and detector `12`'s ball contains detector `40` at
distance `24`. Suppose detector `7` has a boundary distance of `4` and the other
two have none.

Then the output is:

    edges          = [ {i=0, j=1, w_int=6, entry=<index of the (7,12) entry>} ]
    boundary_edges = [ {i=0, w_int=4, det=7} ]

Detector `40` produces nothing: `30` and `24` are both greater than `10`. The
`(7, 12)` edge is looked up from detector `7`, because `7 < 12`.

## What is not in scope

The following steps happen after this block and are not part of it. They are
mentioned only so the boundary of your work is clear.

- Grouping the resulting edges into connected components.
- Running the matching algorithm on those components.
- Turning matched pairs back into a correction.
- Reading the observable mask arrays.

Your block's job ends the moment the edge list and the boundary edge list exist.

## Checklist of what to produce

1. An interface specification: the exact signals or ports, their widths, and the
   handshake, for both input groups and all three output parts.
2. A statement of which of the two modes above your design follows, and
   therefore which input tables it must be able to read.
3. A statement of whether your output is already sorted by `(i, j)` or whether
   sorting is left to a later stage.
4. A statement of what your design does when a shot is larger than expected, or
   when the requested `horizon` exceeds `t_max_int`. The software throws an error
   in the second case (`ball_graph.cc:157`).
5. The HDL or HLS source itself, producing, for a given set of tables and a given
   shot, output identical to what `build_ball_graph` produces for the same
   inputs.
