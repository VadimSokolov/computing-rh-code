# Torus bounds for the truncations (Theorem rz:thm:torus)

The programs in this folder prove the torus bounds of Section ch:rz4 (ch/rz4.tex): the bounds for the torus extremes in Tables rz:tab:torus and rz:tab:torus1, the thresholds printed there, and the second condition of Theorem rz:thm:torus (Remark rz:rem:torus). Corollary rz:cor:Storus uses them above the thresholds; below the thresholds it uses floating point argument principle counts, which are not in this folder.

## What is proved

Notation as in ch/rz4.tex: P_n(x, w) = sum_k beta_k k^(-2x) exp(-i sum_p v_p(k) w_p) on the torus of the primes p <= n, r(x) = sum_k |alpha_k| k^(-2x), r'(x) = sum_k 2 log k |alpha_k| k^(-2x), m_n = min |P_n|, M_n = max -Re(P_n'/P_n) and l_n = max |P_n'/P_n| (the book's \ell_n), the extremes taken over the band and the whole torus.

- Band 1.2 <= x <= 1.3, n = 16, ..., 24, 32, 36: the lower bounds for m_n and the upper bounds for M_n and l_n of Table rz:tab:torus.
- Band 1 <= x <= 3/2 (the single band variant), n = 7, 9, ..., 15: the same bounds for Table rz:tab:torus1.
- The second condition of Theorem rz:thm:torus, B_n(y) < log(y/pi) - 1/y together with |P_n(1+s, w)| > r(1+s)/y for 0.3 <= s <= 1/2, at the heights y = 29, 33, 36, 40, 45, 49, 55, 61, 68, 211, 555 for n = 16, ..., 24, 32, 36. It improves as y increases, so it holds at every larger height, in particular from each threshold on.
- The thresholds of both tables: the smallest integer height Y at which the first condition holds with the proved bounds, with g(Y) > 0 > g(Y-1) decided in ball arithmetic, both with the exact values of r and r' at the lower end of the band and with the values rounded up to three decimals that the tables print.

The bounds are in certified_bounds.json and the thresholds in arb_thresholds.json.

## Provenance

The files come unchanged from the numerical report N3 of the review triage (misc/gpt-review/triage/numerics/N3.md, whose sections "The branch and bound", "Thresholds" and "The second condition" give the method in full; the programs were in misc/gpt-review/triage/numerics/N3/). Every computation ran on Intel nodes of the Hopper cluster with `python-flint` 0.9 (Arb), numpy and scipy, in the Python environment /scratch/vsokolov/rh_book_repro/venv, except that N3 first ran make_tasks.py, which only writes the task list (the pieces, and the targets and heights rounded from the estimates), on a laptop; its rerun on Hopper (job 1313099) gave the same tasks_full.json byte for byte.

The folder holds the programs and the aggregated outputs. The task lists and the per task results (result_*.json, items_*.json and the frontier files of the replay) are not included: the task lists are rebuilt by the steps below, and the per task results of our run stay in the run folders /scratch/vsokolov/rh_book_repro/agents/N3-* on Hopper. Also left out are an abandoned search for the second condition in logarithmic form and development files, which no result uses.

## Method

Functions. The coefficients alpha_k and beta_k are exact rationals (torus_common.py). Since beta_k is real, P_n(x, -w) is the conjugate of P_n(x, w), which leaves every quantity below unchanged, so the half torus 0 <= w_2 <= pi suffices.

Enclosure of a box. On a box with centre (x_c, w_c) and half widths h in x and h_p in w_p, each term moves by at most beta_k k^(-2 x_c) ((k^(2h) - 1) + k^(2h) min(2, sum_p v_p(k) h_p)). Summing gives balls for P_n, P_n' and their partial derivatives, and from them a lower bound for |P_n|^2 and upper bounds for -Re(P_n'/P_n) and |P_n'/P_n|^2 on the box, by a mean value form or a first order form, whichever is sharper. For M_n and l_n the proved lower bound for m_n on the same band serves as the lower bound for |P_n|. The second condition at height y is proved as G > 0 on 0.3 <= s <= 1/2 and the torus, where G(s, w) = a_2 exp((2s - 1/2) T) - a_1 with T = log(y/pi) - 1/y, a_1 = |P_n(3/2 - s, w)| + r(3/2 - s)/y and a_2 = |P_n(1 + s, w)| - r(1 + s)/y; since 2s - 1/2 > 0 and a_1 > 0, this is equivalent to the second condition together with a_2 > 0, and it needs no division.

Search. A depth first branch and bound in numpy double precision discards a box when its bound clears the target by 1e-9, and otherwise halves it along the coordinate with the largest term of the mean value form. The targets are the three decimal outward roundings of multistart estimates (m down, M and l up), and 0 for G. The domain is cut into pieces: one per group for m, 64 for M and for l (w_2, w_3, w_5 in quarters), 16 per n for G (w_2, w_3 in quarters). A box centre on the wrong side of a target would stop the run as a counterexample; none occurred.

Proof in ball arithmetic. The double precision search only chooses the boxes. Every final box is proved again in Arb at 128 bits, with the same formulas, exact rational coefficients and every rounding error enclosed, on the exact box: each final box is identified with an exact dyadic part of the exact piece (x bounds 1, 6/5, 13/10, 3/2; s bounds 3/10, 1/2; phase bounds rational multiples of pi). For m, M and l this is a replay: each piece is cut into sub boxes of at most about 500000 final boxes, the search is rerun on each sub box keeping every final box, and every final box is proved in Arb. For G the run itself keeps and proves every final box. For every sub box the recovered dyadic positions are checked to be consistent and distinct, and the exact volumes of the final boxes to add up to the volume of the sub box, in integer arithmetic; the sub boxes of each piece are checked in the same way to cover the piece. A box failing in Arb would be bisected exactly in Arb; none failed.

Thresholds. arb_thresholds.py encloses r and r' at the lower end of the band as exact finite sums, brackets the root of the increasing function g(y) = log(y/pi) - 2/y - M - (l r + r' + r/y)/(y (m - r/y)) on y > r/m by bisection with sign decisions in Arb, and certifies g(Y) > 0 > g(Y-1) for the threshold Y.

Totals. For m, M and l: 57 groups, 2451 pieces, 6412 sub boxes and 106085082 final boxes, all proved in Arb in 71015 CPU seconds (the double precision searches took 924 CPU seconds). For G: 176 pieces and 2090726 final boxes, all proved in Arb in 1624 CPU seconds (search 27 seconds). No box failed, every coverage check passed, and no piece is missing.

Rigour. The results are proved assuming that the ball arithmetic of `python-flint` and our implementation of the enclosure formulas are correct. The second assumption is the real risk, and two checks support it. (a) Before the replay, in the runs for M and l (bnb4.py), a sample of the final boxes of every piece, for every n (up to 40 of smallest margin and 60 random ones per piece; 139297 checks in all, a box kept both for its margin and at random being checked twice), was proved in Arb with the same formulas (bound_ratio_mlow in arb_disc.py), and every check passed; the largest difference between the double precision bound and the Arb bound was 1.6e-11, against the margin 1e-9. The same boxes went to the independent Arb enclosure check2_mlow of arb_leafcheck.py (interval evaluation of the partial derivatives, with the same lower bound for |P_n|, its own bisection and a budget of 1500 boxes), which succeeded in 135708 of the 139297 checks and ran out of budget in the others (full_status.json). (b) A separate run (bnb3.py) has 2695 tasks, a 240 second search of every piece and a 900 second search of a small box around each extreme point (make_verify.py). N3's array was cancelled after 1831 tasks had finished; the other 864 ran on 28 September 2026 (Hopper job 1316484, from a list of the missing task numbers), so the run now covers every piece of m, M, l and of the abandoned quotient search q, for every n, and the box around every extreme point. For M and l this run uses the enclosure of bnb.py, which bounds |P_n| on each box itself, and not the enclosure of bnb_m.py with the proved lower bound for |P_n| that gives the results; that step is checked by (a) only. Within its time limit the search of a piece either ends below the target or gives up, and it gives up on most pieces of M and l for large n, so (b) checks the enclosures on the boxes it reaches and does not prove the bounds by itself. The final boxes kept by these tasks, chosen as in (a) (242009 checks), were checked three ways: with the Arb version of the enclosure of the search (no failure; largest difference 1.1e-11, relative 1.0e-13), with an independent Arb enclosure by interval evaluation of the partial derivatives (arb_leafcheck.check2, and arb_disc.ArbQNaive for q), with its own bisection and a budget of 1500 boxes (the 4264 checks it could not decide within the budget are counted as naive_fail), and at 300 random points per box in double precision (no point on the wrong side of a target) (verify_aggregate.json). Tasks 0 to 863 also ran a second time, by mistake, in job 1315369: their checks agree with the first run, with no failure in either, and the searches differ only where the time limit stopped them at different points. A rerun of step 6 runs all 2695 tasks. The bounds also agree with values at points: every bound for M_n is within 0.0012 of the value at w = 0 at the lower end of the band, x = 6/5 or x = 1 (arb_points.json).

## Files

Common:

- torus_common.py: the exact coefficients alpha_k and beta_k, the primes and valuations, and the torus functions in double precision.

Estimates and pieces:

- explore.py: multistart optimisation (300000 random samples, then L-BFGS-B from the best 400 samples and from 1500 random starts, with fixed seeds) of m, M and l for every n and band; writes explore.json. Floating point; used only to place the targets.
- thresholds.py: floating point thresholds from the estimates under several readings, the estimated start of the second condition, and the unrounded root y1 of the first condition; writes thresholds.json.
- make_tasks.py: the pieces of the branch and bound, with the targets (three decimal outward roundings of the estimates) and the heights floor(y1) of the second condition; writes tasks_full.json.

Enclosures:

- bnb.py: the double precision enclosures for m, M and l (class Enc, margin SLACK = 1e-9).
- bnb_m.py: the enclosures for M and l with the proved lower bound m_low of |P_n| (class EncM).
- bnb_q.py: the enclosure of the quotient of B_n in logarithmic form (class EncQ). Its search was abandoned, but bnb_g.py builds on its balls, and bnb2.py and bnb3.py import it.
- bnb_g.py: the enclosure of G (class EncG).
- arb_disc.py: the Arb versions, at 128 bits, of the enclosures of bnb.py, bnb_m.py and bnb_q.py.
- arb_g.py: the Arb version of the enclosure of G.
- arb_leafcheck.py: an independent Arb enclosure by interval evaluation of the partial derivatives, with its own bisection, for the checks.

Drivers (one array task each):

- bnb2.py: the search of one piece (run N3-full, of which only the pieces of kind m are used), with an Arb check of a sample of its final boxes.
- bnb4.py: the same for M and l with m_low (run N3-mlow).
- bnb8.py: the search for G, with the Arb proof of every final box and the coverage checks (run N3-g).
- bnb7.py: the replay of m, M and l: reruns the search on each sub box of a task keeping every final box, proves each in Arb and checks the coverage (runs N3-replay, N3-replay2, N3-replay3, N3-replay4).
- bnb3.py: the leaf verification of check (b) (run N3-verify).

Task lists of the replay:

- make_replay.py: cuts every piece of m (from N3-full) and of M and l (from N3-mlow) into a dyadic grid of sub boxes of about 500000 final boxes or fewer, estimated from the final boxes of the piece, and packs them into tasks; writes tasks_replay.json.
- make_replay3.py: for the pieces 2298 (n = 32, l) and 2426 (n = 36, l) of tasks_mlow.json, whose grid sub boxes needed too much memory, sub boxes along the top of the search tree itself: a box is split by the rule of the search until the counted search from it has at most 400000 final boxes; writes items_TASK.json.
- make_replay3b.py: the same for one piece, in parallel. "frontier INDEX" expands the top of the tree breadth first (frontier.json, items_cleared.json); an array task counts below one frontier box (items_TASK.json); "refine FILE SIZE TID..." expands the frontier boxes TID further, for boxes too heavy for one task (frontier2.json, items_cleared2.json), and an array over frontier2_tasks.json counts below them.
- pack_replay3.py: collects the sub boxes of a piece, checks that their exact volumes add up to the piece, and packs them into tasks of at most 500000 counted final boxes (tasks_replay3.json for 2298, tasks_replay4.json for 2426).

Aggregation:

- make_bounds.py: the status of every group in the runs N3-full (m), N3-mlow (M, l) and N3-g (G); writes full_status.json and certified_bounds.json.
- aggregate_replay.py: the replay per group, with the check that the proved sub boxes cover every piece; the sub boxes of the pieces 2298 and 2426 along their search tree replace their grid sub boxes; writes replay_summary.json.
- aggregate_verify.py: the leaf verification per group; writes verify_aggregate.json.

Values at points and thresholds:

- arb_points.py: Arb values of |P_n|, -Re(P_n'/P_n) and |P_n'/P_n| at the optimiser points of explore.json and at the points w = 0 and w_p = pi for all p. A value at a point bounds an extreme from the favourable side only, so it can refute a bound but not prove one; writes arb_points.json.
- arb_thresholds.py: the thresholds in Arb from certified_bounds.json, with the check at the printed three decimal values of r and r' and, from the values at points, a lower bound for the true threshold; writes arb_thresholds.json.

Checks:

- make_verify.py: the tasks of the leaf verification (every piece again for 240 seconds, and a small box around each extreme point for 900 seconds); writes tasks_verify.json. The run itself is bnb3.py, and aggregate_verify.py collects it.

Submission:

- submit_array.sh: submits an sbatch array on the Intel nodes of Hopper, one task per entry of a task list, in /scratch/vsokolov/rh_book_repro/agents/N3-NAME; waits for it, then copies result_*.json and the logs back to ./out_NAME/.

Outputs:

- explore.json: the multistart estimates with their optimising points (estimates, not bounds).
- thresholds.json: the floating point thresholds and heights (estimates).
- arb_points.json: the values at points in Arb.
- full_status.json: per group (n, kind, band or height): the target, the pieces certified or proved, the final boxes, CPU seconds, the Arb sample of check (a), and for G the smallest value at a box centre.
- certified_bounds.json: the proved bounds m, M and l per n and band and, for the band [1.2, 1.3], the height and target of the second condition; the input of arb_thresholds.py.
- replay_summary.json: per group, the pieces, sub boxes and final boxes proved in Arb, failures, CPU seconds and coverage, with the totals.
- arb_thresholds.json: per n and band, r and r' at the lower end of the band, the root and threshold with the proved bounds (certified_3dp), the same with r and r' rounded up to three decimals (table_check), the lower bound from values at points, variants with two decimals, the values printed in the earlier version of Table rz:tab:torus, the heights of the argument principle counts, and the height of the second condition.
- verify_aggregate.json: the checks of (b) per group.

## Rerunning on Hopper

Requirements: ssh access to Hopper under the host alias hopper, jq on the local machine, and the Python environment above. The paths /scratch/vsokolov/rh_book_repro/agents (D in make_replay.py, pack_replay3.py, aggregate_replay.py, make_bounds.py and aggregate_verify.py, FULL in make_verify.py, dir in submit_array.sh) and /scratch/vsokolov/rh_book_repro/venv are written into the scripts; another user changes them there. The aggregation scripts read the per task results from the run folders N3-full, N3-mlow, N3-g, N3-replay, N3-replay2, N3-replay3, N3-replay4, N3-mr3, N3-mr3b-front, N3-mr3b, N3-mr3c-front, N3-mr3c and N3-verify, so a rerun keeps these names.

Single jobs run with the book's misc/tools/hopper_run.sh (options -n run folder, -c cores, -m memory, -t minutes, -g file to fetch; then the script, its input files, and after -- its arguments). Arrays run with submit_array.sh NAME SCRIPT TASKS.json MINUTES MEMORY [MAXIMUM PARALLEL TASKS] -- FILES. The minutes and memory below are those of our run. Every command runs in this folder:

```
H=../../misc/tools/hopper_run.sh
```

1. Estimates and pieces. The two optimisations use fixed seeds; to reproduce our numbers exactly, skip them and start from the archived explore.json and thresholds.json, from which make_tasks.py rebuilds tasks_full.json byte for byte, and the two jq lines rebuild tasks_mlow.json (M and l, with m_low the target of m for the same n and band) and tasks_g.json (G, one task per piece of kind q) byte for byte.

```
bash $H -n N3-explore -c 48 -m 48G -t 180 -g explore.json explore.py torus_common.py -- 48
bash $H -n N3-thresholds -c 24 -m 24G -t 120 -g thresholds.json thresholds.py torus_common.py explore.json -- 24
bash $H -n N3-tasks -c 1 -m 2G -t 10 -g tasks_full.json make_tasks.py torus_common.py explore.json thresholds.json
jq '(map(select(.kind == "m") | {key: "\(.n) \(.band)", value: .target}) | from_entries) as $m | map(select(.kind == "M" or .kind == "l") | .time_limit = 3000 | .m_low = $m["\(.n) \(.band)"] | .label = "mlow piece \(.piece)")' tasks_full.json > tasks_mlow.json
jq -c '[to_entries[] | select(.value.kind == "q") | {q_index: .key, cuts: [], time_limit: 3000, arb: true}]' tasks_full.json > tasks_g.json
```

2. Searches in double precision, and the status of every group. The run full searches every piece of tasks_full.json, but only its pieces of kind m are used (the runs mlow and g supersede the others, some of which do not finish in the time limit).

```
bash submit_array.sh full bnb2.py tasks_full.json 110 6G 500 -- *.py
bash submit_array.sh mlow bnb4.py tasks_mlow.json 60 4G 500 -- *.py
bash submit_array.sh g bnb8.py tasks_g.json 70 8G -- tasks_full.json *.py
bash $H -n N3-makebounds -c 1 -m 8G -t 15 -g certified_bounds.json -g full_status.json make_bounds.py
```

3. The replay in Arb of m, M and l on the grid sub boxes. tasks_replay2.json lists the tasks of the run replay that left no result because they ran out of memory (136 of 323 in our run); they run again with more memory.

```
bash $H -n N3-replay-prep -c 1 -m 8G -t 20 -g tasks_replay.json make_replay.py
bash submit_array.sh replay bnb7.py tasks_replay.json 60 4G -- tasks_full.json tasks_mlow.json *.py
jq -c --slurpfile have <(ls out_replay | sed -n 's/^result_\([0-9]*\)\.json$/\1/p') '[to_entries[] | select(.key as $k | $have | any(. == $k) | not) | .value]' tasks_replay.json > tasks_replay2.json
bash submit_array.sh replay2 bnb7.py tasks_replay2.json 90 12G -- tasks_full.json tasks_mlow.json *.py
```

4. The replay of the pieces 2298 and 2426 of tasks_mlow.json along their own search tree. make_replay3.py builds the sub boxes serially, task 0 for 2298 and task 1 for 2426 (tasks_mr3.json only sets the number of tasks). Every route splits by the rule of the search, so the final boxes are the same whatever the sub boxes (15447067 and 75784949); pack_replay3.py uses the serial result for 2426 when N3-mr3/items_1.json exists, and otherwise the parallel route of make_replay3b.py. In our run the serial tasks took 863 and 3255 seconds, but the piece 2426 was packed while its serial task was still running, so tasks_replay4.json came from the parallel route: a frontier of 302 boxes, of which the two slowest, 124 and 128, were refined into 2550 boxes, counted in 255 tasks of ten. For the serial route, run the second line with tasks_mr3.json set to [2298, 2426] and skip the six lines that follow it. The commands below take the route of our run; pack_replay3.py then requires every frontier box whose task in N3-mr3b left no items file to be among the refined ones.

```
echo '[2298]' > tasks_mr3.json
bash submit_array.sh mr3 make_replay3.py tasks_mr3.json 90 6G -- tasks_mlow.json *.py
bash $H -n N3-mr3b-front -c 1 -m 4G -t 20 -g frontier.json make_replay3b.py tasks_mlow.json *.py -- frontier 2426
jq -c .boxes frontier.json > frontier_tasks.json
bash submit_array.sh mr3b make_replay3b.py frontier_tasks.json 60 4G -- frontier.json tasks_mlow.json *.py
bash $H -n N3-mr3c-front -c 1 -m 4G -t 20 -g frontier2.json make_replay3b.py frontier.json tasks_mlow.json *.py -- refine frontier.json 1024 124 128
jq -c '.boxes | length as $n | [range(0; $n; 10) as $i | [range($i; [$i + 10, $n] | min)]]' frontier2.json > frontier2_tasks.json
bash submit_array.sh mr3c make_replay3b.py frontier2_tasks.json 60 4G -- frontier2.json tasks_mlow.json *.py
bash $H -n N3-pack3 -c 1 -m 2G -t 10 -g tasks_replay3.json pack_replay3.py -- 2298 tasks_replay3.json
bash $H -n N3-pack4 -c 1 -m 4G -t 10 -g tasks_replay4.json pack_replay3.py -- 2426 tasks_replay4.json
bash submit_array.sh replay3 bnb7.py tasks_replay3.json 60 4G -- tasks_full.json tasks_mlow.json *.py
bash submit_array.sh replay4 bnb7.py tasks_replay4.json 60 4G -- tasks_full.json tasks_mlow.json *.py
bash $H -n N3-aggreplay -c 1 -m 8G -t 20 -g replay_summary.json aggregate_replay.py
```

5. Values at points and thresholds.

```
bash $H -n N3-arbpoints -c 1 -m 4G -t 20 -g arb_points.json arb_points.py torus_common.py explore.json
bash $H -n N3-thresholds-final -c 1 -m 4G -t 30 -g arb_thresholds.json arb_thresholds.py torus_common.py arb_points.json certified_bounds.json
```

6. The leaf verification of check (b).

```
bash $H -n N3-mkverify -c 1 -m 4G -t 20 -g tasks_verify.json make_verify.py tasks_full.json explore.json
bash submit_array.sh verify bnb3.py tasks_verify.json 50 4G -- *.py
bash $H -n N3-aggverify -c 2 -m 16G -t 30 -g verify_aggregate.json aggregate_verify.py
```

A rerun succeeds when replay_summary.json reports every group proved with no failure and no missing piece, every group of kind G in full_status.json is proved, and arb_thresholds.json gives the thresholds of the tables, with table_check equal to certified_3dp.
