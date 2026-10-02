# Computing the Riemann Hypothesis: code and data

The scripts and data files behind the tables, figures and certificates of the book *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times* by Nicholas G. Polson and Vadim Sokolov (2026). Appendix A of the book, "Code and reproducibility", names every script, the data file it writes and the tables and figures it produces, and says which computations have no script. The folders here are the ones the appendix names.

## Layout

| Folder | Contents |
|---|---|
| `code/` | The scripts of Chapters 1 to 5, 8 and 12, among them the two Pólya evaluators `polya_flint.py` and `tilted.py` that Appendix A prints. |
| `data/` | The data files those scripts write, and `zeros_1700.txt`, the first 1700 ordinates of the zeros to seventy digits. |
| `verification_reciprocal_zeta/` | Chapters 6 and 7: interval enclosures, certificates and spot checks; `torus/` holds the branch and bound behind the torus bounds, with its own README. |
| `code_bosch/` | The verification script of Section 4.3; `bosch_hcm_xi/` holds an identical copy beside the authors' notes. |
| `basepoint_one/`, `sato_bondesson/`, `riesz_program/` | Chapters 9 to 11. Each folder keeps the layout of the authors' notes, and each script reads and writes its data files in its own folder. |
| `ferro/` | Chapter 13: the certificates of the ferromagnetic laws, each with its data file and the log of its run. |
| `misc/repro/` | The reruns of the whole archive on the Hopper cluster of George Mason University: the drivers `run_all.py`, `run6.py` and `run_riesz.py`, their SLURM files, the comparison reports, and in `reconstruct/` the scripts written for the book where the archive held a data file but no script. |

## Requirements

Python 3.11 with python-flint (the Python interface to FLINT and Arb, for ball arithmetic), mpmath, numpy, scipy, sympy and matplotlib. The reruns of September and October 2026 used the versions in `requirements.txt`: Python 3.11.4, python-flint 0.9.0, mpmath 1.3.0, numpy 2.4.6, scipy 1.17.1, sympy 1.14.0 and matplotlib 3.11.2.

## Running

Each script runs from its own folder and reads and writes its data files in the working directory; Appendix A gives the command lines and the order in which the scripts of `code/` run (`part1.py` and `part2.py` first). Several scripts read the first 1700 ordinates of the zeros from `g_1_400.npy`, `g_401_1000.npy` and `g_1001_1700.npy`, of which copies are in `riesz_program/zeros/` and in several other folders, and `misc/repro/zeros.py` computes them with `mpmath.zetazero`; others read `data/zeros_1700.txt`, which `code/zeros_hp.py` writes. The long computations ran on Hopper: the drivers in `misc/repro/` and the SLURM files name the paths of that cluster and need adapting elsewhere.

Where the book says that a statement is proved in ball arithmetic, the script named beside it is the certificate, and its JSON output records the enclosures. Examples are `code/gsmooth_arb.py` and `code/phasefree_arb.py` (Chapters 5 and 12), `basepoint_one/xi/heat_thresholds.py` (Table 9.7), and `ferro/fm_cutoff.py`, `fm_qcert.py`, `fm_tvcert.py` and `fm_intcert.py` (Chapter 13).

## Not included

* The LaTeX source of the book and of the authors' notes, and the notes compiled to PDF. The READMEs of the notes folders were written for the book's working repository and still name these files.
* `ferro/fm_zeros_hp.txt` (2.3 MB), which `python3 fm_qcert.py zeros 800 9600` writes in `ferro/`.
* The snapshot of six chapters of the book that `misc/repro/reconstruct/spot_checks/compare.py` reads; the quotations it compared are in `book_quotes.tsv` beside it.

## Citation

Nicholas G. Polson and Vadim Sokolov, *Computing the Riemann Hypothesis: Thorin Measures, Pólya Densities and Hitting Times*, 2026. Code and data: https://github.com/VadimSokolov/computing-rh-code

## License

MIT, see `LICENSE`.
