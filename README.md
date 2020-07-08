# A Collection of Algorithms for Relative Pose Estimation of a Calibrated Camera

[![BSD-3-Clause](https://img.shields.io/github/license/prclibo/relative_pose)](https://github.com/prclibo/relative_pose/blob/master/LICENSE)
[![ECCV 2020](https://img.shields.io/badge/ECCV-2020-%231b75bc)]()

## Introduction

This repository contains the following relative pose estimation solvers, in C++ and Matlab API:

* [x] The conventional 5-point algorithm (5P), wrapped from Hartley's well-known implementation.
* [x] 4-point algorithm with a known rotation angle (4P-RA)
* [x] 4-point algorithm under planar motion without knowing the plane direction (4P-ST0)
* [x] 3-point algorithm with a known rotation angle and under planar motion without knowing the plane direction (3P-RA-ST0).

Compared to other relative pose estimation algorithms, 4P-RA, 4P-ST0 and 3P-RA-ST0 leverage extra sensor/motion constraits without requiring extrinsics calibration. This is due to the interesting property of SE(3) invariants.

## Reference

This repository is the source code for the following paper:

```bibtex
@article{li2020relative,
  title={Relative Pose Estimation of Calibrated Cameras with Known SE(3) Invariants},
  author={Li, Bo and Martyushev, Evgeniy and Lee, Gim Hee},
  journal={ECCV},
  year={2020}
}
```

## Contents

* `form` Source code to derive different solver formulations compared in Table 3 in the paper. [Automatic Generator](https://github.com/PavelTrutman/Automatic-Generator) and [Generator for Automatic Polynomial Solvers](https://github.com/prclibo/gaps) are required to run these code.
* `include` C++ API.
* `matlab` Matlab API.
* `perf` Source code for the experiment section in the paper.
