# Invariant SE3 Path Descriptor

This repository contains MATLAB and Python scripts to implement DHB invariant representation of rigid body motions. The original code was written by [Matteo Saveriano](https://github.com/matteosaveriano/DHB_invariant_representation). The original code was written in MATLAB and this repository contains a Python implementation of the same code.

## Components description
The folder `matlab` contains the following constraints:
- ```computeDHB.m```: Compute DHB invariants given a Cartesian trajectory (6D pose or twist).
- ```reconstructTrajectory.m```: Reconstruct a Cartesian trajectory from its DHB invariant representation.
- ```mainDHB.m```: Simple demo script that shows how a simple Cartesian 6D trajectory can be transformed into a set of DHB invariants and vice versa.
- ```tutorial_DHB.mlx```: MATLAB live script that demonstrates some use-cases of the DHB invariant representation.
- ```test_DHB_invariants.m```: Test the DHB invariant representation in various scenarios.

## TODO

- Validate and demonstrate the invariance properties that were claimed in the original paper.
- Show the generalization capability with the DHB invariant representation.
- Create Python scripts and notebook for examples.
- Turn the script into a Python package with pip.

## Software Requirements
The code is developed and tested under _Matlab 2023b_.

## References
Please acknowledge the original authors in any academic publication that used parts of these codes.
```
@article{Lee2018,
    title = "Bidirectional invariant representation of rigid body motions and its application to gesture recognition and reproduction",
    author = "Lee, Dongheui and Soloperto, Raffaele and Saveriano, Matteo",
    journal = "Autonomous Robots",
    year = "2018",
    volume = "42",
    number = "1",
    pages = "125--145"
}

@article{Soloperto2015,
    title = " A Bidirectional Invariant Representation of Motion for Gesture Recognition and Reproduction",
    author = "Soloperto, Raffaele and Saveriano, Matteo and Lee, Dongheui",
    journal = " International Conference on Robotics and Automation (ICRA)",
    year = "2015",
    pages = "6146-6152"
}
```

## Maintainer

- Andy Park <andypark.purdue@gmail.com>
- Sungjoon Choi <sungjoon-choi@korea.ac.kr>

## Licence
This repository contains free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.

This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this code. If not, see http://www.gnu.org/licenses/.
