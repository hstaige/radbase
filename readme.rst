radbase
=======
.. image:: https://codecov.io/github/hstaige/radbase/graph/badge.svg?token=4S0KFLSOZG
  :target: https://codecov.io/github/hstaige/radbase

Overview
---------

Created in preparation for the 2024 IAEA Technical Meeting on Compilation and Evaluation of Nuclear Charge Radii.

**This repostiory is in an alpha stage; all features/syntax are subject to change without notice.**

Contains the following classes:

* ``RadiusDataMapper`` - A database of all nuclear radius related measurements and associated references (in-progress).
* ``RadiusAnalyzer`` - Flexible analysis via a non-linear least squares procedure (in-progress)

  * Built-in statistical outlier detection (not started)
  * Analyzes graph formed by data to reduce computational load (In progress)

* ``RadiusDataVisualizer`` - Interactive visualization of measurements and optimized radii (not started)

  * Aids in visualizing:

    * Primary vs. Secondary nuclides
    * Number of measurements for each connection/nuclide
    * Coverage of individual labs
