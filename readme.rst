radbase
=======
.. image:: https://codecov.io/github/hstaige/radbase/graph/badge.svg?token=4S0KFLSOZG
  :target: https://codecov.io/github/hstaige/radbase

Overview
---------

Created in preparation for the 2025 IAEA Technical Meeting on Compilation and Evaluation of Nuclear Charge Radii.

**This repostiory is in an alpha stage; all features/data/syntax are subject to change without notice.**

Contains the following classes:


* ``RadiusAnalyzer`` - Flexible analysis via a non-linear least squares procedure (in-progress)

  * Built-in statistical outlier detection (not started)
  * Analyzes graph formed by data to reduce computational load (In progress)

* ``RadiusDataVisualizer`` - Interactive visualization of measurements and optimized radii (in-progress)

  * Aids in visualizing:

    * Primary vs. Secondary nuclides
    * Number of measurements for each connection/nuclide
    * Coverage of individual labs

* ``DataEntryVisualizer`` - User interface to enter nuclear radius data (in-progress)

  * Includes different templates for different kinds of information
  * Fields include validation, dynamic resizing, and searches of previously entered data
