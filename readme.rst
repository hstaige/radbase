radbase
=======
.. image:: https://codecov.io/github/hstaige/radbase/graph/badge.svg?token=4S0KFLSOZG
  :target: https://codecov.io/github/hstaige/radbase

Overview
---------

**This repository is in an alpha stage; all features/data/syntax are subject to change without notice.**

As described in the 2026 white paper "Towards better nuclear charge radii", there is currently a need for systematic
tools to store and process data related to nuclear charge radii. ``radbase`` is currently being developed to address this need, with the following long-term objectives:

- Collection of all nuclear charge radius data, including experimental quantities and calculated corrections.
- Development of an automated correlation propagation procedure that converts correlations among input data into correlations between nuclear charge radii.
- Development of modules to search compiled data using user-defined queries.
- Optimization routines that optimize radii given several thousand radius data measurements.

The short-term objectives for ``radbase`` are to:
- Compile all absolute nuclear charge radius data for $Z>2$.
- Finalize the data storage format to support symbolic manipulation and explicit conversion steps.

We expect to complete this stage of the project by the start of Fall 2026. Next, the relative radius data
(experimental/theoretical inputs for radius differences) will be compiled for $Z>2$. At the same time, the automatic
correlation propagation procedure will be developed. Once the project is in a stable state, documentation of the data
and features to support use by those outside the``radbase`` development team will be created.

For questions, please email Hunter Staiger at hstaige@clemson.edu
