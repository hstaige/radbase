Format specifications of data input files. Last updated 12/24/24.

measurements.csv
=================

Description:

Contains measurements of absolute and relative radii measurements. Each row represents an individual measurement. Currently, all measurements are assumed to be independent.
Current supported types are absolute, linear relative, and squared relative (see "Terms" column below for details)

Columns:

1. 'MeasurementId'
    a. dtype: int
    b. Description: unique id of the measurement, primary key
2. 'Iso_Idxs'
    a. dtype: str
    b. Description: comma separated isotope indexes - an isotope index has the format ``f'{Z:03}_{A:03}'`` -
    in other words the Z and A of the isotope leftpadded with zeros to length three.
    EX: Carbon-12 (Z=6, A=12) has an isotope index of 006012
3. 'Term'
    a. dtype: str
    b. Description: String representation of the quantity measured by the experiment. The minimizer will optimize the parameters to make
    'Term' as close to 'Value' as possible.
    c. Supported formats of 'Term'
        1. Absolute - Example: "R006012" -> a measurement of :math: `R_{C12}`
        2. Linear relative - Example: "R006014-R006012" -> a measurement of :math: `R_{C14} - R_{C12}`
        3. Squared relative - Example: "R006014\**2-R006012\**2" -> a measurement of :math: `R_{C14}^2 - R_{C12}^2`
4. 'Value'
    a. dtype: float
    b. Description: Measured value of the quantity described by term. Units are assumed to be in fm or fm^2
5. 'Unc'
    a. dtype: float
    b. Description: Uncertainty of the measurement. Units are assumed to be the same as 'Value'
6. 'RefRemarks'
    a. dtype: str
    b. Description: reference string and other notes
    c. #TODO -> Seperate into ref column (comma seperated ref identifiers) and remarks
7. 'Table'
    a. dtype: str
    b. Description: source table of this data: one of (absolute, noniso, isom, isoe, kalpha, ois)


correlations.csv
================
Description: Contains correlation coefficients between values in measurements.csv

Columns:

1. 'MeasurementId1'
    a. dtype: int
    b. Description: Foreign key of the first measurement
2. 'MeasurementId2'
    a. dtype: int
    b. Description: Foreign key of the second measurement
3. 'Correlation':
    a. dtype: float
    b. Description: correlation coefficient between the first and second measurements; between -1 and 1, inclusive

references.csv
==============
Description: The references for measurements in **measurements.csv**

Columns:
