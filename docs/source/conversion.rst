Ramblings on the conversion process and design.

The general problem is we have a list of potentially correlated input quantities
that we would like to convert to radius data. My number one goal with this system is flexibility;
I should never feel like I am designing my analysis around the system. I also want traceability in
the conversion process so that information can be displayed in dissemination work.

For the latter requirement, it seems like I should be keeping track of the conversion steps - something like
'Measurement.conversion_steps'. These conversion steps should keep track of the full story of how we got this radius.
One example set of conversion information could be

Absolute Radius of C-12, 3.14 fm +/- 0.01 fm
1. Radius calculated using Rka and V2
    a. Rka determined using theortical calculations of refA and experimental energies of refB
        1. theory from refA, includes nuclear polarization corrections, local only.
        2. Exp energy 1 = 3000 keV, Exp energy 2 = 4000 keV
    b. V2 factor calculated from Fourier Bessel coefficients using direct integration
        1. Fourier bessel coefficients and covariance taken from reference C.

Some building blocks of this process are:

Conversion step format - what should be included and how should the
conversion and documentation of that conversion be coupled. Perhaps this is a fairly easy question where we have a conversion step class
which has a convert method and a conversion_doc method that returns a string describing the process?

Some conversion steps are built in to the system - for example, you can convert fermi distribution parameters to a
barrett radius, you can convert FB coefficients to the V2 factor, etc. Others are provided by the theory. For example,
muonic energies to barrett radii, the sensitivity coefficients in the process, etc.

The converter class will take in all the input quantities and the conversion steps and handle the graph that results.
This looks something like starting with the measured muonic energies, which can be turned into a barrett radius quantity,
which can potentially be turned into an rms radius if V2 factors are available. Alternatively, you can just take the Fermi
distribution parameters and calculate a radius directly

In terms of preferences, we would like to have at a minimum, the following options.

 1. Select specific measurements or exclude specific measurements. We should be able to create flexible filters for included measurements.

 2. Pre-processing steps - the big one I am thinking of is the fit over V2 factors / interpolation of missing V2 factors.
  a.This would look something like at the start of the conversion process we get all the V2 factors, do a fit, and then
whenever a V2 factor is requested we return the value of the fit rather than the measured value.

 3. Select preferred conversion paths - for example, if you can, always use the barrett radius recipe rather than a direct
Fermi distribution.
