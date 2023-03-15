# ProbCast

ProbCast is an `R` Package under continuous development by researchers at the University of Strathclyde. It is a collection of functions for probabilistic forecasting (mainly wrappers for qunatile and semi-parametric regression model fitting functions), cross-validation, evaluation and visualisation. Central to ProbCast is the data class `MultiQR`, for storing the results of multiple quantile regression, and methods for working with `MultiQR` objects.

## Set-up

You can install the latest version of ProbCast using:

    # install.packages("devtools")
    library(devtools)
    install_github("jbrowell/ProbCast")
    
The latest release is:

[![DOI](https://zenodo.org/badge/143147931.svg)](https://zenodo.org/badge/latestdoi/143147931)
    
It may be necessary to set the following:

    Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")

## Usage

The package includes a script `Example.R` which demonstrates much of ProbCast's functionallity.

## Guide for Contributors
Contributors should follow the following guidelines:
- Open new branches when adding new functinoality, or making changes to exisig functions
- Raise *issues* when adding functionality, and identifying and fixing bugs
- Use rogygen to include documentation with each new function. Include helpful notes on all inputs, outputs, and wokrings of the function
- Include helpful comments throughout code
- Add yourself as @author to functions that you have "ownership" of
- Invite others to review pull requests (check for documentation, back compatability, confilcts...)

## Acknowledgements

Thanks to everyone who has contributed: Jethro Browell, Ciaran Gilbert, Gordon McFadzean, Rosemary Tawn.

Development of ProbCast has been supported by the following grants and organisations:
- EPSRC Innovation Fellowship "System-winde probabilistic energy forecasting" ([EP/R023484/1](https://gtr.ukri.org/projects?ref=EP/R023484/1) and [EP/R023484/2](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/R023484/2), 2018-2022)
- EPSRC Wind and Marine Energy Systems CDT ([EP/L016680/1](https://gtr.ukri.org/projects?ref=EP/L016680/1))
- Network Innovation Allowance project "Control REACT" ([NIA_NGSO0032](https://smarter.energynetworks.org/projects/nia_ngso0032))
- [University of Strathclyde](https://www.strath.ac.uk/)
- [University of Glasgow](https://www.gla.ac.uk)
- [TNEI Services Ltd](https://www.tneigroup.com/)

## References

ProbCast was introduced in the following paper, and has sinced been used in multiple academic studies and to facilitate training in probabilistic forecasting for researchers and practitioners.
- J. Browell and C. Gilbert, (2020), "ProbCast: Open-source production, evaluation and visualisation of probabilistic forecasts", DOI: [10.1109/PMAPS47429.2020.9183441](https://doi.org/10.1109/PMAPS47429.2020.9183441)
- [Citing articles](https://scholar.google.co.uk/scholar?oi=bibs&hl=en&cites=13378708372874427516)
- [Archived versions on Zenodo](https://doi.org/10.5281/zenodo.3843332)
