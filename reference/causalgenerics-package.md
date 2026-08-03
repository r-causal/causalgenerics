# causalgenerics: Shared Generics for the 'r-causal' Ecosystem

A home for the S3 generics shared across the 'r-causal' ecosystem,
including 'propensity', 'halfmoon', 'positively', and 'balancing'.
Owning the generic definitions in one place lets those packages register
methods without masking one another when several are attached at once.
The generics cover inverse probability weighted estimation, effective
sample size, and the metadata carried by causal weight vectors, and the
package supplies the abstract causal weight class and the result class
those estimates are returned in, together with the methods both classes
carry, so that the ecosystem packages inherit them rather than writing
their own. The design follows that of the 'generics' package, which
provides commonly used S3 generics for the same purpose: letting
packages share a single definition instead of each defining its own.

## See also

Useful links:

- <https://github.com/r-causal/causalgenerics>

- <https://r-causal.github.io/causalgenerics/>

- Report bugs at <https://github.com/r-causal/causalgenerics/issues>

## Author

**Maintainer**: Malcolm Barrett <malcolmbarrett@gmail.com>
([ORCID](https://orcid.org/0000-0003-0299-5825)) \[copyright holder\]

Authors:

- Malcolm Barrett <malcolmbarrett@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-0299-5825)) \[copyright holder\]
