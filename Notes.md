# Chapter 2
## What is missing in Chapter 2

- Introduction is too long. 
- Empasize that gridding is  locations discreizatio **PLUS*** sampling/aggregating.
	- Exception: L2G, which does  **Nearest neigbor sampling** and keeps all values. 

## Comparison with STARE
- STARE: Binning is free at any choosen level
- We can pick a resolution and postpone the resampling with STARE. And we may choose to not do any resampling.
- We can store at different resolution and it is trivial to evaluate coincidence at varying resolutions.
- With gridded data we are stuck with a fix resolution. This means we either have to do sampling or have gaps in our grid (EarthDB)

## Sidenotes
- Analytical frameworks are predicated on projections. This is all fine, except that values have to be sampled. Which get's tricky when the original values vary in size.
- How do we represent iFOVs: Centroid -> Then we can also represent ellipses (Spoiler C3)


# Chapter 3
- What should the stare level for an iFOV be?  The size it represents or the geolocation precision
- Drifting orbits. We will get all continous viewing angles. This will effect e.g. R0 
- How many natural bins are there for a given latitude?