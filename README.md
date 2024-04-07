# MicrolensingAnimation
Animation of Microlensing Event Model

The program implements [Microlensing Event](https://science.nasa.gov/mission/roman-space-telescope/microlensing/)

## A little theory
Microlensing - the lensing effects of small-sized lenses in the universe. Such effects are mostly detectable within our galaxy, so we will focus on microlenses in the Milky Way.

Prameters: 
- **$G$** - the gravitational constant
- **$c$** - the speed of light
- **$D_L$** - angular diameter distance to the lens
- **$D_S$** - angular diameter distance to the source
- **$D_{LS}$** - angular diameter distance between the lens and the source
- **$\beta$** - intrinsic position of the source
- **$\theta$** - observed image position of the source

The Einstein radius: $\theta_E \equiv \sqrt{\frac{4GM_L}{c^2}\frac{D_{LS}}{D_LD_S}}$

The lens equation: $\beta = \theta - \frac{4GM_L}{c^2D_L\theta}\frac{D_{LS}}{D_S}$
we may see, that for a given source position $\beta$ there will always be two images.
