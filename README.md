# MicrolensingAnimation
Animation of Microlensing Event Model

The program implements [Microlensing Event](https://science.nasa.gov/mission/roman-space-telescope/microlensing/), i.e. the lensing effect of small-sized lenses in the universe. Such effects are mostly detectable within our galaxy, so we will focus on microlenses in the Milky Way.

## Necessary equations

Parameters: 
- **$G$** - the gravitational constant
- **$c$** - the speed of light
- **$D_L$** - angular diameter distance to the lens
- **$D_S$** - angular diameter distance to the source
- **$D_{LS}$** - angular diameter distance between the lens and the source ($= D_L - D_S$ in the Milky Way)
- **$\beta$** - intrinsic position of the source
- **$\theta$** - observed image position of the source
- **$u_0$** - impact parameter (angular distance between the lens and the source, normalized by the Einstein radius)
- **$M_L$** - mass of the lens
- **$t_0$** - time of magnification peak
- **$v$** - source relative velocity

The Einstein radius: $\theta_E \equiv \sqrt{\frac{4GM_L}{c^2}\frac{D_{LS}}{D_LD_S}}$

The Einstein radius crossing time (timescale of the microlensing event): $t_E = \frac{D_L\theta_E}{v}$

The lens equation: $\beta = \theta - \frac{4GM_L}{c^2D_L\theta}\frac{D_{LS}}{D_S}$ 

we may see, that for a given source position $\beta$ there will always be two images. In dimensionless form: $x_{\pm} = \frac{1}{2} (1 \pm \frac{\sqrt{y^2 + 4}}{y}) y$

Using the Einstein radius we will get: $\beta = \theta - \frac{\theta_E^2}{\theta}$

Dimensionless distance of the source: $y(t) = \sqrt{u_0^2 + \frac{(t - t_0)^2}{t_E^2}}$

Magnification of the source flux: $A(t) = \frac{y(t)^2 + 2}{y(t)\sqrt{y(t)^2 + 4}}$

## Files
- **INPUT.json** - file with input parameters
- **Microlensing.py** - the main program
- **Microlensing_lc_example.png** - example plot for microlensing light curve and positions of the images (dark circle - the Einstein ring, dashed lines - trajectories of the source and images, yellow star - the lens, colored circles mark the position of the images at a specific time)
- **Animation_example.mp4** - example of animation of such event

To run the program, you can use the command ```python3 Microlensing.py```

## Required python packages
- json
- astropy
- numpy
- matplotlib
- tqdm
- os
- subprocess

You can install them using ```pip``` or ```pip3```

## INPUT-file includes

- Source distant in kpc.
- Baseline flux.
- Source relative velocity in km/s.
- Lens mass in Solar masses.
- Lens distant in kpc.
- Impact parameter.
- Time of magnification peak.
