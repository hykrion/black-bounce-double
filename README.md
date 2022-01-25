<p align="center">
  <img src="/img/sigma-l-wh.PNG">
</p>

# black_bounce_double
Implementation of a black bounce (BB) where we can bounce between Schwarzschild an wormhole (WH) configuration just changing a single parameter.
The two main papers to develop this code were:

- Merce Guerrero, Gonzalo J. Olmo, Diego Rubiera-Garcia, Diego
Sáez-Chillón Gómez, Shadows and optical appearance of black bounces illuminated
by a thin accretion disk arXiv:2105.15073v2 [gr-qc] DOI 10.1088/1475-
7516/2021/08/036
-Adria Delhom, Caio F. B. Macedo, Gonzalo J. Olmo, Luís C. B.
Crispino (2019), Absorption by black hole remnants in metric-affine gravity.
arXiv:1906.06411v1 [gr-qc] DOI 10.1103/PhysRevD.100.024016

This software was developed using [GSL library](https://www.gnu.org/software/gsl/) I've done all the development in Windows so I've used the [MSYS2](https://www.msys2.org/) system. You can use [Chocolatey](https://chocolatey.org/) or install MSYS2 directly.

After that you just have to install the GSL packages using [pacman](https://archlinux.org/pacman/pacman.8.html). pacman is very convenient...

# GUI
The GUI was made with [Tcl/Tk](https://www.tcl.tk/) You need one tcl distribution to use the GUI (but you don't need the GUI to use the program). Another useful tool is [gnuplot](http://www.gnuplot.info/), you can install it with pacman too.
This's how the GUI looks like.

![Black Bounce GUI](/img/gui-1.PNG)

![Black Bounce GUI](/img/gui-2.PNG)

These're the results for a BH, using *aa = 1.5*

![Black Bounce - BH tortoise](/img/tortoise-bh.PNG)

![Black Bounce - BH wave potential](/img/wave-pot-bh.PNG)

![Black Bounce - BH reflexion transmision coefficients](/img/reflexion-transmision-bh.PNG)

![Black Bounce - BH sigma-l](/img/sigma-l-bh.png)

These're the results for a WH, using *aa = 2.5*

![Black Bounce - BH tortoise](/img/tortoise-wh.PNG)

![Black Bounce - BH wave potential](/img/wave-pot-wh.PNG)

![Black Bounce - BH reflexion transmision coefficients](/img/reflexion-transmision-wh.PNG)

![Black Bounce - BH sigma-l](/img/sigma-l-wh.png)
