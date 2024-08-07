# Hartle-Thorne perturbative magnetar model (HTmag) used in my Master's thesis

**WARNING:** This code is no longer developed/maintained!

Main code written for my Master's thesis
>Martin J. Steil, *Structure of slowly rotating magnetized neutron stars in a perturbative approach*,<br>
>Master's thesis, Technische Universität Darmstadt, 2017,<br>
>[`mjsteil.github.io/files/MJSteil_MSc_Thesis_2017.pdf`](https://mjsteil.github.io/files/MJSteil_MSc_Thesis_2017.pdf)

for various numerical computations regarding isolated neutron stars as perturbations to the Tolman-Oppenheimer-Volkoff (TOV) solution.

Includes a Mathematica implementation in [`I-love-Qrot.nb`](https://github.com/MJSteil/HTmag/blob/main/I-love-Q/I-love-Qrot.nb) for the computation of I-love-Q relations (see, e.g. [arxiv.org/abs/1303.1528](https://arxiv.org/abs/1303.1528) for details).

### Dependencies
- [`gsl_wrapper`](https://github.com/MJSteil/gsl_wrapper) – Used for interpolation and numerical (ODE) integration
- [`matplotlib-cpp`](https://github.com/MJSteil/matplotlib-cpp) – Used for plotting (only necessary if one wants to plot them using `matplotlib_cpp`
- [`eos`](https://github.com/MJSteil/eos) – Used for the equation of state
- [`units.hpp`](https://github.com/MJSteil/units) – Used for physical units

- [`LiSK`](https://github.com/DSNJTurtle/LiSK) – Lightweight C++ library for the numerical evaluation of classical polylogarithms Li(n,x) 

<i>Martin Jakob Steil, 2017</i>

