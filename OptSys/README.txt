This README is to keep track of ideas, organizational structure, philosophy, etc.

==================================================
                    GOALS
==================================================

1) Store an Optical System of Arbitrary Complexity
2) Use Fresnel Propagation for Diffraction Analysis
3) Design Optical Elements
4) Polarization?


=================================================
                INSPIRATION
=================================================

1) Draws from the ideas/methods of AOSim2 by JLCodona
2) Draws from the ideas/methods of PIAACMCdesign/Cfits by oguyon


=================================================
                 PHILOSOPHY
=================================================

1) Object Oriented
2) Fresnel Propagation
3) Do things on GPU when possible
    a) Efficiency is a must
4) Write Results to .txt and .fits files
    a) Only display if asked
    b) Create a filesystem
5) Be capable of controlling hardware?
6) GUI?

=================================================
            CLASS STRUCTURE LIST
=================================================

1) Coordinate Grid - Should set up basics of everything
    a) Storage of arrays
    b) pixel-to-real space
    c) simple mathematical methods
    d) interpolation

2) Optical System
    a) Source Size
    b) Wavelength/s
    c) Number of Elements
    d) F-number
    e) Element Structure
        A) Include classes for individual elements?
        B) Only some specific elements?
             i)  DM
            ii)  WFS
           iii)  PIAA
            iv)  FPM
             v)  Aperture
            vi)  LCC

3) Light Field
    a) Propagation Methods
    b) PSF computation
    c) amplitude and phase

4) Hardware Controller

5) Utility Functions
================================================