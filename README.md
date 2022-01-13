# Transmission-Loss-reactive-muffler
In an acoustics course for the sound engineering degree program (UNTREF), I developed this Python script to analytically calculate the transmission loss of a passive reactive silencer.

A silencer is an important noise control element for reduction of machinery exhaust noise, fan noise, and other noise sources involving flow of a gas. In general, a silencer may be defined as an element in the flow duct that acts to reduce the sound transmitted along the duct while allowing free flow of the gas through the flow passage. A silencer may be passive, in which the sound is attenuated by reflection and absorption of the acoustic energy within the element.


Reactive mufflers consist typically of several pipe segments that interconnect a number of a larger diameter chambers. These devices reduce the radiated acoustical power primarily throught impedance mismatch, that is, through the use of acoustical impedance discontinuites to reflect sound back toward the source. The TL is the acoustical power-level difference between the incident and transmitted waves of an anechoically terminated silencer.

With the hypothesis of plane wave propagation, the different variables present between the input and output of a system, through the parameters of a quadrupole. This analytical method is used for the analysis of one-dimensional systems connected in series (such as silencers), making faster in this case to the calculation of the transmission loss. The equivalent transfer matrix of the system (Lumped matrix) gives us allows us to become independent from geometry, and is the result of the product of the individual matrices corresponding to each element.

In this script you can calculate the TL by entering the dimensions of each stage of the silencer and characteristics of the fluid (the values that are assigned in the input variables are by default, modify them). You can also find functions that allow quick calculations of dimensions associated with a chosen frequency, and the calculation of the cutoff frequency to be able to verify if you are working under the hypothesis of plane waves.

Libraries needed:
* `numpy`
* `matplotlib`
