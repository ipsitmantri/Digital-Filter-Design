# Digital Filter Design
My solution for the Filter Design Assignment of the course EE 338: Digital Signal Processing
As a part of this assignment, I had to design 6 different filters based a set of specifications given. The detailed information is available in the [report](./180070032.pdf)
These are the following filters that have been designed:

---
## Infinite Impulse Response Filters
1. Bandpass filter using Butterworth design
2. Bandstop filter using Chebyshev design

## Finite Impulse Response Filters
1. Bandpass filter using Kaiser window
2. Bandstop filter using Kaiser window

## Elliptic Filters
1. Bandpass filter
2. Bandstop filter

---
## Some observations
After realising the bandpass and bandstop filters for a given set of specifications
in various forms, I was able to make the following conclusions –
- Realisation using an elliptic approximation gives the least order, i.e, it
requires mininum resources. The order of resources required is as follows
(for a given set of specifications):  
```Elliptic (3) < Chebyshev (4) < Butterworth (8) < Kaiser FIR (50-60)```
- The elliptic filter has the sharpest transition from passband to stopband or
vice versa for the same specifications. In particular, the decreasing order of
transition sharpness is:  
```Elliptic > Chebyshev > Butterworth > Kaiser FIR```
- The FIR realisation gives a linear phase response in the passband region.
The nonlinearity of the phase response decreases in the following order:  
```Elliptic > Chebyshev > Butterworth > Kaiser FIR```

