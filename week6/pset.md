# Week 6 PSet
Filter Specifications:

(1) Discrete-time passband edge: $$\frac{3}{28} \pi \text{ rad/s }$$ 

(2) Discrete-time stopband edge: $$\frac{5}{28} \pi \text{ rad/s }$$ 

(3) Maximum gain in the passband: 0 dB.

(4) Minimum gain in the passband: -1 dB.

(5) Maximum gain in the stopband: -50 dB.

### Butterworth
The Butterworth filter has the following magnitude response:

![**Figure 1.** Butterworth low pass filter response.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworth.png) 

![**Figure 2.** Butterworth low pass filter passband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworthpb.png) 

![**Figure 3.** Butterworth low pass filter stopband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworthsb.png) 

This filter satisfies the specification as the as the gain in the passband is above -1 dB and the stopband attentuation is above 50 dB. As the Butterworth filter transfer function is given by:

$$G(\Omega) = \frac{1}{\sqrt{1+(\frac{\Omega}{\Omega_c})^{2n}}} $$

The filter magnitude response is monotonically decreasing, thus the stopband attentuation at the maximum frequency in the stopband in Figure 3 is sufficient to show that it meets the requirements for the stopband attentuation. It is worth noting the filter is also maximally flat in the passband. 

## Chebyshev I
The Chebyshev I filter has the following magnitude response:

![**Figure 4.** Chebyshev I low pass filter response.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1.png) 

![**Figure 5.** Chebyshev I low pass filter passband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1pb.png) 

![**Figure 6.** Chebyshev I low pass filter stopband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1sb.png) 

This filter satisfies the specification as the as the gain in the passband is above -1 dB and the stopband attentuation is above 50 dB. As the Butterworth filter transfer function is given by:

$$G(\omega) = \frac{1}{\sqrt{1+\omega^{2n}}} $$

The filter magnitude response is monotonically decreasing, thus the stopband attentuation at the maximum frequency in the stopband in Figure 3 is sufficient to show that it meets the requirements for the stopband attentuation. It is worth noting the filter is also maximally flat in the passband. 