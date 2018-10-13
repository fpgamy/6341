# Week 6 PSet
Filter Specifications:

(1) Discrete-time passband edge: $\frac{3}{28} \pi \text{ rad/s }$ 

(2) Discrete-time stopband edge: $\frac{5}{28} \pi \text{ rad/s }$ 

(3) Maximum gain in the passband: 0 dB.

(4) Minimum gain in the passband: -1 dB.

(5) Maximum gain in the stopband: -50 dB.

### Butterworth
The Butterworth filter has the following magnitude response:

![**Figure 1.** Butterworth low pass filter response.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworth.png) 

![**Figure 2.** Butterworth low pass filter passband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworthpb.png) 

![**Figure 3.** Butterworth low pass filter stopband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/butterworthsb.png) 

This filter satisfies the specification as the as the gain in the passband is above -1 dB and the stopband attentuation is above 50 dB. As the Butterworth filter transfer function is given by:

$$|H_c(j\Omega)|^2 = \frac{1}{\sqrt{1+(\frac{j\Omega}{j\Omega_c})^{2n}}} $$

The filter magnitude response is monotonically decreasing, thus the stopband attentuation at the maximum frequency in the stopband in Figure 3 is sufficient to show that it meets the requirements for the stopband attentuation. It is worth noting the filter is also maximally flat in the passband. 

## Chebyshev I
The Chebyshev I filter has the following magnitude response:

![**Figure 4.** Chebyshev I low pass filter response.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1.png) 

![**Figure 5.** Chebyshev I low pass filter passband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1pb.png) 

![**Figure 6.** Chebyshev I low pass filter stopband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev1sb.png) 

This filter satisfies the specification as the as the gain in the passband equiripple and the maximum deviation from 0 dB is -1dB and the stopband attentuation is above 50 dB. As the Chebyshev I filter transfer function is given by:

$$|H_c(j\Omega)|^2 = \frac{1}{1+(V_N^2(\frac{j\Omega}{j\Omega_c})\epsilon^2 } $$
Where 
$$V_N(x) = \cos(N\cos^{-1}(x))$$
inside the passband, the value of $x$ varies from 0 to 1, and the value inside the cosine varies between 0 and 1. As this is scaled by $\epsilon^2$, the magnitude response squared varies between 1 and $\frac{1}{1+\epsilon^2}$, thus giving the equiripple behaviour. Increasing $\Omega$ beyond $\Omega_c$, the inverse cosine term in $V_N$ is imaginary and positive (refer to Figure 7). For $x \geq 1$ the imaginary part is also monotonically increasing. The cosine expression can be considered as applying an hyperbolic cosine function to the imaginary part of the argument. As the hyperbolic cosine is monotonically increasing for positive values of $x$, the stopband is monotonically decreasing. Thus, as with Butterworth filters, the stopband attentuation at the maximum frequency in the stopband being less than or equal to the maximum stopband gain is sufficient to show that it meets the requirements for the stopband attentuation.

![**Figure 7.** Imaginary part of the function $\text{cos}^{-1}(x)$.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/imacos.png) 


## Chebyshev II
The Chebyshev II filter has the following magnitude response:

![**Figure 8.** Chebyshev II low pass filter response.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev2.png) 

![**Figure 9.** Chebyshev II low pass filter passband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev2pb.png) 

![**Figure 10.** Chebyshev II low pass filter stopband.](https://raw.githubusercontent.com/fpgamy/6341/master/week6/chebyshev2sb.png) 

This filter satisfies the specification as the as the gain in the passband gain is above -1 dB and the maximum gain in the stopband is above -50 dB. As the Chebyshev II filter transfer function is given by:

$$|H_c(j\Omega)|^2 = \frac{1}{1+[(V_N^2(\frac{j\Omega_c}{j\Omega})\epsilon^2]^{-1}} $$

Inside the stopband, the value of $x = \frac{j\Omega_c}{j\Omega}$ varies from 0 to 1, and the value inside the cosine varies between 0 and 1. As this is scaled by $\epsilon^2$, the magnitude response squared varies between 1 and $\frac{1}{1+\epsilon^2}$, thus giving the equiripple behaviour in the stopband. 


