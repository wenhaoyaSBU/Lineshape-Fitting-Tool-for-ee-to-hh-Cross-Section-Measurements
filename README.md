# Cross Section Lineshape Fitting Package


## How to use
* Put your signal event number measurement results, as well as other quantities like the efficiencies, luminosities, and uncertainties, in `data/para.txt`
* Change the vaccum polarization file `vaccFile/vacc_nojpsi.dat` if you are not working near $J/\psi$ resonance energy
* Do `make`, and then do `./fit_XS`.


## Outputs
* `output/calculated_xs_data.txt`: Calculated cross section from your input `para.txt` at each energy point.
* `output/getpoint*.txt`: Fitted cross section (in $\mathrm{pb}^{-1}$) at C.O.M energy grid points. Can be used to draw the fitting result.
* `output/fitted*.txt`: The post-fit parameters and their uncertainties, defined in the construction of `gMinuit` in `fit_XS.cxx`.


## Fitting Model
$$
\sigma^{\text{Born}}(W; M, \Gamma, \Gamma_{ee}^{(0)}, \mathcal{F}_{\omega\pi^0}, C, \phi) = \frac{\pi \alpha^2}{6}\left(\frac{\mathcal{F}_{\omega\pi^0}}{W^4}\right)^2 q_f^3 \bigg |1 + \frac{3W^2\Gamma_{ee}^{(0)}}{M\alpha(W^2-M^2+iM\Gamma)}Ce^{i\phi}\bigg |^2 ,
$$
where $q_f$ is a dimensionless quantity relevant to the momentum of the final state particles (which are in our case $\omega$ and $\pi^0$) defined as
$$
q_f = \sqrt{\left(1+\frac{m_{\omega}^2 - m_{\pi^0}^2}{W^2}\right)^2 - \frac{4m_{\omega}^2}{W^2}}.
$$
After considering vaccum polarization effects, the continuum amplitude gets modified and the width $\Gamma_{ee}$ becomes the "physical" value
$$
\sigma^{\text{dressed}}(W; M, \Gamma, \Gamma_{ee}, \mathcal{F}_{\omega\pi^0}, C, \phi) = \frac{\pi \alpha^2}{6}\left(\frac{\mathcal{F}_{\omega\pi^0}}{W^4}\right)^2 q_f^3 \bigg |\frac{1}{|1-\Pi_0(W)|} + \frac{3W^2\Gamma_{ee}}{M\alpha(W^2-M^2+iM\Gamma)}Ce^{i\phi}\bigg |^2 ,
$$
and then the ISR correction
$$
\sigma^{\text{ISR}}(W; M, \Gamma, \Gamma_{ee}, \mathcal{F}_{\omega\pi^0}, C, \phi) = \int_0^{1-(\frac{W_{\text{min}}}{W})^2}\text{d}x F(x, W)\sigma^{\text{dressed}}(W\sqrt{1-x}),
$$
finally the beam energy spread effects
$$
\sigma^{\text{obs}}(W; M, \Gamma, \Gamma_{ee}, \mathcal{F}_{\omega\pi^0}, C, \phi, S_E) = \int_{W - n S_E}^{W + n S_E}\text{d}W' \sigma^{\text{ISR}}(W') \cdot \frac{1}{S_E\sqrt{2\pi}}\exp
    \bigg (-\frac{(W' - W)^2}{2S_E^2}\bigg ).
$$

* <span style="color:red;">!!!NOTE!!!</span>: We use an analytic formula to calculate the ISR correction, but the effect from $q_f$ is considered by Taylor-expanding it into power series of $(W-W_\mathrm{min})$ and then truncate at 5th order, then calculate the $\sigma^{\text{ISR}}$ as the sum of 0th to 5th order corrections. <span style="color:red;">So please modify the `q0Corr_2body` - `q5Corr_2body` variables in `headers/variables.h` if you use a different final state!!!</span>
* If you want to turn off the correction, simply set `q0Corr_2body` to the $q_f$ value calculated using quantities from your final state and then set others to 0.
* Also, by default the $W^{-8}$ dependence of Born cross section is implemented, as in the first equation. Within the ISR integral bases we currently have, the $W^{-7}$ and $W^{-6}$ dependencies are also available. You can change to models with these dependencies by modifying the `Ana` function in `headers/physicsFuncs.h`, by uncommenting the corresponding lines.
* Besides, by default we use an energy cut $W_\mathrm{min} = 2.8$ GeV. This should also be modified to the threshold you use.
* You can change the cost function by modifying the `fcn` function in `headers/variables.h`.

## Cost Function
The cost function is the function to be minimized in the fit procedure. In our script we use a least-square fit, and the cost function is

$$
\begin{aligned}
&\chi^2_\text{mixed}(M, \Gamma, \Gamma_{ee}, \mathcal{F}_{\omega\pi^0}, C, \phi, S_E, f) \notag \\
    &= 
    \sum_{i = 1}^{24}
    \bigg (
    \frac{\sigma^\text{obs}(W_i; M, \Gamma, \Gamma_{ee}, \mathcal{F}_{\omega\pi^0}, C, \phi, S_E) - f \cdot \sigma^\text{mea}_i}{\Delta \sigma^\text{mea}_i}
    \bigg )^2 \notag \\
    &+ \sum_{i = 1}^{24}
    \bigg (
    \frac{W_i - W_i^\text{Prop.}}{\Delta W_i}
    \bigg )^2 
    + \bigg ( 
    \frac{S_E - \overline{S_E}}{\Delta S_E} 
    \bigg )^2 
    + \bigg ( 
    \frac{M - \overline{M}}{\Delta M} 
    \bigg )^2 \notag \\
    &+ 
    \bigg (
    \frac{\Gamma - \overline{\Gamma}}{\Delta \Gamma}
    \bigg )^2 
    + 
    \bigg (
    \frac{\Gamma_{ee} - \overline{\Gamma_{ee}}}{\Delta \Gamma_{ee}}
    \bigg )^2 
    + 
    \bigg (
    \frac{f - 1}{\Delta f}
    \bigg )^2, 
\end{aligned}
$$

where the subscript `mixed` means the amplitude `C` is set to float. If a process is asuumed to have only EM contribution to its resonance decay amplitude, `C = 1`. Otherwise `C` has to be extracted from the fit. 
* There are nuisance parameters (NP) for the external measurements of quantities $W_i, S_E, M, \Gamma, \Gamma_{ee}, f$. $f$ is the NP reponsible for correlated systematic uncertainty, which can be modified in `para.txt`. The central values and uncertainties of other quantities can be modified in `headers/variables.h`.


## Formula Notes
* In the directory `notes` there are three PDF files, explaining the formula used in the codes. 
* The main note is `ISR_Integral_W_neg8.pdf` which explains how we expand the 2-body phase space factor in VP decays $q_f^3$ into power series in our calculation, and how to decompose the ISR integral into different integral bases. (We do this expansion because it is ensured that the ISR integral with Born cross section having integer $W$ dependence in the prefactor will have a closed-form analytic expression in terms of special functions.)
* The note `ISR_Basis_Extension.pdf` gives analytic formula for the bases integrals.
* The note `PropagatorIntegral.pdf` gives the derivation of propagator-related ISR basis integral, in which the Appell's hypergeometric funcition $F_1$ is inevitable. The evaluation of $F_1$ involves a double series, which takes up most of the computation time. In this note we try to find a way that is the most efficient to evaluate this integral.

## Acknowledgment
Wenhao Ya would like to thank Zhikun Xi, Yiqi Du, Baoxin Liu, Yanning Wang, Guangbao Sun, Xiang Zhou, and Zhenyu Zhang, for their pioneering work in finding the analytic form of ISR integrals and the fitting code construction, as well as the valuable discussions regarding this project.