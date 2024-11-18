# Barnes' double Gamma function

This package provides an implementation of the multiple $\Gamma$ functions $G, \Gamma_2, \Gamma_b$ as defined in
[this wikipedia article](https://en.wikipedia.org/wiki/Multiple_gamma_function).

## Algorithm

Internally, the double gamma functions are computed using a product formula from [this paper](https://arxiv.org/abs/2208.13876v1):

$$G(z ; \tau)=\frac{1}{\tau \Gamma(z)} e^{\tilde{a}(\tau) \frac{z}{\tau}+\tilde{b}(\tau) \frac{z^2}{2 \tau^2}} \prod_{m \geq 1} \frac{\Gamma(m \tau)}{\Gamma(z+m \tau)} e^{z \psi(m \tau)+\frac{z^2}{2} \psi^{\prime}(m \tau)},$$

where

$$\begin{aligned}
C(\tau) & :=\sum_{k=1}^{m-1} \psi(k \tau)+\frac{1}{2} \psi(m \tau)-\frac{1}{\tau} \ln \left(\frac{\Gamma(m \tau)}{\sqrt{2 \pi}}\right)-\frac{\tau}{12} \psi^{\prime}(m \tau) \\
& +\frac{\tau^3}{720} \psi^{(3)}(m \tau)-\frac{\tau^5}{30240} \psi^{(5)}(m \tau)+\frac{\tau^7}{1209600} \psi^{(7)}(m \tau)+O\left(m^{-9}\right) \\
D(\tau) & :=\sum_{k=1}^{m-1} \psi^{\prime}(k \tau)+\frac{1}{2} \psi^{\prime}(m \tau)-\frac{1}{\tau} \psi(m \tau)-\frac{\tau}{12} \psi^{\prime \prime}(m \tau) \\
& +\frac{\tau^3}{720} \psi^{(4)}(m \tau)-\frac{\tau^5}{30240} \psi^{(6)}(m \tau)+\frac{\tau^7}{1209600} \psi^{(8)}(m \tau)+O\left(m^{-10}\right) .
\end{aligned}$$

and

$$\tilde{a}(\tau)=a(\tau)-\gamma \tau=\frac{\tau}{2} \ln (2 \pi \tau)+\frac{1}{2} \ln (\tau)-\tau C(\tau), \quad \tilde{b}(\tau)=b(\tau)+\frac{\pi^2 \tau^2}{6}=-\tau \ln (\tau)-\tau^2 D(\tau),$$

the functions $\Gamma$ and $\psi$ are the Euler gamma and digamma function:

$$\psi(z) = (\log \Gamma)'(z).$$

## Functions exported by the module

```@docs
logdoublegamma
```

```@docs
doublegamma
```

```@docs
log_barnesdoublegamma
```

```@docs
barnesdoublegamma
```

```@docs
loggamma2
```

```@docs
gamma2
```
