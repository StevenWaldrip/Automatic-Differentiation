# Automatic-Differentiation

An automatic differentiation implimentation for [Octave](https://octave.org/).

It can handle the basic operators, most core functions, matrix factorizations, and functions. See [my site](https://mathsfromnothing.au/scalar-functions-vector-matrix-and-tensor-functions) for a complete list and a detailed discription of how it works. The code is in a single file ad.m that can be downloaded from the ‘implementation and how to use it’ section of the page linked above or ad.m from the inst directory. To test it copy adtest.m into the same directory, it should test all of the code paths. It currently does not support complex numbers. They can be used but the results may be incorect.

To use it copy ad.m into a directory. Open that directory in octave. To find the derivative of a scalar function at several locations the following code can be used
```
x=ad(linspace(-2,2));
y=x.^2;
plot(x.value,y.value,x.value,diag(y.diff(x)));
```
For the derivative with respect to a matrix or vector, for example, a multivariable normal distribution with respect to the mean vector and covarence matrix the following code can be used
```
mu=ad([2;3]);
Sigma=ad([4 1;1 5]);
x=[1;4];%not taking the derivative with respect to x, so not ad variable
P=@(x,mu,Sigma)exp(-.5*(x-mu)'*(Sigma\(x-mu)))/(2*pi^(2/2)*det(Sigma)^.5);
p=P(x,mu,Sigma)
dpdmu=full(p.diff(mu))
dpdmucheck=ad.fd(@(mu)P(x,mu,Sigma.value),mu.value)
dpdSigma=full(p.diff(Sigma))
dpdSigmacheck=ad.fd(@(Sigma)P(x,mu.value,Sigma),Sigma.value)
```
