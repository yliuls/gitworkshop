library(io)


## test 

## test 


set.seed(1111);

N <- 500;
J <- 100;
K <- 2;  # tumour, normal


# test 2


# test 2 



s0 <- 0.5;
d0 <- 10;
kappa <- 10;
gamma <- 0.001;


# test 3

# test 3


w <- rbeta(N, 2, 2);
w <- cbind(w, 1-w);

# NB lambda and phi would be difficult to estimate a priori
lambda <- rgamma(J, 5*2, 2);
phi <- c(0.8, 0.4);

# J x K
v <- matrix(unlist(lapply(phi, function(p) rbinom(J, 1, p))), nrow=J, ncol=K);

# N x J x K
z <- array(0, dim = c(N, J, K));

# test 3

for (i in 1:N) {
	z[i, , ] <- ifelse(v == 0, 
		rexp(J*K, kappa),
		rgamma(J*K, lambda/gamma, 1/gamma)
	);
}


## test 4


# N x J
mu <- log(w[ ,1] * exp(z[, , 1]) + w[, 2] * exp(z[, , 2]));

tau <- rgamma(J, d0/2, d0*s0*s0/2);
summary(tau^-0.5)



## test 4
# N x J
# rnorm takes sd instead of var, so we need to sqrt
sigma <- tau^-0.5;
summary(sigma)
y <- matrix(rnorm(N*J, mu, rep(sigma, each=N)), nrow=N, ncol=J);


# visualize

hist(w[, 1]);
hist(lambda)
hist(z, breaks=100)
hist(mu, breaks=100)
hist(y, breaks=100)

j <- 2;
v[j, ]

plot(z[, j, 1], z[, j, 2])
plot(w[, 1], mu[, j])
plot(w[, 1], exp(mu[, j]))
plot(w[, 1], exp(y[, j]))

# test 5

# test 5

# write to file

# test 5
sim <- list(
	data = list(
		N = N, J = J, K = K,
		d0 = d0, s0 = s0,
		kappa = kappa, gamma = gamma,
		lambda = lambda, phi = phi,
		w = w,
		y = y
	),
	params = list(
		v = v, z = z, mu = mu,
		tau = tau
	)
);

qwrite(sim, "sim-data.rds");
qwrite(sim, "sim-data.json");
