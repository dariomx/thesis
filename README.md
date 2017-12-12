This is the thesis of my Msc. in Computer Science. While the program specialized
in parallel and distributed computing, I chose as topic Numerical Linear Algebra.
This was because I was heading towards Applied Mathematics (Data Science).

In concrete, the thesis revolves around proposing a better (faster) algorithm
for computing the Fiedler Vector (second eigenvector of Laplacian matrix). The
context was an application doing spectral clustering, and among the requirements,
they wanted the computation to occur serially.

The proposal was LOBPCG, which can compute the required eigenvector in subsecond
scale; as opposed to the minute scale of the dense-matrix algorithm that the
application was using. I really had a lot of fun with this project.

