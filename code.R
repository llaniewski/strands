library(rgl)
library(Rcpp)

N = 20
k = 4
n = k+(k-1)*4
ord = 2
R = 0.02
L = 0.7


p = matrix(runif(N*3,L/2,1-L/2),N)
v = matrix(rnorm(N*3),N); v[,3] = v[,3]/10
v = v/sqrt(rowSums(v^2))

ret = expand.grid(t=seq(-L/2,L/2,len=k),i=seq_len(nrow(p)))
p = ret$t*v[ret$i,]+p[ret$i,]

Rcpp::sourceCpp("code.cpp")

bs = bsplines(k,n,ord)
plot3d(p, type="n", asp="iso",xlim=0:1,ylim=0:1,zlim=0:1)

v = p;
v[] = 0;
dt = 0.0002;
for (it in 1:100000) {
  p = p + v*dt;
  m = fun(p,k,n,ord,R,0.1,0.1,L);
  m[,3] = m[,3] - 0.01;
  v = v + m*dt;
  v = v-0.2*dt*v;
  print(range(v))
  
  if (it %% 2000 == 0) {
    ret = expand.grid(t=c(0,5),i=seq_len(nrow(p)))
    segments3d(p[ret$i,] + ret$t*v[ret$i,])
    segments3d(p[ret$i,] + ret$t*m[ret$i,],col=3)
    clear3d();
    points3d(expand.grid(c(-0.1,1.1),c(-0.1,1.1),c(-0.1,1.1)))
    quads3d(c(0,1,1,0),c(0,0,1,1),c(0,0,0,0))
    np = matrix(t(bs) %*% matrix(p,k),ncol=3)
    spheres3d(np,radius = R,col=rep(seq_len(nrow(p)/k)+1,each=n))
  }
}


          