# ------------------ description of the algorithm ------------------------ #
# The original objective function is ||y - X * beta||_1 / n + ||y -
# Z * gamma||_1 / n + lam1 * ||beta||_1 + lam2 * ||gamma||_1 + lam3 * 
# ||Eta * beta - gamma||_1. 
# We would like to use linear programming to solve this problem.
# The new objective function is alpha+(1, n) @ alpha-(1, n) @ pi+(1, n) @ pi-(1, n) @
# lam3 * tau+(1, q) @ lam3 * tau-(1, q) @ lam1 * beta+(1, p) @ lam1 * beta-(1, p) @
# lam2 * gamma+(1, q) @ lam2 * gamma-(1, q). We want to minimize this function.
# +- means positive and negative part. @ means sum. (1, q) means it's a q-dimensional
# row vector.
# The conditions are alpha+-, beta+-, pi+-, tau+-, gamma+- >=0.
# alpha+ - alpha- = (y - X * (beta+ - beta-)) / n
# pi+ - pi- = (y - Z * (gamma+ - gamma-)) / n
# tau+ - tau- = Eta * (beta+ - beta-) - gamma+ + gamma-
# ------------------------------------------------------------------ #

# define some functions to simplify the vector/matrix generation process.
library(lpSolve)
ed = function(p, val){
  rep(val, p)
} # vector of length p with all elements equal to val.
dg = function(p, val){
  diag(rep(val, p))
} # p by p diagonal matrix with diagonal elements val.
mx = function(m, n, val){
  matrix(rep(val, m * n), nrow = m, ncol = n)
} # m by n matrix with all elements equal to val. 

# standardization
standardize = function(x)
{
   p = ncol(x)
   n = nrow(x)
   center = colMeans(x)
   x.mean = x - matrix(rep(center, n), n, p, byrow = T)
   scale = sqrt(colSums(x.mean ^ 2) / n)
   xx = t(t(x.mean) / scale)
   list(xx = xx, center = center, scale = scale)
}

# main function to find the ARMI regression coefficients
ARMI = function(y, X, Z, Eta, lam1, lam2, lam3){
  n = length(y)
  p = ncol(X)
  q = ncol(Z)
  if (n != nrow(X) | n != nrow(Z) | p != ncol(Eta) | q != nrow(Eta)) 
    stop("Dimensions do not match!")
  
  # normalization
  Y = y - mean(y)
  dat.x = standardize(X)
  dat.z = standardize(Z)
  x = dat.x$xx
  z = dat.z$xx
  scale.x = dat.x$scale
  scale.z = dat.z$scale
    
  # Specify the coefficients for the objective function, it should be of length
  # 4 * n + 2 * q + 2 * p + 2 * q.
  c.obj = c(ed(n, 1), ed(n, 1), ed(n, 1), ed(n, 1), 
            ed(q, lam3), ed(q, lam3), 
            ed(p, lam1), ed(p, lam1), 
            ed(q, lam2), ed(q, lam2)) 
  
  #coef for alpha+ - alpha- = (y - x * beta) / n
  cmp1 = cbind(dg(n, 1), dg(n, -1), mx(n, 2 * n + 2 * q, 0), x / n, -x / n, 
                mx(n, 2 * q, 0))
  
  #coef for pi+ - pi- = (y - z * gamma) / n
  cmp2 = cbind(mx(n, 2 * n, 0), dg(n, 1), dg(n, -1), mx(n, 2 * q + 2 * p, 0), 
                z / n, -z / n)
  
  #coef for tau+ - tau- = Eta * beta+ - Eta * beta- - gamma+ + gamma-
  cmp3 = cbind(mx(q, 4 * n, 0), dg(q, 1), dg(q, -1), -Eta %*% diag(1 / scale.x), 
               Eta %*% diag(1 / scale.x), diag(1 / scale.z), -diag(1 / scale.z))
  
  c.con = rbind(cmp1, cmp2, cmp3)
  c.dir = rep("=", 2 * n + q)
  c.rhs = c(Y / n, Y / n, rep(0, q))
  lps = lp(direction = "min", objective.in = c.obj, const.mat = c.con, 
           const.dir = c.dir, const.rhs = c.rhs)
  if (lps$status != 0) cat("no feasible solution\n")
  bp = lps$solution[(4 * n + 2 * q + 1):(4 * n + 2 * q + p)]
  bm = lps$solution[(4 * n + 2 * q + p + 1):(4 * n + 2 * q + 2 * p)]
  gp = lps$solution[(4 * n + 2 * q + 2 * p + 1):(4 * n + 3 * q + 2 * p)]
  gm = lps$solution[(4 * n + 3 * q + 2 * p + 1):(4 * n + 4 * q + 2 * p)]
  bet = bp - bm
  gamm = gp - gm
  betahat = bet / scale.x
  gammahat = gamm / scale.z
  return(list(betahat = betahat, gammahat = gammahat,
              bet = bet, gamm = gamm))
}
# -------------- end of the algorithm --------------------- #


