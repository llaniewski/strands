#include <Rcpp.h>
#include "spline.h"

using namespace Rcpp;

typedef std::array<size_t,2> edge_t;
typedef std::array<double,3> point_t;

template <typename T, size_t N>
T dot(const std::array<T,N>& x, const std::array<T,N>& y) {
  T sum = 0;
  for (size_t i=0; i<N; i++) sum += x[i]*y[i];
  return sum;
}

double clamp(double x) {
  if (x<0) return 0;
  if (x>1) return 1;
  return x;
}

std::pair<double, double> touch(point_t p0, point_t p1, point_t p2, point_t p3) {
  double t = 0, s = 0;
  point_t a0,a1,a2;
  for (size_t i=0; i<3; i++) {
    a0[i] = p2[i] - p0[i];
    a1[i] = p1[i] - p0[i];
    a2[i] = p3[i] - p2[i];
  }
  double M11 =  dot(a1,a1);
  double M12 = -dot(a1,a2);
  double M21 =  dot(a2,a1);
  double M22 = -dot(a2,a2);
  double RHS1 = dot(a0,a1);
  double RHS2 = dot(a0,a2);
  double det = M11*M22 - M12*M21;
  t = (  M22*RHS1 - M12*RHS2)/det;
  s = (- M21*RHS1 + M11*RHS2)/det;
  return std::make_pair(clamp(t),clamp(s));
}



// [[Rcpp::export]]
NumericMatrix bsplines(int k, int n, int ord) {
  NumericMatrix bs(k,n);
  for (size_t i=0; i<k; i++)
    for (size_t j=0; j<n; j++)
      bs(i,j) = bspline_b(j/(n-1.0), k, i, ord, false);
  return bs;
}
  


// [[Rcpp::export]]
NumericMatrix fun(NumericMatrix points, int k, int n, int ord, double R, double K, double C, double L) {
  const size_t DIM = 3;
  NumericMatrix ret(points.nrow(),points.ncol());
  assert(points.ncol() == DIM);
  std::vector<double> tab1(k*n);
  auto bs = [&](size_t i, size_t j) -> double& { return tab1[i+k*j]; };
//  inline double bspline_b(double x, int n, int w, int k, bool per)
  for (size_t i=0; i<k; i++)
    for (size_t j=0; j<n; j++)
      bs(i,j) = bspline_b(j/(n-1.0), k, i, ord, false);
  
  size_t N = points.nrow()/k;
  assert(N*k == points.nrow());
  
  

  auto interp = [&](size_t s, size_t i) -> point_t {
    point_t x{};
    for (size_t g=0; g<k; g++) {
      double v = bs(g,i);
      for (size_t d=0; d<3; d++) {
        x[d] += points(s*k+g, d) * v;
      }
    }
    return x;
  };
  
  auto force = [&](point_t p1, point_t p2) -> point_t {
    point_t d{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
    double dist_2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
    double dist = sqrt(dist_2);
    if (dist < 2*R) {
      double w = 2*R - dist;
      double sc = -w/dist;
      return {d[0]*sc,d[1]*sc,d[2]*sc}; 
    } else {
      return {0,0,0};
    }
  };
  
  auto addforce = [&](size_t s, size_t i, point_t f) {
    for (size_t g=0; g<k; g++) {
      double v = bs(g,i);
      for (size_t d=0; d<3; d++) {
        ret(s*k+g, d) += f[d] * v;
      }
    }
  };
  

  {
    std::vector<double> tab2(n*n);
    auto dist2 = [&](size_t i, size_t j) -> double& { return tab2[i+n*j]; };
  
    for (size_t s1=0; s1<N; s1++) {
      for (size_t s2 = s1+1; s2<N; s2++) {
        for (size_t i=0; i<n; i++) {
          point_t p1 = interp(s1,i);
          for (size_t j=0; j<n; j++) {
            point_t p2 = interp(s2,j);
            point_t d{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
            dist2(i,j) = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
          }
        }
        for (size_t i=0; i<n; i++) {
          for (size_t j=0; j<n; j++) {
            double d = 1.0e9;
            if (i > 0) if (d > dist2(i-1,j)) d = dist2(i-1,j);
            if (i + 1  < n) if (d > dist2(i+1,j)) d = dist2(i+1,j);
            if (j > 0) if (d > dist2(i,j-1)) d = dist2(i,j-1);
            if (j + 1  < n) if (d > dist2(i,j+1)) d = dist2(i,j+1);
            if (d > dist2(i,j)) {
              point_t p1 = interp(s1,i);
              point_t p2 = interp(s2,j);
              point_t f = force(p1,p2);
              addforce(s1,i,f);
              for (size_t d=0; d<3; d++) f[d] = -f[d];
              addforce(s2,j,f);
            }
          }
        }
      }
    }
  }
  {
    
    const size_t W = 6;
    std::array< std::pair<size_t, double>, W > wall = {{{0,-R},{0,1+R},{1,-R},{1,1+R},{2,-R},{2,1+R}}};
    
    std::vector<double> tab2(n*W);
    auto dist2 = [&](size_t i, size_t j) -> double& { return tab2[i+n*j]; };
    auto wall_point = [&](point_t p, size_t j) -> point_t { p[wall[j].first] = wall[j].second; return p; };
    
    for (size_t s1=0; s1<N; s1++) {
      for (size_t i=0; i<n; i++) {
        point_t p1 = interp(s1,i);
        for (size_t j=0; j<W; j++) {
          point_t p2 = wall_point(p1,j);
          point_t d{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
          dist2(i,j) = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
        }
      }

      for (size_t i=0; i<n; i++) {
        for (size_t j=0; j<W; j++) {
          double d = 1.0e9;
          if (i > 0) if (d > dist2(i-1,j)) d = dist2(i-1,j);
          if (i + 1  < n) if (d > dist2(i+1,j)) d = dist2(i+1,j);
          if (d > dist2(i,j)) {
            point_t p1 = interp(s1,i);
            point_t p2 = wall_point(p1,j);
            point_t f = force(p1,p2);
            addforce(s1,i,f);
          }
        }
      }
    }
  }
  
  
  double l = L/n;
  auto force2 = [&](point_t p1, point_t p2) -> point_t {
    point_t d{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
    double dist_2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
    double dist = sqrt(dist_2);
    double w = l - dist;
    double sc = -w/dist;
    return {d[0]*sc,d[1]*sc,d[2]*sc}; 
  };
  
  
  {
    for (size_t s1=0; s1<N; s1++) {
      for (size_t i=0; i<n-1; i++) {
        point_t p1 = interp(s1,i);
        point_t p2 = interp(s1,i+1);
        point_t f = force2(p1,p2);
        addforce(s1,i,f);
        for (size_t d=0; d<3; d++) f[d] = -f[d];
        addforce(s1,i+1,f);
      }
    }
  }
  {
    for (size_t s1=0; s1<N; s1++) {
      for (size_t i=0; i<n-2; i++) {
        point_t p1 = interp(s1,i);
        point_t p2 = interp(s1,i+1);
        point_t p3 = interp(s1,i+2);
        point_t d{2*p2[0]-p1[0]-p3[0],2*p2[1]-p1[1]-p3[1],2*p2[2]-p1[2]-p3[2]};
        point_t f = d;
        addforce(s1,i,f);
        addforce(s1,i+2,f);
        for (size_t d=0; d<3; d++) f[d] = -2*f[d];
        addforce(s1,i+1,f);
      }
    }
  }
  return ret;
}
