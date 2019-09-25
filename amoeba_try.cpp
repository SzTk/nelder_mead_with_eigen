#include <Eigen/Dense>
#include <iostream>
#define NDIM 4
#define NPOINT NDIM + 1
#define NMAX 1e+2

using namespace Eigen;

extern MatrixXd readMatrix(const char *filename);

double fun(VectorXd &p) {
  // return pQ(p-2*eta)
  MatrixXd Q = MatrixXd::Identity(NDIM, NDIM);
  VectorXi _eta(NDIM);
  VectorXd eta(NDIM);
  srand((unsigned int)time(0));
  eta = 10 * VectorXd::Random(NDIM);
  _eta = eta.cast<int>();
  eta = _eta.cast<double>();
  return p.transpose() * Q * (p - 2 * eta);
}

double amoeba_try(MatrixXd *p, VectorXd *y, VectorXd *psum, int worstidx,
                  double alpha, double (*fun)(VectorXd &)) {
  std::cout << "compute amoeba_try" << std::endl;
  double ytry;
  Eigen::VectorXd ptry(NDIM);
  Eigen::VectorXd _psum(NDIM);
  int ndim = (*p).cols();

  // compute new poit to try.
  std::cout << "alpha:" << alpha << std::endl;
  std::cout << "p:" << (*p) << std::endl;
  std::cout << "psum:" << (*psum) << std::endl;
  double alpha1 = (1.0 - alpha) / ndim;
  double alpha2 = alpha1 - alpha;
  ptry = alpha1 * (*psum) - alpha2 * (*p).row(worstidx).transpose();

  // compute the value at new point.
  ytry = (*fun)(ptry);
  std::cout << "try point is\n" << ptry << std::endl;
  std::cout << "new value is :" << ytry << std::endl;
  // if new point is better than worst, swap worst and new.
  if (ytry < (*y)(worstidx)) {
    (*y)(worstidx) = ytry;
    (*psum) += ptry - (*p).row(worstidx).transpose();
    (*p).row(worstidx) = ptry;
  }
  return ytry;
}

int main() {
  MatrixXd p(NPOINT, NDIM); // 4dimension 5 point
  VectorXd y(NPOINT), _temp(NDIM), psum(NDIM);
  double ytry, ysave;
  double rtol, ftol = 1e-7;
  double (*Fun)(VectorXd &) = fun;
  int ihi, ilo, inhi;

  // Initialize Symplex.

  p = readMatrix("test.csv");

  // clang-format on

  // compute value at each point.
  for (int i = 0; i < p.rows(); i++) {
    _temp = p.row(i);
    y(i) = (*Fun)(_temp);
  }

  psum = p.colwise().sum();
  std::cout << "p:\n" << p << std::endl;

  for (int k = 0; k < NMAX; k++) {
    // compute best,worst and next worst point.
    ilo = 0;
    ihi = y(0) > y(1) ? (inhi = 1, 0) : (inhi = 0, 1);
    for (int i; i < y.size(); i++) {
      if (y(i) <= y(ilo))
        ilo = i;
      if (y(i) > y(ilo)) {
        inhi = ihi;
        ihi = i;
      } else if (y(i) > y(inhi) && i != ihi) {
        inhi = i;
      }
    }
    std::cout << "ihi:" << ihi << "\ninhi:" << inhi << "\nilo:" << ilo
              << std::endl;
    rtol = 2.0 * fabs(y(ihi) - y(ilo)) / (fabs(y(ihi)) + fabs(y(ilo)));
    if (rtol < ftol) {
      std::cout << rtol << std::endl;
      std::cout << "optimal solution found" << std::endl;
      break;
    }
    std::cout << "compute next point, first try opposite" << std::endl;

    // First, try opposide point.
    ytry = amoeba_try(&p, &y, &psum, ihi, -1.0, (*Fun));
    if (ytry < y(ilo)) {
      // if opposide point is good, more scceed.
      std::cout << "new value is better than current best." << std::endl;
      ytry = amoeba_try(&p, &y, &psum, ihi, 2.0, (*Fun));
    } else if (ytry >= y(inhi)) {
      // if before try is all fail, try compressiong.
      std::cout << "new value is worse than inhi, try compressing" << std::endl;
      ysave = y(ihi);
      ytry = amoeba_try(&p, &y, &psum, ihi, 0.5, (*Fun));
      if (ytry >= ysave) {
        std::cout << "new value is worse than previous try. try compress to "
                     "best point"
                  << std::endl;
        for (int i = 0; i < y.size(); i++) {
          if (i != ilo) {
            p.row(i) = 0.5 * (p.row(i) + p.row(ilo));
            _temp = p.row(i);
            y(i) = (*Fun)(_temp);
          }
        }
        psum = p.colwise().sum();
      }
    }
    std::cout << "p:\n" << p << std::endl;
    std::cout << "y:\n" << y << std::endl;
  }
  std::cout << "p:\n" << p << std::endl;
  std::cout << "y:\n" << y << std::endl;
  std::cout << "result:\ny:\n" << y(ilo) << "\np:\n" << p.row(ilo) << std::endl;
}