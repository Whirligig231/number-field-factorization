#include <Eigen/Core>

typedef mpz_class Z;
typedef mpq_class Q;
typedef mod ZN;
typedef poly<Z> Z_X;
typedef poly<Q> Q_X;
typedef poly<ZN> ZN_X;

template <typename T>
using mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;