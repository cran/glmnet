#pragma once
#include <Eigen/Dense>
#include <string>

namespace glmnetpp {

Eigen::MatrixXd load_csv_direct(const std::string& path);
Eigen::MatrixXd load_csv(const std::string& argv0, const std::string& path);

} // namespace glmnetpp
