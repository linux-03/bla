#include <sstream>
#include <pybind11/pybind11.h>

#include "expression.h"
#include "vector.h"
#include "matrix.h"

using namespace ASC_bla;
namespace py = pybind11;

PYBIND11_MODULE(bla, m) {
  #include "bind_bla_obj.h"
}
