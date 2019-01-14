/*cppimport
<%
setup_pybind11(cfg)
cfg['include_dirs'] = ['../../ext/pybind11/include']
cfg['compiler_args'] = ['-std=c++14','-march=native']

cfg['dependencies'] = ['pybind.h']
%>
*/

#define PD_MODULE_NAME py_power_diagram_2d_double
#define PD_TYPE        double
#define PD_DIM         2


#include "../../ext/power_diagram/src/PowerDiagram/Visitors/internal/ZCoords.cpp"
#include "../../ext/power_diagram/src/PowerDiagram/system/ThreadPool.cpp"
#include "../../ext/power_diagram/src/PowerDiagram/system/Assert.cpp"

#include "pybind.h"
