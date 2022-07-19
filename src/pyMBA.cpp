#include <vector> 
#include <memory>
#include <algorithm>
#include <cstdint>

// #include <iostream>
// #include <string>
// #include <sstream>
// #include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/functional.h>
// #include <pybind11/stl.h>

#include <MBA.h>
#include <UCButils.h>
#include <PointAccessUtils.h>

namespace py = pybind11;
using float64 = double;
using uint32 = std::uint32_t;

template<typename T>
auto inline begin(std::shared_ptr<T> ptr) -> typename T::iterator { return ptr->begin(); }

template<typename T>
auto inline end(std::shared_ptr<T> ptr) -> typename T::iterator { return ptr->end(); }

struct python_mba {
    
    std::shared_ptr<std::vector<float64>> x_arr;
    std::shared_ptr<std::vector<float64>> y_arr;
    std::shared_ptr<std::vector<float64>> z_arr;

    std::shared_ptr<MBA> mba;
    UCBspl::SplineSurface surf;

    python_mba(
        py::array_t<float64> _x_arr,
        py::array_t<float64> _y_arr,
        py::array_t<float64> _z_arr,
        float64 extension_u,
        float64 extension_v,
        uint32 level
        )
    {
        if (_x_arr.ndim() != 1 || _y_arr.ndim() != 1 || _z_arr.ndim() != 1)
            throw std::runtime_error("Number of dimensions must be one");

        if (_x_arr.size() != _y_arr.size() || _x_arr.size() != _z_arr.size() ||  _y_arr.size() != _z_arr.size())
            throw std::runtime_error("Input shapes must match");

        init(_x_arr, _y_arr, _z_arr, extension_u, extension_v, level);
    }

    python_mba(
        py::array_t<float64> _x_arr,
        py::array_t<float64> _y_arr,
        py::array_t<float64> _z_arr,
        float64 extension,
        uint32 level
        )
    {
        if (_x_arr.ndim() != 1 || _y_arr.ndim() != 1 || _z_arr.ndim() != 1)
            throw std::runtime_error("Number of dimensions must be one");

        if (_x_arr.size() != _y_arr.size() || _x_arr.size() != _z_arr.size() ||  _y_arr.size() != _z_arr.size())
            throw std::runtime_error("Input shapes must match");

        init(_x_arr, _y_arr, _z_arr, extension, extension, level);
    }

    std::vector<float64> make_available(py::array_t<float64> numpy_array)
    {
        std::vector<float64> data(numpy_array.size());
        std::copy(numpy_array.data(), numpy_array.data()+numpy_array.size(), data.begin());
        return data;
    }

    void init(
        py::array_t<float64> _x_arr,
        py::array_t<float64> _y_arr,
        py::array_t<float64> _z_arr,
        float64 extension_u,
        float64 extension_v,
        uint32 level
        )
    {
        x_arr = std::make_shared<std::vector<float64>>(make_available(_x_arr));
        y_arr = std::make_shared<std::vector<float64>>(make_available(_y_arr));
        z_arr = std::make_shared<std::vector<float64>>(make_available(_z_arr));

        mba = std::make_shared<MBA>(x_arr, y_arr, z_arr);

        float64 max_x = *std::max_element(begin(x_arr), end(x_arr));
	    float64 min_x = *std::min_element(begin(x_arr), end(x_arr));
	    float64 max_y = *std::max_element(begin(y_arr), end(y_arr));
	    float64 min_y = *std::min_element(begin(y_arr), end(y_arr));
	    
	std::cout << min_x << " " << max_x << " " << min_y << " " << max_y << std::endl;

        // following lines allow one to extend the generated surface of coef %
        float64 coef_u = extension_u/100.;
        float64 coef_v = extension_v/100.;

        float64 size_x = std::fabs(max_x - min_x);
        float64 size_y = std::fabs(max_y - min_y);
	    
	std::cout << size_x << " " << size_y << std::endl;

	    mba->setDomain(min_x-coef_u*size_x, min_y-coef_v*size_y, max_x+coef_u*size_x, max_y+coef_v*size_y);

      	// creation of the spline surface
        uint32 m0, n0;
        if (size_x >= size_y) {
            n0 = 1;
            m0 = std::floor((size_x / size_y) + 0.5);
        } else {
            m0 = 1;
            n0 = std::floor((size_x / size_y) + 0.5);
        }
        if (n0 == 0)
            n0 = 1;
        if (m0 == 0)
            m0 = 1;
        mba->MBAalg(m0, n0, level);

       surf = mba->getSplineSurface();
	std::cout << surf.umin() << " " << surf.umax() << " " << surf.vmin() << " " << surf.vmax() << std::endl;
    }

    float64 umin()
    {
        return surf.umin();
    }

	float64 vmin()
    {
        return surf.vmin();
    }

    float64 umax()
    {
        return surf.umax();
    }
    
    float64 vmax()
    {
        return surf.vmax();
    }

    float64 f(float64 u, float64 v)
    {
        return surf.f(u, v);
    }

    // TODO check memory cleaning
};

void register_mba(py::module &m) {
    std::string name = "MBA";
    std::string desc = "Multilevel B-spline Approximation";

    py::class_<python_mba>(m, name.c_str(), desc.c_str())
        .def(py::init<
                    py::array_t<float64>,
                    py::array_t<float64>,
                    py::array_t<float64>,
                    float64,
                    uint32
                    >(), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("extension"), py::arg("level")
            )
        .def("u_min", &python_mba::umin)
        .def("v_min", &python_mba::vmin)
        .def("u_max", &python_mba::umax)
        .def("v_max", &python_mba::vmax)
        .def("f", &python_mba::f, py::arg("u"), py::arg("v"));


}   

PYBIND11_MODULE(pyMBA, m) {
    register_mba(m);
}
