#include <vector> 
#include <memory>
#include <algorithm>


#include <iostream>
// #include <string>
// #include <sstream>
// #include <functional>

#include "algos.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
// #include <pybind11/functional.h>
// #include <pybind11/stl.h>

#include <MBA.h>
#include <UCButils.h>
#include <PointAccessUtils.h>

namespace py = pybind11;


template<typename T>
auto inline begin(std::shared_ptr<T> ptr) -> typename T::iterator { return ptr->begin(); }

template<typename T>
auto inline end(std::shared_ptr<T> ptr) -> typename T::iterator { return ptr->end(); }

struct VEC3 {
    float64 x;
    float64 y;
    float64 z;
};


struct python_mba {
    
    py::array_t<float64> b_x_arr;
    py::array_t<float64> b_y_arr;
    py::array_t<float64> b_z_arr;

    py::array_t<VEC3> values_;


    std::shared_ptr<MBA> mba;
    UCBspl::SplineSurface surf;

    float64 level_;
    float64 extension_u_;
    float64 extension_v_;

    // python_mba(
    //     py::array_t<float64> _x_arr,
    //     py::array_t<float64> _y_arr,
    //     py::array_t<float64> _z_arr,
    //     float64 extension_u,
    //     float64 extension_v,
    //     uint32 level
    //     )
    // {
    //     if (_x_arr.ndim() != 1 || _y_arr.ndim() != 1 || _z_arr.ndim() != 1)
    //         throw std::runtime_error("Number of dimensions must be one");

    //     if (_x_arr.size() != _y_arr.size() || _x_arr.size() != _z_arr.size() ||  _y_arr.size() != _z_arr.size())
    //         throw std::runtime_error("Input shapes must match");

    //     // init(_x_arr, _y_arr, _z_arr, extension_u, extension_v, level);
    //     b_x_arr = _x_arr;
    //     b_y_arr = _y_arr;
    //     b_z_arr = _z_arr;
    //     level_ = level;
    //     extension_u_ = extension_u;
    //     extension_v_ = extension_v;
    // }

    // python_mba(
    //     py::array_t<float64> _x_arr,
    //     py::array_t<float64> _y_arr,
    //     py::array_t<float64> _z_arr,
    //     float64 extension,
    //     uint32 level
    //     )
    // {
    //     if (_x_arr.ndim() != 1 || _y_arr.ndim() != 1 || _z_arr.ndim() != 1)
    //         throw std::runtime_error("Number of dimensions must be one");

    //     if (_x_arr.size() != _y_arr.size() || _x_arr.size() != _z_arr.size() ||  _y_arr.size() != _z_arr.size())
    //         throw std::runtime_error("Input shapes must match");

    //     //init(_x_arr, _y_arr, _z_arr, extension, extension, level);
    //     b_x_arr = _x_arr;
    //     b_y_arr = _y_arr;
    //     b_z_arr = _z_arr;
    //     level_ = level;
    //     extension_u_ = extension;
    //     extension_v_ = extension;
    // }

    python_mba(
        py::array_t<VEC3> values,
        float64 extension_u,
        float64 extension_v,
        uint32 level
        )
    {
        // if (_x_arr.ndim() != 1 || _y_arr.ndim() != 1 || _z_arr.ndim() != 1)
        //     throw std::runtime_error("Number of dimensions must be one");

        // if (_x_arr.size() != _y_arr.size() || _x_arr.size() != _z_arr.size() ||  _y_arr.size() != _z_arr.size())
        //     throw std::runtime_error("Input shapes must match");

        // init(_x_arr, _y_arr, _z_arr, extension_u, extension_v, level);
        values_ = values;
        level_ = level;
        extension_u_ = extension_u;
        extension_v_ = extension_v;
    }

    std::vector<float64> make_available(py::array_t<float64>& numpy_array)
    {
        std::vector<float64> data(numpy_array.size());
        std::copy(numpy_array.data(), numpy_array.data()+numpy_array.size(), data.begin());
        return data;
    }

    py::array do_slice(py::array a, py::int_ start, py::int_ stop) 
    {
        auto res = a[py::make_tuple(py::slice(start, stop, 1), py::slice(start, stop, 1))];
        return res;
    }


    void compute_fault()
    {
        // Eigen::MatrixXd values;

        // for (int i = 0; i < 100; i++)
        //     for (int j = 0; j < 100; j++)
        //         values(i,j,0) = b_x_arr
    }

    void compute_horizon()
    {        
        auto x_arr = std::make_shared<std::vector<float64>>(make_available(b_x_arr));
        auto y_arr = std::make_shared<std::vector<float64>>(make_available(b_y_arr));
        auto z_arr = std::make_shared<std::vector<float64>>(make_available(b_z_arr  ));

        compute(x_arr, y_arr, z_arr);
    }

    void compute(std::shared_ptr<std::vector<float64>> x_arr,
                std::shared_ptr<std::vector<float64>> y_arr,
                std::shared_ptr<std::vector<float64>> z_arr)
    {
        mba = std::make_shared<MBA>(x_arr, y_arr, z_arr);

        float64 max_x = *std::max_element(begin(x_arr), end(x_arr));
	    float64 min_x = *std::min_element(begin(x_arr), end(x_arr));
	    float64 max_y = *std::max_element(begin(y_arr), end(y_arr));
	    float64 min_y = *std::min_element(begin(y_arr), end(y_arr));

        // following lines allow one to extend the generated surface of coef %
        float64 coef_u = extension_u_/100.;
        float64 coef_v = extension_v_/100.;

        float64 size_x = std::fabs(max_x - min_x);
        float64 size_y = std::fabs(max_y - min_y);

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
        mba->MBAalg(m0, n0, level_);

        surf = mba->getSplineSurface();
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
        // .def(py::init<
        //             py::array_t<float64>,
        //             py::array_t<float64>,
        //             py::array_t<float64>,
        //             float64,
        //             uint32
        //             >(), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("extension"), py::arg("level")
        //     )
        .def(py::init<
                    py::array_t<VEC3>,
                    float64,
                    float64,
                    uint32
                    >(), py::arg("values"), py::arg("extension_u"), py::arg("extension_v"), py::arg("level")
            )
        // .def(py::init<
        //             py::array_t<float64>,
        //             py::array_t<float64>,
        //             py::array_t<float64>,
        //             float64,
        //             float64,
        //             uint32
        //             >(), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("extension_u"), py::arg("extension_v"), py::arg("level")
        //     )
        .def("compute_horizon", &python_mba::compute_horizon)
        .def("compute_fault", &python_mba::compute_fault)
        .def("u_min", &python_mba::umin)
        .def("v_min", &python_mba::vmin)
        .def("u_max", &python_mba::umax)
        .def("v_max", &python_mba::vmax)
        .def("f", &python_mba::f, py::arg("u"), py::arg("v"));


}   

PYBIND11_MODULE(pyMBA, m) {
    PYBIND11_NUMPY_DTYPE(VEC3, x, y, z);
    register_mba(m);
}