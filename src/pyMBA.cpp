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

struct python_mba {
    
    py::array_t<float64> values_;
    float64 extension_u_;
    float64 extension_v_;
    float64 level_;

    std::shared_ptr<MBA> mba;

    Eigen::MatrixXd vertices_;
    Eigen::MatrixXi triangles_;

    uint32 vertex_count_;
    uint32 triangle_count_;

    python_mba(
        py::array_t<float64> values,
        float64 extension_u,
        float64 extension_v,
        uint32 level
        )
    {
        if (values.ndim() != 2)
            throw std::runtime_error("Number of dimensions must be 2");

        values_ = values;
        extension_u_ = extension_u;
        extension_v_ = extension_v;
        level_ = level;
        vertex_count_ = 0;
        triangle_count_ = 0;
    }

    void compute_fault(float64 nb_u, float64 nb_v, float64 scale)
    {
        auto r = values_.unchecked<2>(); // x must have ndim = 2; can be non-writeable

        std::vector<Eigen::Vector3d> points(r.shape(0));
        
        for (py::ssize_t i = 0; i < r.shape(0); i++)
            points[i] = Eigen::Vector3d(r(i, 0), r(i, 1), r(i, 2));

        Eigen::Vector4d pca_centroid;
        compute_3d_centroid(points, pca_centroid);
        Eigen::Matrix3d covariance;
        compute_covariance_matrix_normalized(points, pca_centroid, covariance);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(covariance, Eigen::ComputeEigenvectors);
        Eigen::Matrix3d eigen_vectors_pca = eigen_solver.eigenvectors();

        eigen_vectors_pca.col(0) = eigen_vectors_pca.col(2);
        eigen_vectors_pca.col(2) = eigen_vectors_pca.col(0).cross(eigen_vectors_pca.col(1));


        // Transform the original cloud to the origin where the principal components correspond to the axes.
        Eigen::Matrix4d projectionTransform(Eigen::Matrix4d::Identity());
        projectionTransform.block<3,3>(0,0) = eigen_vectors_pca.transpose();
        projectionTransform.block<3,1>(0,3) = -1. * (projectionTransform.block<3,3>(0,0) * pca_centroid.head<3>());

        Eigen::Matrix4d projectionTransformInv = projectionTransform.inverse();

        std::vector<double> z_save(points.size());
        std::vector<Eigen::Vector3d> proj_points(points.size());
        for(uint32 i = 0 ; i < points.size(); ++i)
        {
            Eigen::Vector4d p = projectionTransform * points[i].homogeneous();
            z_save[i] = p[2];
            proj_points[i][0] = p[0];
            proj_points[i][1] = p[1];
            proj_points[i][2] = 0;
        }

        std::vector<uint32> polygon = convex_hull(proj_points);
        auto [center, axis0, axis1, extent0, extent1, area] = min_area_rectangle_of_hull<Eigen::Vector3d>(polygon, proj_points);

        Eigen::Matrix3d rot;
        rot.row(0) = axis0;
        rot.row(1) = axis1;
        rot.row(2) = Eigen::Vector3d(0.,0.,1.);

        for(uint32 i = 0 ; i < proj_points.size(); ++i)
            proj_points[i][2] = z_save[i];


        auto x_arr = std::make_shared<std::vector<float64>>();
        auto y_arr = std::make_shared<std::vector<float64>>();
        auto z_arr = std::make_shared<std::vector<float64>>();

        for(auto &p : proj_points)
        {
            Eigen::Vector3d pp = rot * p;

            x_arr->push_back(pp[0]);
            y_arr->push_back(pp[1]);
            z_arr->push_back(pp[2]);
        }

        compute(x_arr, y_arr, z_arr);


    }

    void compute_horizon(uint32 nb_u, uint32 nb_v, float64 scale)
    {   
        auto r = values_.unchecked<2>(); // x must have ndim = 3; can be non-writeable

        auto x_arr = std::make_shared<std::vector<float64>>(r.shape(0));
        auto y_arr = std::make_shared<std::vector<float64>>(r.shape(0));
        auto z_arr = std::make_shared<std::vector<float64>>(r.shape(0));

        for(py::ssize_t i = 0; i < r.shape(0); i++)
        {
            (*x_arr)[i] = r(i,0);
            (*y_arr)[i] = r(i,1);
            (*z_arr)[i] = r(i,2);
        }

        compute(x_arr, y_arr, z_arr);

        build_surface(nb_u, nb_v, scale);
    }
    
    void build_surface(uint32 nb_u, uint32 nb_v, float64 scale)
    {
        UCBspl::SplineSurface surf = mba->getSplineSurface();

        float64 u_min = surf.umin();
	    float64 v_min = surf.vmin();
	    float64 u_max = surf.umax();
	    float64 v_max = surf.vmax();
	    float64 du = (u_max - u_min)/(float64(nb_u-1));
	    float64 dv = (v_max - v_min)/(float64(nb_v-1));

        //compute the vertices
        vertex_count_ = nb_u * nb_v;
        vertices_ = Eigen::MatrixXd(vertex_count_, 3);
        for(uint32 i = 0 ; i < nb_v ; ++i)
        {
            float64 v = v_min + i * dv;
            for(uint32 j = 0 ; j < nb_u; ++j)
            {
                float64 u = u_min + j * du;
                vertices_.row(i*nb_u+j) = Eigen::RowVector3d(u, v, surf.f(u,v));
            }
        }

        //compute the faces
        triangle_count_ = 2 * (nb_u - 1) * (nb_v - 1);
        triangles_ = Eigen::MatrixXi(triangle_count_, 3);
        uint32 idx = 0;
        for(uint32 y = 0 ; y < nb_v - 1 ; ++y)
        {
            for(uint32 x = 0 ; x < nb_u - 1; ++x)
            {
                triangles_.row(idx) = Eigen::RowVector3i(
                            (x + 0) + (y + 0) * nb_u,
                            (x + 1) + (y + 0) * nb_u,
                            (x + 0) + (y + 1) * nb_u);
                ++idx;
                triangles_.row(idx) = Eigen::RowVector3i(
                            (x + 1) + (y + 0) * nb_u,
                            (x + 1) + (y + 1) * nb_u,
                            (x + 0) + (y + 1) * nb_u);
                ++idx;
            }
        }
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
    }

    const Eigen::MatrixXd &vertices() { return vertices_; }
    const Eigen::MatrixXi &triangles() { return triangles_; }


    uint32 triangle_count()
    {
        return triangle_count_;
    }
    
    uint32 vertex_count()
    {
        return vertex_count_;
    }
	// float64 vmin()
    // {
    //     return surf.vmin();
    // }

    // float64 umax()
    // {
    //     return surf.umax();
    // }
    
    // float64 vmax()
    // {
    //     return surf.vmax();
    // }

    // float64 f(float64 u, float64 v)
    // {
    //     return surf.f(u, v);
    // }

    // TODO check memory cleaning
};

void register_mba(py::module &m) {
    std::string name = "MBA";
    std::string desc = "Multilevel B-spline Approximation";

    py::class_<python_mba>(m, name.c_str(), desc.c_str())
        .def(py::init<
                    py::array_t<float64>,
                    float64,
                    float64,
                    uint32
                    >(), py::arg("values"), py::arg("extension_u"), py::arg("extension_v"), py::arg("level")
            )
        .def("compute_horizon", &python_mba::compute_horizon, py::arg("nb_u"), py::arg("nb_v"), py::arg("scale"))
        .def("compute_fault", &python_mba::compute_fault, py::arg("nb_u"), py::arg("nb_v"), py::arg("scale"))
        .def("vertex_count", &python_mba::vertex_count)
        .def("triangle_count", &python_mba::triangle_count)
        ;
        // .def("u_max", &python_mba::umax)
        // .def("v_max", &python_mba::vmax)
        //.def("f", &python_mba::f, py::arg("u"), py::arg("v"));
        


}   

PYBIND11_MODULE(pyMBA, m) {
    register_mba(m);
}