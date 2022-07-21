#include <vector>
#include <cstdint>
#include <limits>
#include <Eigen/Dense>

using float64 = double;
using uint32 = std::uint32_t;
using int32 = std::int32_t;

enum Orientation2D
{
	ALIGNED = 0,
	RIGHT,
	LEFT
};

template <class Scalar>
inline auto almost_equal_absolute(Scalar x, Scalar y, const Scalar epsilon = std::numeric_limits<Scalar>::epsilon() ) -> typename std::enable_if<std::is_floating_point<Scalar>::value, bool>::type
{
	// static_assert(epsilon > 0);
	return std::fabs(y - x) < epsilon;
}


/**
 * return the side of point P w.r.t. the vector (Pb-Pa)
 * Tells if P is on/right/left of the line (Pa, Pb)
 * @param P the point
 * @param Pa origin point
 * @param Pb end point
 * @return the orientation
 */
template <typename VEC2a, typename VEC2b, typename VEC2c>
Orientation2D side(const Eigen::MatrixBase<VEC2a>& P, const Eigen::MatrixBase<VEC2b>& Pa, const Eigen::MatrixBase<VEC2c>& Pb)
{
	//static_assert(is_same_vector<VEC2a,VEC2b,VEC2c>::value, "parameters must have same type");
	//static_assert(is_dim_of<VEC2a, 2>::value, "The size of the vector must be equal to 2.");

	using Scalar = typename VEC2a::Scalar;

	Scalar p = (P[0] - Pa[0]) * (Pb[1] - Pa[1]) - (Pb[0] - Pa[0]) * (P[1] - Pa[1]) ;

	if (almost_equal_absolute(p, Scalar(0)))
		return Orientation2D::ALIGNED;
	else if (p > Scalar(0))
		return Orientation2D::RIGHT;
	
	return Orientation2D::LEFT;
}

// Prints convex hull of a set of n points.
// https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/
template <typename ATTR>
std::vector<uint32> convex_hull(const ATTR& positions)
{

	uint32 n = positions.size();

	// There must be at least 3 points
//	if (n < 3)
//		return;

	std::vector<uint32> hull_;

	// Find the leftmost point
	uint32 l = 0;
	for (uint32 i = 1; i < n; i++)
		if (positions[i][0] < positions[l][0])
			l = i;

	// Start from leftmost point, keep moving counterclockwise
	// until reach the start point again
	uint32 p = l, q;
	do
	{
		hull_.push_back(p);

		// Search for a point 'q' such that orientation(p, i, q) is
		// counterclockwise for all points 'i'
		q = (p + 1) % n;
		for (uint32 i = 0; i < n; i++)
			if(side(positions[p].template head<2>(), positions[i].template head<2>(), positions[q].template head<2>()) == Orientation2D::LEFT)
				q = i;

		p = q; // Set p as q for next iteration
	}
	while (p != l);

	//check if points are counterclockwise

	float64 signed_area;
	for(uint32 i = 0 ; i < hull_.size() ; i++)
	{
		uint32 elt1 = hull_[i];
		uint32 elt2 = hull_[i+1];
		if(i == hull_.size()-1) //last one
			elt2 = hull_[0];

		signed_area += positions[elt1][0] * positions[elt2][1] - positions[elt2][0] * positions[elt1][1];
	}

	signed_area /= 2.;

	if(signed_area < 0)
		std::cout << "not counterclokwise\n";

	return hull_;
}


//https://www.geometrictools.com/Documentation/MinimumAreaRectangle.pdf
template <typename VEC3>
VEC3 perp(VEC3 p)
{
	return VEC3(p[1], -p[0], p[2]);
}

template <typename VEC3, typename ATTR>
auto min_area_rectangle_of_hull(std::vector<uint32>& polygon, const ATTR& points)
{
	float64 r_area = std::numeric_limits<float64>::max();
	VEC3 center, axis[2], extent;

	for(std::size_t i0 = polygon.size() - 1, i1 = 0; i1 < polygon.size() ; i0 = i1++)
	{
		VEC3 origin = points[polygon[i0]];
		VEC3 U0 = points[polygon[i1]] - origin;
		U0.normalize();
		VEC3 U1 = -perp(U0);

		float64 min0(0), max0(0);
		float64 max1(0);

		for(std::size_t j = 0 ; j < polygon.size() ; ++j)
		{
			VEC3 D = points[polygon[j]] - origin;
			float64 dot = U0.dot(D);
			if(dot < min0)
				min0 = dot;
			else if(dot > max0)
				max0 = dot;
			dot = U1.dot(D);
			if(dot > max1)
				max1 = dot;
		}

		float64 area = (max0 - min0) * max1;

		if(area < r_area)
		{
			center = origin + ((min0 + max0) / 2.) * U0 + (max1 / 2.) * U1;
			axis[0] = U0 ;
			axis[1] = U1 ;
			extent[0] = (max0 - min0) / 2.;
			extent[1] = max1 /  2;
			r_area = area;
		}
	}

	return std::make_tuple(center, axis[0], axis[1], extent[0], extent[1], r_area);
}

template <typename ATTR, typename Scalar>
uint32 compute_3d_centroid(const ATTR& attr,
                           Eigen::Matrix<Scalar, 4, 1>& centroid)
{
    // Initialize to 0
    centroid.setZero();
    uint32 count = 0;

    for (const auto& p : attr)
    {
        centroid += p.homogeneous();
        ++count;
    }

    centroid /= static_cast<Scalar>(count);
    centroid[3] = Scalar(1.);

    return count;
}

template <typename ATTR, typename Scalar>
uint32 compute_covariance_matrix(const ATTR& attr,
                                 const Eigen::Matrix<Scalar, 4, 1> &centroid,
                                 Eigen::Matrix<Scalar, 3, 3> &covariance_matrix)
{
    // Initialize to 0
    covariance_matrix.setZero();
    uint32 point_count = 0;

    for (const auto& p : attr)
    {
        ++point_count;
        Eigen::Matrix<Scalar, 4, 1> pt = p.homogeneous() - centroid;

        covariance_matrix (1, 1) += pt[1] * pt[1];
        covariance_matrix (1, 2) += pt[1] * pt[2];

        covariance_matrix (2, 2) += pt[2] * pt[2];

        pt *= pt[0];
        covariance_matrix (0, 0) += pt[0];
        covariance_matrix (0, 1) += pt[1];
        covariance_matrix (0, 2) += pt[2];
    }

    covariance_matrix (1, 0) = covariance_matrix (0, 1);
    covariance_matrix (2, 0) = covariance_matrix (0, 2);
    covariance_matrix (2, 1) = covariance_matrix (1, 2);

    return point_count;
}

template <typename ATTR, typename Scalar>
uint32 compute_covariance_matrix_normalized(const ATTR& attr,
                                            const Eigen::Matrix<Scalar, 4, 1> &centroid,
                                            Eigen::Matrix<Scalar, 3, 3> &covariance_matrix)
{
    uint32 count = compute_covariance_matrix(attr, centroid, covariance_matrix);
    covariance_matrix /= static_cast<Scalar>(count);
    return count;
}