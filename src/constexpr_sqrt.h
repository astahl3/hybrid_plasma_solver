#include "types.h"

#include <limits>

namespace Detail
{
	D_REAL constexpr sqrtNewtonRaphson(D_REAL x, D_REAL curr, D_REAL prev)
	{
		return curr == prev
			? curr
			: sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
	}
}

/*
* Constexpr version of the square root
* Return value:
*	- For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
* From https://gist.github.com/alexshtf/eb5128b3e3e143187794
*/
constexpr D_REAL ce_sqrt(D_REAL x)
{
	return x >= 0 && x < std::numeric_limits<D_REAL>::infinity()
		? Detail::sqrtNewtonRaphson(x, x, 0)
		: std::numeric_limits<D_REAL>::quiet_NaN();
}