#pragma once
#include <type_traits>

namespace nsConcepts
{
	template <unsigned N, unsigned int M>
	concept EqualNumbers = (N == M);

	template <unsigned N, unsigned int M>
	concept GreaterThanNumber = (N > M);

	template <unsigned int N>
	concept GreaterThanZero = GreaterThanNumber<N, 0>;
}
