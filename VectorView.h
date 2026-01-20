#pragma once

#include <cassert>
#include <stdlib.h>
#include <concepts>
#include "Concepts.h"

enum class ViewType
{ 
	Row, 
	Column 
};

template<typename T, ViewType viewType>
class VectorView
{
	public:
		VectorView() noexcept = delete;
		VectorView(T* array, unsigned int size, unsigned int stride) noexcept :
			viewData(array),
			viewSize(size),
			viewStride(stride)
		{
		}

		VectorView& operator=(const std::initializer_list<T>& array) noexcept
		{
			assert(array.size() == viewSize);

			unsigned int r = 0;
			for (const T& value : array) { viewData[r * viewStride] = value; r++; }

			return *this;
		}

		T& operator()(unsigned int i) noexcept
		{
			assert(i < viewSize);
			return viewData[i * viewStride];
		}

		const T& operator()(unsigned int i) const noexcept
		{
			assert(i < viewSize);
			return viewData[i * viewStride];
		}

		unsigned int size() const noexcept
		{
			return viewSize;
		}

		template<typename P, ViewType Q>
		friend void Print(const VectorView<P, Q>&);

	private:

		T* viewData;
		unsigned int viewSize;
		unsigned int viewStride;
};