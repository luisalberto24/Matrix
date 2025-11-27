#pragma once

#include <memory>

struct TypeFactory
{
		template<typename T, typename ...Args>
		static std::unique_ptr<T> Create(Args&& ...args)
		{
			return std::make_unique<T>(std::forward<Args>(args)...);
		}
};