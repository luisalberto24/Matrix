#pragma once

#include <memory>

struct TypeFactory
{
		template<typename T, typename ...Args>
		static std::unique_ptr<T> Create_Unique_Ptr(Args&& ...args)
		{
			return std::make_unique<T>(std::forward<Args>(args)...);
		}

		template<typename T, typename ...Args>
		static std::shared_ptr<T> Create_Shared_Ptr(Args&& ...args)
		{
			return std::make_shared<T>(std::forward<Args>(args)...);
		}

		template<typename T, typename ...Args>
		static T* Create(Args&& ...args)
		{
			return new T(std::forward<Args>(args)...);
		}
};