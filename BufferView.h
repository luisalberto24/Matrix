#pragma once 

#include <cassert>
#include "Concepts.h"

template<typename T, unsigned int N>
requires 
    nsConcepts::GreaterThanZero<N>
struct BufferView 
{
    constexpr T* begin() noexcept { return  &data[0]; }
    constexpr const T* begin() const noexcept { return  &data[0]; }
    constexpr const T* cbegin() const noexcept { return  &data[0]; }
    constexpr T* end() noexcept { return &data[N]; }
    constexpr const T* end() const noexcept { return &data[N]; }
    constexpr const T* cend() const noexcept { return &data[N]; }
    constexpr unsigned int size() const noexcept { return N; }
    const T& operator[](unsigned int i) const
    {
        assert(i < N);
        return data[i];
    }
    T& operator[](unsigned int i)
    {
        assert(i < N);
        return data[i];
    }
private:
    T data[N]{};
};