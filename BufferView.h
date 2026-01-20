#pragma once 

#include <cassert>
#include "Concepts.h"


template<typename T, unsigned int N>
requires 
    nsConcepts::GreaterThanZero<N>
struct BufferView 
{
    using BufferViewArray = T[N];

    BufferView() noexcept = delete;
    explicit BufferView(BufferViewArray& begin) noexcept : begPtr(&begin[0]), endPtr(&begin[N]){}
    explicit BufferView(T* begin, T* end, unsigned int size) noexcept : begPtr(begin), endPtr(end)
    { 
        assert
            (
                size <= N && 
                size > 0 && 
                begin != nullptr && 
                end != nullptr && 
                (begin + size) <= (end + 1)
            ); 
        endPtr = begin + size; 
    }
    constexpr T* begin() noexcept { return begPtr; }
    constexpr const T* begin() const noexcept { return begPtr; }
    constexpr const T* cbegin() const noexcept { return begPtr; }
    constexpr const T* end() const noexcept { return endPtr; }
    constexpr const T* cend() const noexcept { return endPtr; }
    constexpr unsigned int size() const noexcept { return N; }
    const T& operator[](unsigned int i) const
    {
        assert(i < N);
        return begPtr[i];
    }
    T& operator[](unsigned int i)
    {
        assert(i < N);
        return begPtr[i];
    }
private:
    T* begPtr = nullptr;
    T* endPtr = nullptr;
};