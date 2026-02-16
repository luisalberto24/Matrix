#pragma once
template<typename T>
class Node
{
public:
	Node() noexcept = default;
	Node(T pVvalue) : value(pVvalue), next(nullptr), prev(nullptr) {}
	const T& GetData() const noexcept { return value; }
	const void SetData(T pValue) const noexcept { value = pValue; }
	operator T& () noexcept { return value; }
	operator T& () const noexcept { return value; }
protected:
	T value;
	Node* next;
	Node* prev;
	

};

