#pragma once
#include <cstdint>
#include <iterator>

namespace glmnetpp {
namespace util {

// forward declaration
template <class IntType>
struct one_to_zero_iterator;

template <class IntType>
inline constexpr bool 
    operator==(const one_to_zero_iterator<IntType>& it1,
               const one_to_zero_iterator<IntType>& it2)
    { return it1.curr_ == it2.curr_; }

template <class IntType>
inline constexpr bool 
    operator!=(const one_to_zero_iterator<IntType>& it1,
               const one_to_zero_iterator<IntType>& it2)
    { return it1.curr_ != it2.curr_; }

template <class IntType=uint32_t>
struct one_to_zero_iterator
{
    using difference_type = int32_t;
    using value_type = IntType;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::bidirectional_iterator_tag;
    
    one_to_zero_iterator(const value_type* begin)
        : curr_(begin)
    {}
    
    one_to_zero_iterator& operator++() { ++curr_; return *this; }
    one_to_zero_iterator& operator--() { --curr_; return *this; }
    one_to_zero_iterator operator++(int) { auto tmp = *this; ++curr_; return tmp; }
    one_to_zero_iterator operator--(int) { auto tmp = *this; --curr_; return tmp; }
    value_type operator*() { return *curr_ - 1; }
    
    friend constexpr bool operator==<>(const one_to_zero_iterator&,
                                       const one_to_zero_iterator&);
    friend constexpr bool operator!=<>(const one_to_zero_iterator&,
                                       const one_to_zero_iterator&);
private:
    const value_type* curr_;
};

} // namespace util
} // namespace glmnetpp
