#pragma once
#include <cstdint>
#include <iterator>

namespace glmnetpp {
namespace util {

// forward declaration
template <class IntType>
struct counting_iterator;

template <class IntType>
inline constexpr bool 
    operator==(const counting_iterator<IntType>& it1,
               const counting_iterator<IntType>& it2)
    { return it1.curr_ == it2.curr_; }

template <class IntType>
inline constexpr bool 
    operator!=(const counting_iterator<IntType>& it1,
               const counting_iterator<IntType>& it2)
    { return it1.curr_ != it2.curr_; }

template <class IntType=uint32_t>
struct counting_iterator
{
    using difference_type = int32_t;
    using value_type = IntType;
    using pointer = value_type*;
    using reference = IntType&;
    using iterator_category = std::bidirectional_iterator_tag;
    
    constexpr counting_iterator(value_type begin)
        : curr_(begin)
    {}
    
    constexpr counting_iterator& operator++() { ++curr_; return *this; }
    constexpr counting_iterator& operator--() { --curr_; return *this; }
    constexpr counting_iterator operator++(int) { auto tmp = *this; ++curr_; return tmp; }
    constexpr counting_iterator operator--(int) { auto tmp = *this; --curr_; return tmp; }
    constexpr reference operator*() { return curr_; }
    
    friend constexpr bool operator==<>(const counting_iterator&,
                                       const counting_iterator&);
    friend constexpr bool operator!=<>(const counting_iterator&,
                                       const counting_iterator&);
private:
    value_type curr_;
};

} // namespace util
} // namespace glmnetpp
