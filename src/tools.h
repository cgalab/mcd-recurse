/* CG:SHOP 2020: Minimum Convex Decomposition -- Recursor Tool
*
*  Copyright 2019, 2020 Peter Palfraader
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "geom.h"

#include <random>
#include <vector>

/** A vector that does not change its data array due to reallaction.
 *
 * This means that pointers to it are usually not invalidated due
 * to adding elements.  This means we also need to know the size in advance.
 */
template <class T>
struct FixedVector
  : private std::vector<T> {
  private:
    using Base = std::vector<T>;
    using size_type = typename Base::size_type;
    using value_type = typename Base::value_type;

    using Base::capacity;
  public:
    FixedVector() = default;
    /* No copy construction or copy assignment */
    FixedVector(const FixedVector&) = delete;
    FixedVector& operator=(const FixedVector&) = delete;
    /* Moving is fine, however. */
    FixedVector(FixedVector&&) noexcept = default;
    FixedVector& operator=(FixedVector&&) noexcept = default;

    using const_iterator = typename Base::const_iterator;

    void reserve(size_type new_size) {
      assert(size() == 0);
      Base::reserve(new_size);
    }
    void resize(size_type new_size, const value_type& val) {
      assert(size() == 0);
      Base::resize(new_size, val);
    }
    using Base::at;
    using Base::back;
    using Base::begin;
    using Base::data;
    using Base::end;
    using Base::size;
    using Base::operator[];

    void emplace_back (value_type&& val) {
      assert(size() < capacity());
      Base::emplace_back(std::forward<value_type>(val));
    }
    void push_back (const value_type& val) {
      assert(size() < capacity());
      Base::push_back(val);
    }
    void push_back (value_type&& val) {
      assert(size() < capacity());
      Base::push_back(std::forward<value_type>(val));
    }
};



/** The signum function.
 *
 * From https://stackoverflow.com/a/4609795 by user79758.
 */
template <typename T> inline constexpr
int signum(T x, std::false_type is_signed) {
    return T(0) < x;
}

template <typename T> inline constexpr
int signum(T x, std::true_type is_signed) {
    return (T(0) < x) - (x < T(0));
}

template <typename T> inline constexpr
int signum(T x) {
    return signum(x, std::is_signed<T>());
}

/** Get a random element from an iterator
 */
template<typename Iter, typename RNG>
Iter random_element(Iter begin, Iter end, RNG& rng) {
    std::uniform_int_distribution<> rnd_dist(0, std::distance(begin, end) - 1);
    std::advance(begin, rnd_dist(rng));
    return begin;
}

/** allocate memory.
 *
 * Allocates size bytes of memory, exiting the program if allocation
 * fails.
 */
inline
void *
my_malloc_c(size_t size) {
  void *result;

  result = malloc(size);

  if (result == NULL) {
    fprintf(stderr, "Out of memory on malloc()\n.");
    exit(1);
  }

  return result;
}
#define my_free_c(p) STMT_BEGIN \
  free(p); \
  p = NULL; \
  STMT_END


/** compiler hints */
#define LIKELY(condition) __builtin_expect(static_cast<bool>(condition), 1)
#define UNLIKELY(condition) __builtin_expect(static_cast<bool>(condition), 0)


/** Return a pair in sorted order */
inline
std::pair<unsigned,unsigned>
sorted_pair(unsigned u, unsigned v) {
  if (u < v) {
    return std::make_pair(u, v);
  } else {
    return std::make_pair(v, u);
  }
}
