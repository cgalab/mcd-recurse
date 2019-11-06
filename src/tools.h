#include "geom.h"

#include <vector>

template <class T>
struct FixedVector
  : private std::vector<T> {
  private:
    using Base = std::vector<T>;
    using size_type = typename Base::size_type;
    using value_type = typename Base::value_type;

    using Base::capacity;
  public:
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
