/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2015        by Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_GRAPH_DETAIL_SHARED_QUEUE_HPP
#define ALPS_GRAPH_DETAIL_SHARED_QUEUE_HPP

#include <vector>
#include <cassert>

namespace alps {
namespace graph {
namespace detail {

/**
  * shared_queue_view is a queue which shares the data storage (shared_queue_data).
  * Each view can only:
  *     -> remove elements from the top
  *     -> append elements to the back
  *     -> check if an element was queued before
  * Elements of the queue can not be changed.
  * Elements are pairs of a Key type (must be an unsigned type) and a Mapped type.
  * Each key can only be queued once.
  *
  * All copies of a shared_queue_view share the same data storage.
  * If a shared_queue_view gets out of scope its changes to the shared_queue_data object get rewinded in a simple way.
  * The class is simple, but not very robust (I don't have time for a more elaborate design right now):
  * for the rewind to work correctly the views *must* be destroyed in inverse order (like a stack) - otherwise e.g. push_back malfunctions.
  * TODO: build in a safeguard/check for this.
  */
template <typename Key, typename Mapped>
class shared_queue_data
{
  public:
    typedef Key    key_type;
    typedef Mapped mapped_type;
    typedef std::pair<key_type, mapped_type> value_type;

    explicit shared_queue_data(key_type max_element)
        : queue_(), seen_elements_(max_element+1)
    {
        BOOST_STATIC_ASSERT(( boost::is_unsigned<key_type>::value ));
        queue_.reserve(max_element+1);
    }

    void reset()
    {
        seen_elements_.reset();
        queue_.clear();
    }

    std::size_t size() const
    {
        return queue_.size();
    }

    value_type const& operator[](std::size_t i) const
    {
        return queue_[i];
    }

    void push_back(value_type const& t)
    {
        if(was_queued(t.first))
            return;
        queue_.push_back(t);
        seen_elements_.set(t.first);
    }

    bool was_queued(key_type k) const
    {
        return seen_elements_.test(k);
    }

    void mark_queued(key_type k)
    {
        seen_elements_.set(k);
    }

    void rewind(std::size_t position) // no throw
    {
        assert( position <= queue_.size() );
        // This should always be true, but we want to ensure no exceptions
        // for use in destructors
        if( position < queue_.size() )
        {
            for(typename std::vector<value_type>::const_iterator it = queue_.begin()+position; it != queue_.end(); ++it)
                seen_elements_.reset(it->first);
            queue_.resize(position);
        }
    }
  private:
    std::vector<value_type> queue_;
    boost::dynamic_bitset<> seen_elements_;
};

template <typename Key, typename Mapped>
class shared_queue_view
{
  private:
    typedef shared_queue_data<Key, Mapped> shared_storage;
  public:
    typedef typename shared_storage::key_type    key_type;
    typedef typename shared_storage::mapped_type mapped_type;
    typedef typename shared_storage::value_type  value_type;

    explicit shared_queue_view(shared_storage & queue)
    : queue_(queue), begin_(0), original_end_(queue.size())
    {
    }

    shared_queue_view(shared_queue_view const& v)
    : queue_(v.queue_), begin_(v.begin_), original_end_(queue_.size())
    {
    } 

    ~shared_queue_view()
    {
        // Rewind changes made by this view
        queue_.rewind(original_end_);
    }

    std::size_t size() const
    {
        return queue_.size() - begin_;
    }

    value_type const& front() const
    {
        assert( queue_.size() > begin_ );
        return queue_[begin_];
    }

    void pop_front()
    {
        ++begin_;
    }

    void push_back(value_type const& t)
    {
        queue_.push_back(t);
    }

    bool was_queued(key_type k) const
    {
        return queue_.was_queued(k);
    }
  private:
    shared_queue_view operator = (shared_queue_view const&);

    shared_storage &    queue_;
    std::size_t         begin_;
    std::size_t const   original_end_;
};


} // end namespace detail
} // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_DETAIL_SHARED_QUEUE_HPP
