// from <boost/functional.hpp>

// In this header file, binder1st_void, bind1st_void, etc, and
// patially-specialized version of mem_fun_t, mem_ref_t, etc are
// defined, in order to support compilers with the `void-return bug'.

#ifndef BOOST_FUNCTIONAL_VOID_HPP
#define BOOST_FUNCTIONAL_VOID_HPP

#include <boost/functional.hpp>

namespace boost
{
    // --------------------------------------------------------------------------
    // binder1st, bind1st
    // --------------------------------------------------------------------------
    template <class Operation>
    class binder1st_void
        : public std::unary_function<typename binary_traits<Operation>::second_argument_type,
                                     typename binary_traits<Operation>::result_type>
    {       
      public:
        binder1st_void(typename binary_traits<Operation>::param_type x,
                  typename call_traits<typename binary_traits<Operation>::first_argument_type>::param_type y)
            :
            op(x), value(y)
        {}
        
        typename binary_traits<Operation>::result_type
        operator()(typename call_traits<typename binary_traits<Operation>::second_argument_type>::param_type x) const
        {
            op(value, x);
        }
        
      protected:
        typename binary_traits<Operation>::function_type op;
        typename binary_traits<Operation>::first_argument_type value;
    };

    template <class Operation>
    inline binder1st_void<Operation> bind1st_void(const Operation &op,
                                        typename call_traits<
                                                    typename binary_traits<Operation>::first_argument_type
                                        >::param_type x)
    {
        // The cast is to placate Borland C++Builder in certain circumstances.
        // I don't think it should be necessary.
        return binder1st_void<Operation>((typename binary_traits<Operation>::param_type)op, x);
    }

    template <class Operation>
    inline binder1st_void<Operation> bind1st_void(Operation &op,
                                        typename call_traits<
                                                    typename binary_traits<Operation>::first_argument_type
                                        >::param_type x)
    {
        return binder1st_void<Operation>(op, x);
    }

    // --------------------------------------------------------------------------
    // binder2nd, bind2nd
    // --------------------------------------------------------------------------
    template <class Operation>
    class binder2nd_void
        : public std::unary_function<typename binary_traits<Operation>::first_argument_type,
                                     typename binary_traits<Operation>::result_type>
    {
      public:
        binder2nd_void(typename binary_traits<Operation>::param_type x,
                  typename call_traits<typename binary_traits<Operation>::second_argument_type>::param_type y)
            :
            op(x), value(y)
        {}
        
        typename binary_traits<Operation>::result_type
        operator()(typename call_traits<typename binary_traits<Operation>::first_argument_type>::param_type x) const
        {
            op(x, value);
        }               
        
      protected:
        typename binary_traits<Operation>::function_type op;
        typename binary_traits<Operation>::second_argument_type value;
    };

    template <class Operation>
    inline binder2nd_void<Operation> bind2nd_void(const Operation &op,
                                        typename call_traits<
                                                    typename binary_traits<Operation>::second_argument_type
                                        >::param_type x)
    {
        // The cast is to placate Borland C++Builder in certain circumstances.
        // I don't think it should be necessary.
        return binder2nd_void<Operation>((typename binary_traits<Operation>::param_type)op, x);
    }

    template <class Operation>
    inline binder2nd_void<Operation> bind2nd_void(Operation &op,
                                        typename call_traits<
                                                    typename binary_traits<Operation>::second_argument_type
                                        >::param_type x)
    {
        return binder2nd_void<Operation>(op, x);
    }

    // --------------------------------------------------------------------------
    // mem_fun, etc
    // --------------------------------------------------------------------------
    template <class T>
    class mem_fun_t<void, T> : public std::unary_function<T*, void>
    {
      public:
        explicit mem_fun_t(void (T::*p)())
            :
            ptr(p)
        {}
        void operator()(T* p) const
        {
            (p->*ptr)();
        }
      private:
        void (T::*ptr)();
    };

    template <class T, class A>
    class mem_fun1_t<void, T, A> : public std::binary_function<T*, A, void>
    {
      public:   
        explicit mem_fun1_t(void (T::*p)(A))
            :
            ptr(p)
        {}
        void operator()(T* p, typename call_traits<A>::param_type x) const
        {
            (p->*ptr)(x);
        }
      private:
        void (T::*ptr)(A);
    };

    template <class T>
    class const_mem_fun_t<void, T> : public std::unary_function<const T*, void>
    {
      public:
        explicit const_mem_fun_t(void (T::*p)() const)
            :
            ptr(p)
        {}
        void operator()(const T* p) const
        {
            (p->*ptr)();
        }
      private:
        void (T::*ptr)() const;        
    };

    template <class T, class A>
    class const_mem_fun1_t<void, T, A> : public std::binary_function<const T*, A, void>
    {
      public:
        explicit const_mem_fun1_t(void (T::*p)(A) const)
            :
            ptr(p)
        {}
        void operator()(const T* p, typename call_traits<A>::param_type x) const
        {
            (p->*ptr)(x);
        }
      private:
        void (T::*ptr)(A) const;
    };
    
    // --------------------------------------------------------------------------
    // mem_fun_ref, etc
    // --------------------------------------------------------------------------
    template <class T>
    class mem_fun_ref_t<void, T> : public std::unary_function<T&, void>
    {
      public:
        explicit mem_fun_ref_t(void (T::*p)())
            :
            ptr(p)
        {}
        void operator()(T& p) const
        {
            (p.*ptr)();
        }
      private:
        void (T::*ptr)();
    };

    template <class T, class A>
    class mem_fun1_ref_t<void, T, A> : public std::binary_function<T&, A, void>
    {
      public:
        explicit mem_fun1_ref_t(void (T::*p)(A))
            :
            ptr(p)
        {}
        void operator()(T& p, typename call_traits<A>::param_type x) const
        {
            (p.*ptr)(x);
        }
      private:
        void (T::*ptr)(A);
    };
    
    template <class T>
    class const_mem_fun_ref_t<void, T> : public std::unary_function<const T&, void>
    {
      public:
        explicit const_mem_fun_ref_t(void (T::*p)() const)
            :
            ptr(p)
        {}
        
        void operator()(const T &p) const
        {
            (p.*ptr)();
        }
      private:
        void (T::*ptr)() const;
    };

    template <class T, class A>
    class const_mem_fun1_ref_t<void, T, A> : public std::binary_function<const T&, A, void>
    {
      public:
        explicit const_mem_fun1_ref_t(void (T::*p)(A) const)
            :
            ptr(p)
        {}

        void operator()(const T& p, typename call_traits<A>::param_type x) const
        {
            (p.*ptr)(x);
        }
      private:
        void (T::*ptr)(A) const;
    };
    
} // namespace boost

#endif // BOOST_FUNCTIONAL_VOID_HPP
