#ifndef __galileo_predictive_visitor_fusion_hpp__
#define __galileo_predictive_visitor_fusion_hpp__

#define BOOST_FUSION_INVOKE_MAX_ARITY 12

#include <boost/variant/static_visitor.hpp>
#include <boost/fusion/include/invoke.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>

namespace galileo
{
  namespace fusion
  {

    namespace bf = boost::fusion;
    using NoArg = boost::blank;

  } // namespace fusion

} // namespace galileo

namespace boost
{
  namespace fusion
  {

    // Append the element T at the front of boost fusion vector V.
    template <typename T, typename V>
    typename result_of::push_front<V const, T>::type append(T const &t, V const &v)
    {
      return push_front(v, t);
    }

    // Append the elements T1 and T2 at the front of boost fusion vector V.
    template <typename T1, typename T2, typename V>
    typename result_of::push_front<typename result_of::push_front<V const, T2>::type const, T1>::
        type
        append(T1 const &t1, T2 const &t2, V const &v)
    {
      return push_front(push_front(v, t2), t1);
    }

    // Append the elements T1, T2 and T3 at the front of boost fusion vector V.
    template <typename T1, typename T2, typename T3, typename V>
    typename result_of::push_front<
        typename result_of::push_front<typename result_of::push_front<V const, T3>::type const, T2>::
            type const,
        T1>::type
    append(T1 const &t1, T2 const &t2, T3 const &t3, V const &v)
    {
      return push_front(push_front(push_front(v, t3), t2), t1);
    }

    // Append the elements T1, T2, T3 and T4 at the front of boost fusion vector V.
    template <typename T1, typename T2, typename T3, typename T4, typename V>
    typename result_of::push_front<
        typename result_of::push_front<
            typename result_of::
                push_front<typename result_of::push_front<V const, T4>::type const, T3>::type const,
            T2>::type const,
        T1>::type
    append(T1 const &t1, T2 const &t2, T3 const &t3, T4 const &t4, V const &v)
    {
      return push_front(push_front(push_front(push_front(v, t4), t3), t2), t1);
    }

    // Append the elements T1, T2, T3, T4 and T5 at the front of boost fusion vector V.
    template <typename T1, typename T2, typename T3, typename T4, typename T5, typename V>
    typename result_of::push_front<
        typename result_of::push_front<
            typename result_of::push_front<
                typename result_of::
                    push_front<typename result_of::push_front<V const, T5>::type const, T4>::type const,
                T3>::type const,
            T2>::type const,
        T1>::type
    append(T1 const &t1, T2 const &t2, T3 const &t3, T4 const &t4, T5 const &t5, V const &v)
    {
      return push_front(push_front(push_front(push_front(push_front(v, t5), t4), t3), t2), t1);
    }

  } // namespace fusion

} // namespace boost

#endif // __galileo_predictive_visitor_fusion_hpp__