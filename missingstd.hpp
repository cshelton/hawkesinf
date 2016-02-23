#ifndef MISSINGSTD_HPP
#define MISSINGSTD_HPP

namespace mystd {
	// g++ 4.8.1 seems to be missing this version of lower_bound
	template<typename _ForwardIterator, typename _Tp, typename CMP>
		_ForwardIterator
		lower_bound(_ForwardIterator __first, _ForwardIterator __last,
				const _Tp& __val, CMP cmp)
		{
			typedef typename std::iterator_traits<_ForwardIterator>::difference_type
				_DistanceType;

			_DistanceType __len = std::distance(__first, __last);
			while (__len > 0) {
				_DistanceType __half = __len >> 1;
				_ForwardIterator __middle = __first;
				std::advance(__middle, __half);
				if (cmp(*__middle,__val)) {
					__first = __middle;
					++__first;
					__len = __len - __half - 1;
				}
				else __len = __half;
			}
			return __first;
		}
}

#endif
