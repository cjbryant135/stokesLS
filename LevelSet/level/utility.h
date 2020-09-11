#ifndef __UTILITY_H__
#define __UTILITY_H__

template <class T>
inline T max(const T& x, const T& y) 
{return x>y ? x : y;}

template <class T>
inline T max(const T& x, const T& y, const T& z) 
{return max(max(x,y),z);}

template <class T>
inline T max(const T& w, const T& x, const T& y, const T& z) 
{return max(max(w,x),max(y,z));}

template <class T>
inline T max(const T* v, const int n)
{T ans = v[0]; for (int i=1; i<n; ++i) ans = max(ans,v[i]); return(ans);}

template <class T>
inline T min(const T& x, const T& y) 
{return x<y ? x : y;}

template <class T>
inline T min(const T& x, const T& y, const T& z) 
{return min(min(x,y),z);}

template <class T>
inline T min(const T& w, const T& x, const T& y, const T& z) 
{return min(min(w,x),min(y,z));}

template <class T>
inline T min(const T* v, const int n)
{T ans = v[0]; for (int i=1; i<n; ++i) ans = min(ans,v[i]); return(ans);}

template <class T>
inline void swap(T& x, T& y) {T t = x; x = y; y = t;}

inline double frac(const double v, const double y0, const double y1)
{return (v-y0)/(y1-y0);} 

#endif


