#ifdef _WIN32 

#include <windows.h>
#define os_ticks GetTickCount
#if _MSC_VER < 1900
#define snprintf(buf,len,format, ...) _snprintf_s(buf,len,len,format,__VA_ARGS__)
#endif

#if _MSC_VER < 1800
#error C++11 required. please use Visual Studio 2013 or higher
#endif

#else //GCC

#ifdef _GNUC_
#if !defined(__cplusplus)  || __cplusplus < 201103L
#error C++11 or higher required. please use g++ 4.8.1 with -std=c++14 or -std=gnu++14
#endif
#endif

#endif