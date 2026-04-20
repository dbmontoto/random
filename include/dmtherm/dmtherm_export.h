#ifndef DMTHERM_DMTHERM_EXPORT_H
#define DMTHERM_DMTHERM_EXPORT_H

// Symbol visibility / export control for the DMTherm Engine DLL.
//
// Consumers should include `dmtherm/dmtherm.h` (which includes this header).

#if defined(_WIN32) || defined(__CYGWIN__)
  #if defined(DMTHERM_ENGINE_BUILD)
    #define DMTHERM_API __declspec(dllexport)
  #elif defined(DMTHERM_ENGINE_USE_DLL)
    #define DMTHERM_API __declspec(dllimport)
  #else
    #define DMTHERM_API
  #endif
#else
  #if defined(DMTHERM_ENGINE_BUILD)
    #define DMTHERM_API __attribute__((visibility("default")))
  #else
    #define DMTHERM_API
  #endif
#endif

#endif // DMTHERM_DMTHERM_EXPORT_H

