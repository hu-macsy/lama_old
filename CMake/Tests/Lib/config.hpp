
#ifdef DLL_LINKING

    // Microsoft Visual C++ compiler requires dllimport / dllexport

    #ifdef COMPILING_DLL

        #define DLL_IMPORTEXPORT   __declspec( dllexport )

    #else // COMPILING_DLL is defined

        #define DLL_IMPORTEXPORT   __declspec( dllimport )

    #endif //COMPILING_DLL

#else

    // ignore DLL_IMPORTEXPORT for other compilers

    #define DLL_IMPORTEXPORT 

#endif

