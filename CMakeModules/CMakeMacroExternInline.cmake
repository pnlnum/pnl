
function(check_extern_inline result)
    #
    # Checks whether the following construction is valid
    #
    # extern int foo (int);
    # extern inline int foo (int x) { return x+1;}
    # int foo (int x) { return x+1;}


    file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_extern_inline.c"
        "extern int foo (int);
        extern inline int foo (int x) { return x+1;}
        int foo (int x) { return x+1;}
        int main (int arg, char *argv[]) { foo(1); return 0;}
        "
        )

    try_compile(HAS_EXTERN_INLINE 
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_extern_inline.c")

    if(HAS_EXTERN_INLINE)
        set(PNL_HAVE_INLINE true cache internal "have 'extern inline' definition")
        set (PNL_INLINE_DECL "extern" PARENT_SCOPE)
        set (PNL_INLINE_FUNC "extern inline" PARENT_SCOPE)
        message (STATUS "Check extern inline...yes")
    else ()
        message (STATUS "Check extern inline...no")
    endif(HAS_EXTERN_INLINE)
    set (${result} ${HAS_EXTERN_INLINE} PARENT_SCOPE)

endfunction(check_extern_inline)


function(check_c99_inline result)
    #
    # Checks whether the following construction is valid
    #
    # extern inline void* foo () { foo(); return &foo; }
    # int main (int argc, char *argv[]) { return foo() != 0; }
    #


    file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_c99_inline.c"
        "extern inline void* foo () { foo(); return &foo; }
        int main (int argc, char *argv[]) { return foo() != 0; }
        "
        )

    try_compile(HAS_C99_INLINE 
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_c99_inline.c")

    if(HAS_C99_INLINE)
        set (PNL_INLINE_DECL "inline" PARENT_SCOPE)
        set (PNL_INLINE_FUNC "inline" PARENT_SCOPE)
        message (STATUS "Check C99 inline...yes")
    else ()
        message (STATUS "Check C99 inline...no")
    endif(HAS_C99_INLINE)
    set (${result} ${HAS_C99_INLINE} PARENT_SCOPE)

endfunction(check_c99_inline)


set(PNL_INLINE_DECL )
set(PNL_INLINE_FUNC )

if (MINGW)
    # If we are cross-compiling, detection of inline may be wrong. 
    # We would rather deactivating it
    message (STATUS "Check inline...not possible when cross-compiling")
    set(PNL_HAVE_INLINE "" CACHE INTERNAL "Cannot build inline functions")
else ()
    check_extern_inline (HAS_EXTERN_INLINE)
    if (NOT HAS_EXTERN_INLINE)
        check_c99_inline (HAS_C99_INLINE)
    endif (NOT HAS_EXTERN_INLINE)

    if (HAS_EXTERN_INLINE OR HAS_C99_INLINE)
        set(PNL_HAVE_INLINE 1 CACHE INTERNAL "Can build inline functions")
    else ()
        set(PNL_HAVE_INLINE "" CACHE INTERNAL "Cannot build inline functions")
    endif (HAS_EXTERN_INLINE OR HAS_C99_INLINE)
endif (MINGW)

