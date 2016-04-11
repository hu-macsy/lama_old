# OpenMP usage
heading3 ( "OpenMP" "USE_OPENMP" )
    found_message ( "OpenMP" "OPENMP_VERSION" "OPTIONAL" "Version ${OPENMP_VERSION}" )
    found_message ( "compile flag" "OpenMP_CXX_FLAGS" "OPTIONAL" "${OpenMP_CXX_FLAGS}" )
    found_message ( "schedule type" "SCAI_OMP_SCHEDULE" "OPTIONAL" "set to \"${SCAI_OMP_SCHEDULE}\"" )

# LAMA CUDA
heading3 ( "CUDA" "CUDA_ENABLED" )
    found_message ( "CUDA" "CUDA_FOUND" "OPTIONAL" "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
    found_message ( "Compute Capability" "CUDA_HAVE_GPU" "OPTIONAL" "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
heading3 ( "MIC" "USE_MIC" )