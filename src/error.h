/**
 * \file error.h
 * \brief Exception handling library
 * 
 * \details This header file contains error codes and prototypes of error-related functions. 
 * 
 * \author Ching-Yin Ng
 * \date March 2025
 */

#ifndef ERROR_H
#define ERROR_H

/* Error codes */
#define GRAV_UNKNOWN_ERROR -1
#define GRAV_SUCCESS 0
#define GRAV_FAILURE 1
#define GRAV_VALUE_ERROR 2
#define GRAV_POINTER_ERROR 3
#define GRAV_MEMORY_ERROR 4
#define GRAV_OS_ERROR 5
#define GRAV_NOT_IMPLEMENTED_ERROR 6


/* Traceback code */
#define GRAV_TRACEBACK_NOT_INITIALIZED -1
#define GRAV_TRACEBACK_SUCCESS 0
#define GRAV_TRACEBACK_MALLOC_FAILED 1
#define GRAV_TRACEBACK_TRUNCATED 2
#define GRAV_TRACEBACK_SNPRINTF_FAILED 3

typedef struct ErrorStatus
{
    int return_code;
    char *traceback;
    int traceback_code_;
} ErrorStatus;

/**
 * \brief Wrapper for raise_warning function.
 * 
 * \param error_msg Error message.
 */
#define WRAP_RAISE_WARNING(error_msg) \
    raise_warning(error_msg, __FILE__, __LINE__, __func__)

/**
 * \brief Wrapper for raise_error function.
 * 
 * \param error_status ErrorStatus struct.
 * \param error_code Error code.
 * \param error_msg Error message.
 * 
 * \return ErrorStatus struct.
 */
#define WRAP_RAISE_ERROR(error_code, error_msg) \
    raise_error(error_code, error_msg, __FILE__, __LINE__, __func__)

/**
 * \brief Wrapper for traceback function.
 * 
 * \param function_call Function call to be traced.
 * 
 * \return ErrorStatus struct.
 */
#define WRAP_TRACEBACK(function_call) \
    traceback(function_call, #function_call, __FILE__, __LINE__, __func__)

/**
 * \brief Make a error status struct with return code set to SUCCESS.
 * 
 * \return ErrorStatus struct with return code set to SUCCESS.
 */
ErrorStatus make_success_error_status(void);

/**
 * \brief Raise a warning and print to stderr.
 * 
 * \param warning_msg Warning message.
 * \param warning_file File where the warning occurs.
 * \param warning_line Line number where the warning occurs.
 * \param warning_func Function where the warning occurs.
 */
void raise_warning(
    const char *__restrict warning_msg,
    const char *__restrict warning_file,
    const int warning_line,
    const char *__restrict warning_func
);

/**
 * \brief Raise an error.
 * 
 * \param error_code Error code.
 * \param error_msg Error message.
 * \param error_file File where the error occurs.
 * \param error_line Line number where the error occurs.
 * \param error_func Function where the error occurs.
 * 
 * \return ErrorStatus struct.
 */
ErrorStatus raise_error(
    const int error_code,
    const char *__restrict error_msg,
    const char *__restrict error_file,
    const int error_line,
    const char *__restrict error_func
);

/**
 * \brief Stack traceback if error occurs.
 * 
 * \param error_status Pointer to the error status struct.
 * \param function_call_source_code Source code of the function call.
 * \param error_file File where the error occurs.
 * \param error_line Line number where the error occurs.
 * \param error_func Function where the error occurs.
 * 
 * \return ErrorStatus struct.
 */
ErrorStatus traceback(
    ErrorStatus error_status,
    const char *__restrict function_call_source_code,
    const char *__restrict error_file,
    const int error_line,
    const char *__restrict error_func
);

/**
 * \brief Free the memory allocated for the traceback string.
 */
void free_traceback(ErrorStatus *__restrict error_status);

/**
 * \brief Print the traceback string to stderr and free the memory.
 * 
 * \param error_status Pointer to the error status struct.
 */
void print_and_free_traceback(ErrorStatus *__restrict error_status);
 
#endif
