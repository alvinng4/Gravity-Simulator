/**
 * \file error.c
 * \brief Exception handling functions
 * 
 * \author Ching-Yin Ng
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "error.h"

#define RESET "\033[0m"
#define BRIGHT_RED "\033[1;31m"
#define DIM_RED "\033[5;31m"
#define YELLOW_BOLD "\033[1;33m"
#define CYAN_REGULAR "\033[0;36m"
#define PURPLE_REGULAR "\033[0;35m"
#define PURPLE_BRIGHT "\033[0;95m"
#define PURPLE_BRIGHT_BOLD "\033[1;95m"

WIN32DLL_API ErrorStatus make_success_error_status(void)
{
    ErrorStatus error_status = {
        .return_code = GRAV_SUCCESS,
        .traceback = NULL,
        .traceback_code_ = GRAV_TRACEBACK_NOT_INITIALIZED
    };
    return error_status;
}

WIN32DLL_API void raise_warning(
    const char *restrict warning_file,
    const int warning_line,
    const char *restrict warning_func,
    const char *restrict warning_msg
)
{
    fprintf(
        stderr,
        "%sWarning:%s In %s\"%s\"%s, line %s%d%s in %s%s%s:\n    %s%s%s\n",
        YELLOW_BOLD,
        RESET,
        CYAN_REGULAR,
        warning_file,
        RESET,
        CYAN_REGULAR,
        warning_line,
        RESET,
        CYAN_REGULAR,
        warning_func,
        RESET,
        PURPLE_REGULAR,
        warning_msg,
        RESET
    );
}

WIN32DLL_API ErrorStatus raise_warning_fmt(
    const char *restrict warning_file,
    const int warning_line,
    const char *restrict warning_func,
    const char *restrict format,
    ...
)
{
    /* Determine the message size */
    va_list args1;
    va_start(args1, format);
    int warning_msg_size = vsnprintf(NULL, 0, format, args1) + 1;
    va_end(args1);

    /* Allocate memory for the warning message */
    char *restrict warning_msg = malloc(warning_msg_size);
    if (!warning_msg)
    {
        return WRAP_RAISE_ERROR(
            GRAV_MEMORY_ERROR,
            "Failed to allocate memory for warning message"
        );
    }

    va_list args2;
    va_start(args2, format);
    int actual_msg_size = vsnprintf(warning_msg, warning_msg_size, format, args2);
    va_end(args2);

    if (actual_msg_size < 0)
    {
        free(warning_msg);
        return WRAP_RAISE_ERROR(
            GRAV_UNKNOWN_ERROR,
            "Failed to encode warning message"
        );
    }
    else if (actual_msg_size >= warning_msg_size)
    {
        free(warning_msg);
        return WRAP_RAISE_ERROR(
            GRAV_UNKNOWN_ERROR,
            "Warning message is truncated"
        );
    }

    raise_warning(
        warning_file,
        warning_line,
        warning_func,
        warning_msg
    );

    /* Free memory */
    free(warning_msg);

    return make_success_error_status();
}

WIN32DLL_API ErrorStatus raise_error(
    const char *restrict error_file,
    const int error_line,
    const char *restrict error_func,
    const int error_code,
    const char *restrict error_msg
)
{
    ErrorStatus error_status = {
        .traceback = NULL,
        .traceback_code_ = GRAV_TRACEBACK_NOT_INITIALIZED
    };

    char *error_type;
    switch (error_code)
    {
        case GRAV_FAILURE:
            error_status.return_code = error_code;
            error_type = "Failure";
            break;
        case GRAV_VALUE_ERROR:
            error_status.return_code = error_code;
            error_type = "ValueError";
            break;
        case GRAV_POINTER_ERROR:
            error_status.return_code = error_code;
            error_type = "PointerError";
            break;
        case GRAV_MEMORY_ERROR:
            error_status.return_code = error_code;
            error_type = "MemoryError";
            break;
        case GRAV_OS_ERROR:
            error_status.return_code = error_code;
            error_type = "OSError";
            break;
        case GRAV_NOT_IMPLEMENTED_ERROR:
            error_status.return_code = error_code;
            error_type = "NotImplementedError";
            break;
        default:
            error_status.return_code = GRAV_UNKNOWN_ERROR;
            error_type = "UnknownError";
            break;
    }

    const int traceback_size = (
        strlen(error_file)
        + strlen(error_func)
        + strlen(error_msg)
        + strlen(error_type)
        + 3 * strlen(CYAN_REGULAR)
        + strlen(PURPLE_BRIGHT_BOLD)
        + strlen(PURPLE_REGULAR)
        + 5 * strlen(RESET)
        + strlen("    File \"\", line  in \n: \n")
        + snprintf(NULL, 0, "%d", error_line)   // Number of digits in error_line
        + 1  // Null terminator
    );
    error_status.traceback = malloc(traceback_size * sizeof(char));
    if (!error_status.traceback)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_MALLOC_FAILED;
        free(error_status.traceback);
        error_status.traceback = NULL;

        goto err_memory_alloc;
    }

    int actual_traceback_size = snprintf(
        error_status.traceback,
        traceback_size,
        "    File %s\"%s\"%s, line %s%d%s in %s%s%s\n%s%s%s: %s%s%s\n",
        CYAN_REGULAR,
        error_file,
        RESET,
        CYAN_REGULAR,
        error_line,
        RESET,
        CYAN_REGULAR,
        error_func,
        RESET,
        PURPLE_BRIGHT_BOLD,
        error_type,
        RESET,
        PURPLE_REGULAR,
        error_msg,
        RESET
    );

    if (actual_traceback_size < 0)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SNPRINTF_FAILED;
        free(error_status.traceback);
        error_status.traceback = NULL;
    }
    else if (actual_traceback_size >= traceback_size)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_TRUNCATED;
    }
    else
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SUCCESS;
    }

    return error_status;

err_memory_alloc:
    return error_status;
}

WIN32DLL_API ErrorStatus raise_error_fmt(
    const char *restrict error_file,
    const int error_line,
    const char *restrict error_func,
    const int error_code,
    const char *restrict format,
    ...
)
{
    ErrorStatus error_status = {
        .traceback = NULL,
        .traceback_code_ = GRAV_TRACEBACK_NOT_INITIALIZED
    };

    char *error_type;
    switch (error_code)
    {
        case GRAV_FAILURE:
            error_status.return_code = error_code;
            error_type = "Failure";
            break;
        case GRAV_VALUE_ERROR:
            error_status.return_code = error_code;
            error_type = "ValueError";
            break;
        case GRAV_POINTER_ERROR:
            error_status.return_code = error_code;
            error_type = "PointerError";
            break;
        case GRAV_MEMORY_ERROR:
            error_status.return_code = error_code;
            error_type = "MemoryError";
            break;
        case GRAV_OS_ERROR:
            error_status.return_code = error_code;
            error_type = "OSError";
            break;
        case GRAV_NOT_IMPLEMENTED_ERROR:
            error_status.return_code = error_code;
            error_type = "NotImplementedError";
            break;
        default:
            error_status.return_code = GRAV_UNKNOWN_ERROR;
            error_type = "UnknownError";
            break;
    }

    /* Determine the message size */
    va_list args1;
    va_start(args1, format);
    const int formatted_string_size = vsnprintf(NULL, 0, format, args1) + 1;
    va_end(args1);

    const int error_msg_size = (
        strlen(error_file)
        + strlen(error_func)
        + (formatted_string_size - 1)
        + strlen(error_type)
        + 3 * strlen(CYAN_REGULAR)
        + strlen(PURPLE_BRIGHT_BOLD)
        + strlen(PURPLE_REGULAR)
        + 5 * strlen(RESET)
        + strlen("    File \"\", line  in \n: \n")
        + snprintf(NULL, 0, "%d", error_line)   // Number of digits in error_line
        + 1  // Null terminator
    );

    /* Allocate memory for the error message */
    char *restrict formatted_string = malloc(formatted_string_size);
    error_status.traceback = malloc(error_msg_size * sizeof(char));
    if (!error_status.traceback || !formatted_string)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_MALLOC_FAILED;
        error_status.traceback = NULL;
        goto err_malloc;
    }

    /* Format the error message */
    va_list args2;
    va_start(args2, format);
    int actual_msg_size = vsnprintf(formatted_string, formatted_string_size, format, args2);
    va_end(args2);

    if (actual_msg_size < 0)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SNPRINTF_FAILED;
        error_status.traceback = NULL;
        goto err_fmt;
    }
    else if (actual_msg_size >= formatted_string_size)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_TRUNCATED;
        error_status.traceback = NULL;
        goto err_fmt;
    }

    /* Format the error message */
    int actual_error_msg_size = snprintf(
        error_status.traceback,
        error_msg_size,
        "    File %s\"%s\"%s, line %s%d%s in %s%s%s\n%s%s%s: %s%s%s\n",
        CYAN_REGULAR,
        error_file,
        RESET,
        CYAN_REGULAR,
        error_line,
        RESET,
        CYAN_REGULAR,
        error_func,
        RESET,
        PURPLE_BRIGHT_BOLD,
        error_type,
        RESET,
        PURPLE_REGULAR,
        formatted_string,
        RESET
    );

    if (actual_error_msg_size < 0)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SNPRINTF_FAILED;
        error_status.traceback = NULL;
        goto err_fmt;
    }
    else if (actual_error_msg_size >= error_msg_size)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_TRUNCATED;
    }
    else
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SUCCESS;
    }

    free(formatted_string);

    return error_status;

err_fmt:
err_malloc:
    free(formatted_string);
    free(error_status.traceback);
    return error_status;
}

WIN32DLL_API ErrorStatus traceback(
    ErrorStatus error_status,
    const char *restrict function_call_source_code,
    const char *restrict error_file,
    const int error_line,
    const char *restrict error_func
)
{
    if (
        error_status.return_code == GRAV_SUCCESS
        || error_status.traceback_code_ != GRAV_TRACEBACK_SUCCESS
    )
    {
        return error_status;
    }

    const int traceback_size = (
        strlen(error_file)
        + strlen(error_func)
        + strlen(function_call_source_code)
        + 3 * strlen(CYAN_REGULAR)
        + strlen(PURPLE_REGULAR)
        + 4 * strlen(RESET)
        + strlen("    File \"\", line  in \n: \n        \n")
        + snprintf(NULL, 0, "%d", error_line)   // Number of digits in error_line
        + strlen(error_status.traceback)
        + 1  // Null terminator
    );

    char *new_traceback = malloc(traceback_size * sizeof(char));
    if (!new_traceback)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_MALLOC_FAILED;
        error_status.traceback = NULL;
        free(new_traceback);
        free(error_status.traceback);
        return error_status;
    }

    int actual_traceback_size = snprintf(
        new_traceback,
        traceback_size,
        "    File %s\"%s\"%s, line %s%d%s in %s%s%s\n        %s%s%s\n%s",
        CYAN_REGULAR,
        error_file,
        RESET,
        CYAN_REGULAR,
        error_line,
        RESET,
        CYAN_REGULAR,
        error_func,
        RESET,
        PURPLE_REGULAR,
        function_call_source_code,
        RESET,
        error_status.traceback
    );

    if (actual_traceback_size < 0)
    {
        error_status.traceback_code_ = GRAV_TRACEBACK_SNPRINTF_FAILED;
        free(new_traceback);
        free(error_status.traceback);
        error_status.traceback = NULL;
    }
    else if (actual_traceback_size >= traceback_size)
    {
        free(error_status.traceback);
        error_status.traceback = new_traceback;
        error_status.traceback_code_ = GRAV_TRACEBACK_TRUNCATED;
    }
    else
    {
        free(error_status.traceback);
        error_status.traceback = new_traceback;
        error_status.traceback_code_ = GRAV_TRACEBACK_SUCCESS;
    }

    return error_status;
}

WIN32DLL_API void free_traceback(ErrorStatus *restrict error_status)
{
    if (error_status->traceback)
    {
        free(error_status->traceback);
        error_status->traceback = NULL;
    }
}

WIN32DLL_API void print_and_free_traceback(ErrorStatus *restrict error_status)
{
    fprintf(stderr, "\n%sTraceback%s %s(most recent call last):%s\n", BRIGHT_RED, RESET, DIM_RED, RESET);
    switch (error_status->traceback_code_)
    {
        case GRAV_TRACEBACK_NOT_INITIALIZED:
            fputs("    Something went wrong. Traceback not initialized.\n", stderr);
            break;
        case GRAV_TRACEBACK_SUCCESS:
            fputs(error_status->traceback, stderr);
            free(error_status->traceback);
            error_status->traceback = NULL;
            break;
        case GRAV_TRACEBACK_MALLOC_FAILED:
            fputs("    Something went wrong. Failed to allocate memory for traceback.\n", stderr);
            break;
        case GRAV_TRACEBACK_TRUNCATED:
            fputs(error_status->traceback, stderr);
            fputs("\n    Something went wrong. Traceback are truncated.\n", stderr);
            free(error_status->traceback);
            error_status->traceback = NULL;
            break;
        case GRAV_TRACEBACK_SNPRINTF_FAILED:
            fputs("    Something went wrong. Failed to write to traceback.\n", stderr);
            break;
    }
}
