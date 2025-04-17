/**
 * \file progress_bar.c
 * 
 * \brief Functions for displaying a progress bar in the terminal.
 * 
 * \author Ching-Yin Ng
 */

#include <stdbool.h>
#include <stdio.h>
#include <time.h>

#include "common.h"
#include "utils.h"
#include "progress_bar.h"

#define PROGRESS_BAR_LENGTH 40
#define MIN_UPDATE_INTERVAL_SECOND 0.1

#ifdef _WIN32
    #define BULLET " "
    #define BAR "-"
#else
    #define BULLET "\u2022"
    #define BAR "\u2501"
#endif

/* Color codes */
#define RESET "\033[0m"
#define DEEP_GREEN "\033[0;32m"
#define MOCHA_GREEN "\033[38;5;106m"
#define BRIGHT_RED "\033[38;5;197m"
#define GREY "\033[0;90m"
#define YELLOW "\033[0;33m"
#define CYAN "\033[0;36m"
#define MAGENTA "\033[0;35m"

/**
 * \brief Print the progress bar to stdout.
 * 
 * \param progress_bar_param Pointer to progress bar parameters
 * \param percent Percentage of completion
 * \param estimated_time_remaining Estimated time remaining in seconds
 * \param is_end Whether the progress bar is at the end
 */
IN_FILE void print_progress_bar(
    const ProgressBarParam *restrict progress_bar_param,
    double percent,
    const double estimated_time_remaining,
    const bool is_end
)
{
    if (percent < 0.0)
    {
        percent = 0.0;
    }
    else if (percent > 1.0)
    {
        percent = 1.0;
    }

    const int num_red_bar = (int) (percent * PROGRESS_BAR_LENGTH);
    const int num_dark_bar = PROGRESS_BAR_LENGTH - num_red_bar;

    /* Elapsed time */
    const time_t time_elapsed_time_t = grav_get_current_time() - progress_bar_param->start;

    time_t hours_elapsed = time_elapsed_time_t / 3600;
    time_t minutes_elapsed = (time_elapsed_time_t % 3600) / 60;
    time_t seconds_elapsed = time_elapsed_time_t % 60;

    if (hours_elapsed < 0)
    {
        hours_elapsed = 0;
    }
    if (minutes_elapsed < 0)
    {
        minutes_elapsed = 0;
    }
    if (seconds_elapsed < 0)
    {
        seconds_elapsed = 0;
    }

    if (hours_elapsed > 99999)
    {
        hours_elapsed = 99999;
        minutes_elapsed = 59;
        seconds_elapsed = 59;
    }

    /* Remaining time */
    char remaining_time_str[15];
    const time_t estimated_time_remaining_time_t = (time_t) estimated_time_remaining;
    time_t hours_remaining = estimated_time_remaining_time_t / 3600;
    time_t minutes_remaining = (estimated_time_remaining_time_t % 3600) / 60;
    time_t seconds_remaining = estimated_time_remaining_time_t % 60;

    if (hours_remaining > 99999)
    {
        hours_remaining = 99999;
        minutes_remaining = 59;
        seconds_remaining = 59;
    }

    if (hours_remaining < 0 || minutes_remaining < 0 || seconds_remaining < 0)
    {
        remaining_time_str[0] = '-';
        remaining_time_str[1] = '-';
        remaining_time_str[2] = ':';
        remaining_time_str[3] = '-';
        remaining_time_str[4] = '-';
        remaining_time_str[5] = ':';
        remaining_time_str[6] = '-';
        remaining_time_str[7] = '-';
        remaining_time_str[8] = '\0';
    }
    else
    {
        snprintf(
            remaining_time_str,
            15,
            "%02d:%02d:%02d",
            (int) hours_remaining,
            (int) minutes_remaining,
            (int) seconds_remaining
        );
    }

    /* Print progress bar */
    fputs("\r\033[?25l", stdout); // Start of line and hide cursor
    if (!is_end)
    {
        fputs(BRIGHT_RED, stdout);
    }
    else
    {
        fputs(MOCHA_GREEN, stdout);
    }

    for (int i = 0; i < num_red_bar; i++)
    {
        fputs(BAR, stdout);
    }
    fputs(GREY, stdout);
    for (int i = 0; i < num_dark_bar; i++)
    {
        fputs(BAR, stdout);
    }
    printf(
        "%s %3d%%%s %s %s%02d:%02d:%02d%s %s %s%s%s\033[K",
        DEEP_GREEN,
        (int) (percent * 100),
        RESET,
        BULLET,
        YELLOW,
        (int) hours_elapsed,
        (int) minutes_elapsed,
        (int) seconds_elapsed,
        RESET,
        BULLET,
        CYAN,
        remaining_time_str,
        RESET
    );

    if (is_end)
    {
        fputs("\n", stdout); // New line
        fputs("\033[?25h", stdout); // Show cursor
    }

    fflush(stdout);
}

ErrorStatus start_progress_bar(
    ProgressBarParam *restrict progress_bar_param,
    const double total
)
{
    ErrorStatus error_status;

    progress_bar_param->start = grav_get_current_time();
    progress_bar_param->time_last_update = progress_bar_param->start;
    progress_bar_param->current_progress = 0.0;
    progress_bar_param->total = total;
    if (progress_bar_param->total <= 0.0)
    {
        error_status = WRAP_RAISE_ERROR(GRAV_VALUE_ERROR, "Total must be greater than 0.");
        goto error;
    }

    progress_bar_param->last_five_progress_percent[0] = 0.0;
    progress_bar_param->diff_time_last_five_update[0] = 0.0;
    progress_bar_param->at_least_five_count = 0;

    print_progress_bar(progress_bar_param, 0.0, 0, false);

    return make_success_error_status();

error:
    return error_status;
}

IN_FILE time_t least_squares_regression_remaining_time(
    const ProgressBarParam *restrict progress_bar_param,
    const double diff_now_start
)
{
    const double target_x = 1.0;
    const double x[5] = {
        progress_bar_param->last_five_progress_percent[0],
        progress_bar_param->last_five_progress_percent[1],
        progress_bar_param->last_five_progress_percent[2],
        progress_bar_param->last_five_progress_percent[3],
        progress_bar_param->last_five_progress_percent[4]
    };
    const double y[5] = {
        progress_bar_param->diff_time_last_five_update[0],
        progress_bar_param->diff_time_last_five_update[1],
        progress_bar_param->diff_time_last_five_update[2],
        progress_bar_param->diff_time_last_five_update[3],
        progress_bar_param->diff_time_last_five_update[4]
    };

    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_x_squared = 0.0;
    double sum_xy = 0.0;

    for (int i = 0; i < 5; i++)
    {
        sum_x += x[i];
        sum_y += y[i];
        sum_x_squared += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    const double m = (5.0 * sum_xy - sum_x * sum_y) / (5.0 * sum_x_squared - sum_x * sum_x);
    const double b = (sum_y - m * sum_x) / 5.0;

    const time_t estimated_time_remaining = m * target_x + b - diff_now_start;

    return estimated_time_remaining;
}

void update_progress_bar(
    ProgressBarParam *restrict progress_bar_param,
    double current_progress,
    bool is_end
)
{
    const double current_time = grav_get_current_time();
    if (
        (current_time - progress_bar_param->time_last_update) < MIN_UPDATE_INTERVAL_SECOND
        && !is_end
    )
    {
        return;
    }
    
    progress_bar_param->time_last_update = current_time;
    progress_bar_param->current_progress = current_progress;
    double percent = progress_bar_param->current_progress / progress_bar_param->total;
    if (percent > 1.0)
    {
        percent = 1.0;
    }

    const double diff_now_start = current_time - progress_bar_param->start;

    if (progress_bar_param->at_least_five_count < 5)
    {
        progress_bar_param->last_five_progress_percent[progress_bar_param->at_least_five_count] = percent;
        progress_bar_param->diff_time_last_five_update[progress_bar_param->at_least_five_count] = diff_now_start;
        (progress_bar_param->at_least_five_count)++;
    }
    else
    {
        progress_bar_param->last_five_progress_percent[0] = progress_bar_param->last_five_progress_percent[1];
        progress_bar_param->last_five_progress_percent[1] = progress_bar_param->last_five_progress_percent[2];
        progress_bar_param->last_five_progress_percent[2] = progress_bar_param->last_five_progress_percent[3];
        progress_bar_param->last_five_progress_percent[3] = progress_bar_param->last_five_progress_percent[4];
        progress_bar_param->last_five_progress_percent[4] = percent;

        progress_bar_param->diff_time_last_five_update[0] = progress_bar_param->diff_time_last_five_update[1];
        progress_bar_param->diff_time_last_five_update[1] = progress_bar_param->diff_time_last_five_update[2];
        progress_bar_param->diff_time_last_five_update[2] = progress_bar_param->diff_time_last_five_update[3];
        progress_bar_param->diff_time_last_five_update[3] = progress_bar_param->diff_time_last_five_update[4];
        progress_bar_param->diff_time_last_five_update[4] = diff_now_start;
    }

    if (is_end)
    {
        // We still need to estimate the remaining time
        // even when the progress bar ends since it could
        // be triggered by an error rather than the actual
        // completion of the task.
        if (progress_bar_param->at_least_five_count < 5)
        {
            print_progress_bar(progress_bar_param, percent, -1.0, true);
        }
        else 
        {
            time_t estimated_time_remaining = (least_squares_regression_remaining_time(progress_bar_param, diff_now_start));
            print_progress_bar(progress_bar_param, percent, estimated_time_remaining, true);
        }
    }
    else
    {
        progress_bar_param->time_last_update = current_time;

        if (progress_bar_param->at_least_five_count < 5)
        {
            print_progress_bar(progress_bar_param, percent, -1.0, false);
        }
        else
        {
            time_t estimated_time_remaining = (least_squares_regression_remaining_time(progress_bar_param, diff_now_start));
            print_progress_bar(progress_bar_param, percent, estimated_time_remaining, false);
        }
    }
}
