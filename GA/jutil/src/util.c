
#if defined (__cplusplus)
extern "C" {
#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/stat.h>

#include "util.h"

Log_Level log_level         = UTIL_DEFAULT_LOG_LEVEL;
int print_leading_spaces    = 0;

void mkdir_ifn_exists(const char * dirname)
{
    if (!file_exists(dirname))
    {
        mkdir(dirname, 0700);
    }
}

bool file_exists(const char * filename)
{
    struct stat buffer;
    return (stat (filename, &buffer) == 0);
}

bool float_approximates(
    const float actual,
    const float ref)
{
    return float_approximates_err(actual, ref, 1e-5);
}

bool float_approximates_err(
    const float actual,
    const float ref,
    const float tolerance)
{
    const float d = fabs(actual - ref);

    const bool approx = d <= tolerance;
    // if (!approx)
    //     wrn("Difference exceeds tolerance (%.10f>%.10f)!\n",
    //         d, tolerance);
    return approx;
}

void write_binary(
    const char* filename,
    const char* ptr,
    size_t write_size)
{
    FILE *outFile = fopen(filename, "wb");

    if (!outFile)
        die("Could not open output file \"%s\"\n", filename);

    while (write_size > 0)
    {
        const size_t wrote = fwrite(ptr, write_size, 1, outFile);

        if (wrote == 0)
            break;

        write_size -= (write_size * wrote);
    }

    fclose(outFile);
}


void print_arg_title(const char *title)
{
    set_leading_spaces(4);
    raw("%s\n", title);
    reset_leading_spaces();
}

void print_arg_bundles(const argument_bundle **argbv, const int n)
{
    for (int i = 0; i < n; i++)
    {
        const argument_bundle *ab = &((*argbv)[i]);
        set_leading_spaces(8);
        raw("%s, %s\n", ab->short_form, ab->long_form);
        set_leading_spaces(12);
        raw("%s\n", ab->description);
    }
    (*argbv) += n;
    reset_leading_spaces();
}

bool arg_matches(const char *arg_in, const argument_bundle ab)
{
    // anything shorter than 2 chars cannot match
    if (arg_in[0] == '\0' || arg_in[1] == '\0')
        return false;

    // the +1 for short form is to skip "-"
    // the +2 for long form is to skip "--"
    return 0 == strcmp((char *) (arg_in + 1), ab.short_form) ||
           0 == strcmp((char *) (arg_in + 2), ab.long_form);
}

void parse_args(const int argc,
                char const *argv[],
                const int argbc,
                const argument_bundle *argbv,
                void (*const failure_callback) (const char *failure_arg_in))
{
    for (int i = 1; i < argc; i++)
    {
        bool found_arg = false;

        // iterate through argument bundle to match argument
        // call callback function if matched
        for (int j = 0; j < argbc; j++)
        {
            if ((found_arg = arg_matches(argv[i], argbv[j])))
            {
                if (argbv[j].call_back != NULL)
                {
                    if (argbv[j].exp_val )
                    {
                        if (i + 1 > argc - 1)
                        {
                            err("Argument %s expects a pairing value\n",
                                argv[i]);
                            failure_callback(argv[i]);
                            break;
                        }
                        else
                        {
                            argbv[j].call_back(argv[++i]);
                        }
                    }
                    else
                    {
                        argbv[j].call_back(argv[i]);
                    }
                    break;
                }
            }
        }

        if (!found_arg)
        {
            err("Unknown argument: %s\n", argv[i]);
            failure_callback(argv[i]);
        }
    }
}

int get_leading_spaces()
{
    return print_leading_spaces;
}

void reset_leading_spaces()
{
    print_leading_spaces = 0;
}

void set_leading_spaces(int n)
{
    print_leading_spaces = n;
}

float parse_float(const char *str)
{
    char *next;
    float val = strtof(str, &next);
    if (strlen(next))
        die("Failed to parse float from \"%s\"\n", str);
    return val;
}

long parse_long(const char *str)
{
    char *next;
    long val = strtol(str, &next, 0); // base 0 allows format detection
    if (strlen(next))
        die("Failed to parse long integer from \"%s\"\n", str);

    return val;
}

const char *get_error_string()
{
    return strerror(errno);
}

Log_Level get_log_level()
{
    return log_level;
}

void set_log_level(Log_Level lvl)
{
    if (lvl < 0 || lvl >= LOG_ERROR)
    {
        err("Invalid log level: %d\n", lvl);
        err("Must be between 0 and %d\n", LOG_RAW);
        raw("   LOG_PROOF       : %d\n", LOG_PROOF);
        raw("   LOG_DEBUG       : %d\n", LOG_DEBUG);
        raw("   LOG_WARN        : %d\n", LOG_WARN);
        raw("   LOG_MESSAGE     : %d\n", LOG_MESSAGE);
        raw("   LOG_RAW         : %d\n", LOG_RAW);
        // raw("   LOG_ERROR  : %d\n", LOG_ERROR);
        // raw("   LOG_DEATH  : %d\n", LOG_DEATH);
        exit(1);
    }
    log_level = lvl;
    msg("Log level set to %s\n", Log_Level_String[log_level]);
}

void _log(const char *filename, const int line, const Log_Level lvl, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    char space_buffer[print_leading_spaces + 1];
    memset(space_buffer, ' ', print_leading_spaces);
    space_buffer[print_leading_spaces] = '\0';

    char *new_fmt;
    FILE *fd;
    switch (lvl)
    {
    case LOG_PROOF:
        fd = stdout;
        asprintf(&new_fmt, "%s[PRF] %s%s%s",
                 CLR_CYN, CLR_NRM, space_buffer, fmt);
        break;
    case LOG_DEBUG:
        fd = stdout;
        asprintf(&new_fmt, "%s[DBG] %s%s%s",
                 CLR_BLU, CLR_NRM, space_buffer, fmt);
        break;
    case LOG_WARN:
        fd = stdout;
        asprintf(&new_fmt, "%s[WRN] %s%s%s",
                 CLR_YEL, CLR_NRM, space_buffer, fmt);
        break;
    case LOG_MESSAGE:
        fd = stdout;
        asprintf(&new_fmt, "%s[MSG] %s%s%s",
                 CLR_MAG, CLR_NRM, space_buffer, fmt);
        break;
    case LOG_RAW:
        fd = stdout;
        asprintf(&new_fmt, "%s%s", space_buffer, fmt);
        break;
    case LOG_ERROR:
        fd = stderr;
        asprintf(&new_fmt, "%s[ERR] %s%s%s",
                 CLR_RED, CLR_NRM, space_buffer, fmt);
        break;
    case LOG_DEATH:
        fd = stderr;
        asprintf(&new_fmt, "\n%s[DIE %s:%d] %s%s%s",
                 CLR_RED, filename,  line, CLR_NRM, space_buffer, fmt);
        break;
    default:
        fd = stdout;
        asprintf(&new_fmt, "%s[???] %s%s%s",
                 CLR_RED, CLR_NRM, space_buffer, fmt);
    }

    vfprintf(fd, new_fmt, args);
    fflush(fd);

    va_end(args);

    if (lvl == LOG_DEATH)
    {
        if (errno != 0)
            fprintf(fd, "\nError before death: code %d (%s)\n", errno, get_error_string());
        exit(1);
    }
}

void nsleep(long ns)
{
    struct timespec sleep_time;
    sleep_time.tv_sec = ns / 1000000000;
    sleep_time.tv_nsec = ns % 1000000000;
    nanosleep(&sleep_time, NULL);
}

double get_timestamp_us()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_usec + tv.tv_sec * 1e6;
}

void test_logs()
{
    prf("Proof\n");

    dbg("Debug\n");

    wrn("Warning\n");

    msg("Message\n");

    raw("Raw\n");

    set_log_level(LOG_DEBUG);
    raw_at(LOG_DEBUG, "Raw at Debug (should show)\n");

    set_log_level(LOG_WARN);
    raw_at(LOG_DEBUG, "Raw at Debug (should not show)\n");

    raw_at(LOG_WARN, "Raw at Warning (should show)\n");

    set_log_level(LOG_MESSAGE);
    raw_at(LOG_WARN, "Raw at Warning (should NOT show)\n");

    err("Error\n");

    msg("Test done - now exit with die()...\n");
    die("Exit successful\n");
}

#if defined(_TEST_JUTIL)
int main(void)
{

    test_logs();
    return 0;
}
#endif

#if defined (__cplusplus)
}
#endif