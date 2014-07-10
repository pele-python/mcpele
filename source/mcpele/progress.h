#ifndef _MCPELE_PROGRESS_H
#define _MCPELE_PROGRESS_H

#include <ctime>
#include <utility>
#include <iostream>

namespace mcpele{

/*
 * For a loop of total_iterations steps (assumed to take on average the same time),
 * progress keeps track of the time that elapsed, computes the likely total time
 * the whole loop will take, prints the local time at which the loop will likely
 * terminate, and prints how much time is left until loop terimnation.
 * Output is printed 100 times, for each percent of the loop that is done.
 *
 * The usage is as follows:
 *
 * const int iterations = 100000;
 * progress status(iterations);
 * size_t niter = 0;
 * while (niter < iterations) {
 *      do_something;
 *      ++niter;
 *      status.next(niter);
 * }
 *
 * Example output:
 *
 * ---
 * percentage done
 * 16 %
 * ---
 * time elapsed
 * 10 seconds
 * ---
 * estimated time to completion
 * 53 seconds
 * ---
 * estimated total run time
 * 1 minute  4 seconds
 * ---
 * estimated completion local time
 * Thu Jul 10 11:29:38 2014
 * ---
 *
 */

class progress{
public:
    typedef size_t index_t;
    typedef double float_t;
    typedef long long int long_t;
private:
    const float_t inverse_of_total_iterations;
    const index_t total_iterations;
    index_t curr;
    index_t prev;
public:

    progress(const index_t totalin)
        : inverse_of_total_iterations(1.0/static_cast<float_t>(totalin)),
        total_iterations(totalin),
        curr(0),
        prev(0)
    {}

    void next(const index_t idx, std::ostream& stm = std::cout)
    {
        curr = static_cast<index_t>(static_cast<float_t>(idx) * inverse_of_total_iterations * 100);
        if (curr != prev) {
            print_time_percentage(idx - 1, stm);
        }
    }

    void print_time_percentage(const index_t smp, std::ostream& stm)
    {
        // percentage done
        update_and_print_percentage_complete(stm);
        // time elapsed
        get_and_print_elapsed_time(stm);
        // estimated time to completion
        estimate_and_print_time_to_complete(smp, stm);
        // estimated total time
        estimate_and_print_total_time(smp, stm);
        // estimated completion time in local time
        estimate_and_print_completion_local_time(smp, stm);
    }

    void update_and_print_percentage_complete(std::ostream& stm)
    {
        stm << "---" <<  "\n";
        stm << "percentage done" <<  "\n";
        stm << curr << " %" <<  "\n";
        prev = curr;
        stm << "---" <<  "\n";
    }

    void get_and_print_elapsed_time(std::ostream& stm)
    {
        stm << "time elapsed" <<  "\n";
        print_estimated_time(clock(), stm);
    }

    void estimate_and_print_time_to_complete(const index_t smp, std::ostream& stm)
    {
        stm << "---" <<  "\n";
        stm << "estimated time to completion" <<  "\n";
        print_estimated_time(((float_t)(total_iterations - smp - 1) / (float_t)(smp + 1)) * (float_t)clock(), stm);
    }

    void estimate_and_print_total_time(const index_t smp, std::ostream& stm)
    {
        stm << "---" <<  "\n";
        stm << "estimated total run time" <<  "\n";
        print_estimated_time(((float_t)total_iterations / (float_t)(smp + 1)) * (float_t)clock(), stm);
        stm << "---" <<  "\n";
    }

    void estimate_and_print_completion_local_time(const index_t smp, std::ostream& stm)
    {
        stm << "estimated completion local time" <<  "\n";
        time_t timer = time(NULL);
        timer+=(((float_t)(total_iterations - smp - 1) / (float_t)(smp + 1)) * (float_t)clock()) / CLOCKS_PER_SEC;
        stm << ctime(&timer);
        stm << "---" <<  "\n";
    }

    void print_estimated_time(const long_t inp, std::ostream& stm)
    {
        long_t tm = inp;
        long_t days = tm / ((long_t)CLOCKS_PER_SEC * (long_t)60 * (long_t)60 * (long_t)24);
        tm %= ((long_t)CLOCKS_PER_SEC * (long_t)60 * (long_t)60 * (long_t)24);
        long_t hours = tm / ((long_t)CLOCKS_PER_SEC * (long_t)60 * (long_t)60);
        tm %= ((long_t)CLOCKS_PER_SEC * (long_t)60 * (long_t)60);
        long_t minutes = tm / ((long_t)CLOCKS_PER_SEC * (long_t)60);
        tm %= ((long_t)CLOCKS_PER_SEC * (long_t)60);
        long_t seconds = tm / ((long_t)CLOCKS_PER_SEC);
        if (days) {
            stm << days << " " << ((days > 1) ? "days" : "day") << " ";
        }
        if (hours) {
            stm << hours << " " << ((hours > 1) ? "hours" : "hour") << " ";
        }
        if (minutes) {
            stm << minutes << " " << ((minutes > 1) ? "minutes" : "minute") << " ";
        }
        if (seconds) {
            stm << seconds << " " << ((seconds > 1) ? "seconds" : "second");
        }
        stm <<  "\n";
    }

};


}//namespace mcpele

#endif//#ifndef _MCPELE_PROGRESS_H
