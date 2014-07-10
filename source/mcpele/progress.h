#ifndef _MCPELE_PROGRESS_H
#define _MCPELE_PROGRESS_H

#include <ctime>
#include <utility>
#include <iostream>

namespace mcpele{


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
        stm << "---" <<  "\n";
        stm << "percentage done" <<  "\n";
        stm << curr << " %" <<  "\n";
        prev = curr;
        stm << "---" <<  "\n";
        // time elapsed
        stm << "time elapsed" <<  "\n";
        print_elapsed_time(stm);
        // estimated time to completion
        stm << "---" <<  "\n";
        stm << "estimated time to completion" <<  "\n";
        print_estimated_time(((float_t)(total_iterations - smp - 1) / (float_t)(smp + 1)) * (float_t)clock(), stm);
        // estimated total time
        stm << "---" <<  "\n";
        stm << "estimated total run time" <<  "\n";
        print_estimated_time(((float_t)total_iterations / (float_t)(smp + 1)) * (float_t)clock(), stm);
        stm << "---" <<  "\n";
        // estimated completion time in local time
        stm << "estimated completion local time" <<  "\n";
        time_t timer = time(NULL);
        timer+=(((float_t)(total_iterations - smp - 1) / (float_t)(smp + 1)) * (float_t)clock()) / CLOCKS_PER_SEC;
        stm << ctime(&timer);
        stm << "---" <<  "\n";
    }

    void print_elapsed_time(std::ostream& stm)
    {
        long_t tm = clock();
        long_t days = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
        long_t hours = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
        long_t minutes = tm/((long_t)CLOCKS_PER_SEC*(long_t)60);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60);
        long_t seconds = tm/((long_t)CLOCKS_PER_SEC);
        if (days){
            if (days>1) stm << days << " days  ";
            else stm << days << " day  ";
        }
        if (hours){
            if (hours>1) stm << hours << " hours  ";
            else stm << hours << " hour  ";
        }
        if (minutes){
            if (minutes>1) stm << minutes << " minutes  ";
            else stm << minutes << " minute  ";
        }
        if (seconds){
            if (seconds>1) stm << seconds << " seconds";
            else stm << seconds << " second";
        }
        stm <<  "\n";
    }

    void print_estimated_time(const long_t inp, std::ostream& stm)
    {
        long_t tm = inp;
        long_t days = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
        long_t hours = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
        long_t minutes = tm/((long_t)CLOCKS_PER_SEC*(long_t)60);
        tm%=((long_t)CLOCKS_PER_SEC*(long_t)60);
        long_t seconds = tm/((long_t)CLOCKS_PER_SEC);
        if (days){
            if (days>1) stm << days << " days  ";
            else stm << days << " day  ";
        }
        if (hours){
            if (hours>1) stm << hours << " hours  ";
            else stm << hours << " hour  ";
        }
        if (minutes){
            if (minutes>1) stm << minutes << " minutes  ";
            else stm << minutes << " minute  ";
        }
        if (seconds){
            if (seconds>1) stm << seconds << " seconds";
            else stm << seconds << " second";
        }
        stm <<  "\n";
    }

};


}//namespace mcpele

#endif//#ifndef _MCPELE_PROGRESS_H
