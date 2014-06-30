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
    const float_t oom;
    const index_t m;
    index_t curr;
    index_t prev;
public:

    progress(const index_t totalin)
    :oom(1.0/static_cast<float_t>(totalin))
    ,m(totalin)
    ,curr(0)
    ,prev(0)
	{}

    void next(const index_t idx, std::ostream& stm = std::cout){
	curr = static_cast<index_t>(static_cast<float_t>(idx)*oom*100);
	if (curr!=prev){
	    TimePercentage(idx-1);
	}
    }

    void TimePercentage(const index_t smp){
	// percentage done
	std::cout << "---" << std::endl;
	std::cout << "percentage done" << std::endl;
	std::cout << curr << " %" << std::endl;
	prev = curr;
	std::cout << "---" << std::endl;
	// time elapsed
	std::cout << "time elapsed" << std::endl;
	PrintTime();
	// estimated time to completion
	std::cout << "---" << std::endl;
	std::cout << "estimated time to completion" << std::endl;
	IntToTime(((float_t)(m-smp-1)/(float_t)(smp+1))*(float_t)clock());
	// estimated total time
	std::cout << "---" << std::endl;
	std::cout << "estimated total run time" << std::endl;
	IntToTime(((float_t)m/(float_t)(smp+1))*(float_t)clock());
	std::cout << "---" << std::endl;
	// estimated completion time in local time
	std::cout << "estimated completion local time" << std::endl;
	time_t timer = time(NULL);
	timer+=(((float_t)(m-smp-1)/(float_t)(smp+1))*(float_t)clock())/CLOCKS_PER_SEC;
	std::cout << ctime(&timer);
	std::cout << "---" << std::endl;
    }

    void PrintTime(){
	long_t tm = clock();
	long_t days = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
	long_t hours = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
	long_t minutes = tm/((long_t)CLOCKS_PER_SEC*(long_t)60);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60);
	long_t seconds = tm/((long_t)CLOCKS_PER_SEC);
	if (days){
		if (days>1) std::cout << days << " days  ";
		else std::cout << days << " day  ";
	}
	if (hours){
		if (hours>1) std::cout << hours << " hours  ";
		else std::cout << hours << " hour  ";
	}
	if (minutes){
		if (minutes>1) std::cout << minutes << " minutes  ";
		else std::cout << minutes << " minute  ";
	}
	if (seconds){
		if (seconds>1) std::cout << seconds << " seconds";
		else std::cout << seconds << " second";
	}
	std::cout << std::endl;
    }

    void IntToTime(const long_t inp){
	long_t tm = inp;
	long_t days = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60*(long_t)24);
	long_t hours = tm/((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60*(long_t)60);
	long_t minutes = tm/((long_t)CLOCKS_PER_SEC*(long_t)60);
	tm%=((long_t)CLOCKS_PER_SEC*(long_t)60);
	long_t seconds = tm/((long_t)CLOCKS_PER_SEC);
	if (days){
		if (days>1) std::cout << days << " days  ";
		else std::cout << days << " day  ";
	}
	if (hours){
		if (hours>1) std::cout << hours << " hours  ";
		else std::cout << hours << " hour  ";
	}
	if (minutes){
		if (minutes>1) std::cout << minutes << " minutes  ";
		else std::cout << minutes << " minute  ";
	}
	if (seconds){
		if (seconds>1) std::cout << seconds << " seconds";
		else std::cout << seconds << " second";
	}
	std::cout << std::endl;
    }

};


}//namespace mcpele

#endif//#ifndef _MCPELE_PROGRESS_H
