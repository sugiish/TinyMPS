#include "timer.h"

Timer::Timer()
{
	current_time = 0.0;
    initial_time = 0.0;
    finish_time = 0.0;

    current_time_step = 0.0;
    initial_time_step = 0.0;

    cfl_condition = 0.0;
}

Timer::~Timer()
{

}