#ifndef MPS_TIMER_H_INCLUDED
#define MPS_TIMER_H_INCLUDED

class Timer
{
public:
    double current_time;
    double initial_time;
    double finish_time;

    double current_time_step;
    double initial_time_step;

    double cfl_condition;

    Timer();
    virtual ~Timer();

};

#endif //MPS_TIMER_H_INCLUDED