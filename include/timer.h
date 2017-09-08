// Copyright (c) 2017 Shota SUGIHARA
// Distributed under the MIT License.
#ifndef MPS_TIMER_H_INCLUDED
#define MPS_TIMER_H_INCLUDED

namespace tiny_mps {

// Holds data on time
class Timer {
public:
    Timer() {
        initialize(0, 0, 0);
    }

    Timer(Condition & condition) {
        initialize(condition.initial_time, condition.finish_time, condition.delta_time);
    }

    Timer(double initial_time, double finish_time, double delta_time) {
        initialize(initial_time, finish_time, delta_time);
    }

    virtual ~Timer(){}

    inline void initialize(double initial_time, double finish_time, double delta_time) {
        this->current_time = initial_time;
        this->initial_time = initial_time;
        this->finish_time = finish_time;
        this->current_delta_time = delta_time;
        this->initial_delta_time = delta_time;
        this->loop_count = 0;
    }

    inline bool hasNextLoop() {
        if (current_time > finish_time) return false;
        else return true;
    }

    inline void update() {
        current_time += current_delta_time;
        loop_count++;
    }

    inline void limitCurrentDeltaTime(double max_speed, const Condition& condition) {
        if (max_speed <= 0) return;
        
        current_delta_time = initial_delta_time;
        double dt = condition.average_distance * condition.courant_number / max_speed;
        if (dt < current_delta_time) {
            current_delta_time = dt;
        }

        if (condition.viscosity_calculation == false) return;
        dt = condition.diffusion_number * condition.average_distance * condition.average_distance 
                / condition.kinematic_viscosity;
        if (dt < current_delta_time) {
            current_delta_time = dt;
        }
    }

    inline double getCurrentTime() const { return current_time; }
    inline double getInitialTime() const { return initial_time; }
    inline double getFinishTime() const { return finish_time; }
    inline double getCurrentDeltaTime() const { return current_delta_time; }
    inline double getInitialDeltaTime() const { return initial_delta_time; }
    inline void setInitialDeltaTime(double delta_time) { initial_delta_time = delta_time; }

private:
    double current_time;
    double initial_time;
    double finish_time;
    double current_delta_time;
    double initial_delta_time;
    int loop_count;
};

} // namespace tiny_mps
#endif //MPS_TIMER_H_INCLUDED